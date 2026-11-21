"""M9: the ECU firmware on the Vortex-V SoC RTL (Verilator), behind the same interface as
hil.bridge.HostEcu, plus a lockstep comparator.

  Rv64Ecu       ecu/targets/rv64/build/ecu_rv64.elf running on hil/build/libvortex_soc.so. Every
                call becomes a mailbox command (ecu/targets/rv64/hil_mbox.h): poke the arguments
                and the command word into SoC RAM, clock the core until the firmware acknowledges,
                peek the results. Plant time stands still while the core computes (lockstep);
                the cycles each command took are recorded for the budget report.
  RefEcu        the host build of the same core with the same portable libm (bit-exact reference).
  LockstepEcu   forwards every call to a primary and to shadows and compares events and slow-loop
                outputs bit for bit, call by call.

Vortex-V is consumed read-only from $VORTEX_V_ROOT. All results are MODELLED."""
from __future__ import annotations

import ctypes as C
import os
import re
import shutil
import struct
import subprocess
from functools import lru_cache
from pathlib import Path

import yaml

from hil import bridge
from redline.spec import ROOT

VORTEX_V_ROOT = Path(os.environ.get("VORTEX_V_ROOT", "/home/user/vortex-v"))
HIL_DIR = ROOT / "hil"
RV64_DIR = ROOT / "ecu" / "targets" / "rv64"
MODEL_LIB = HIL_DIR / "build" / "libvortex_soc.so"
ELF = RV64_DIR / "build" / "ecu_rv64.elf"
IMAGE = RV64_DIR / "build" / "ecu_rv64.bin"
REF_LIB = RV64_DIR / "build" / "libecu_host_ref.so"
MBOX_H = RV64_DIR / "hil_mbox.h"
COSIM_PATH = HIL_DIR / "cosim.yaml"
NM = "riscv64-unknown-elf-nm"
MBOX_SYMBOL = "g_mbox"
RAM_BASE_SYMBOL = "RAM_BASE"
U32 = 4
U64 = 8


@lru_cache(maxsize=None)
def cosim_spec() -> dict:
    with open(COSIM_PATH) as f:
        return yaml.safe_load(f)


def toolchain_missing() -> list:
    return [t for t in ("riscv64-unknown-elf-gcc", "verilator", "make") if shutil.which(t) is None]


def build(target: str = "all"):
    """make -C hil model firmware ref (each is a no-op when up to date)."""
    env = dict(os.environ, VORTEX_V_ROOT=str(VORTEX_V_ROOT))
    targets = ["model", "firmware", "ref"] if target == "all" else [target]
    subprocess.run(["make", "-s", "-C", str(HIL_DIR), *targets], check=True, env=env)


@lru_cache(maxsize=None)
def mbox_enums() -> dict:
    """Every enumerator of hil_mbox.h (HIL_CMD_*, HIL_L_*) and its #define integers, parsed from C."""
    text = re.sub(r"/\*.*?\*/", "", MBOX_H.read_text(), flags=re.S)
    out = {}
    for body in re.findall(r"enum\s+\w+\s*\{(.*?)\}", text, flags=re.S):
        nxt = 0
        for item in filter(None, (s.strip() for s in body.split(","))):
            name, _, val = (p.strip() for p in item.partition("="))
            nxt = int(val, 0) if val else nxt
            out[name] = nxt
            nxt += 1
    for name, val in re.findall(r"#define\s+(HIL_\w+)\s+(0x[0-9A-Fa-f]+|\d+)u?\b", text):
        out[name] = int(val, 0)
    return out


def elf_symbols(elf: Path = ELF) -> dict:
    out = subprocess.run([NM, str(elf)], check=True, capture_output=True, text=True).stdout
    return {parts[2]: int(parts[0], 16) for parts in (ln.split() for ln in out.splitlines()) if len(parts) == 3}


class RefEcu(bridge.HostEcu):
    """Host build of the unchanged core with ecu_libm.c: the bit-exact reference for Rv64Ecu."""

    def __init__(self):
        build("ref")
        super().__init__(REF_LIB)


class Rv64Ecu:
    """The firmware on the Vortex-V SoC model. Same methods as HostEcu."""

    def __init__(self, do_build: bool = True):
        if do_build:
            build()
        self.cfg = cosim_spec()
        self.E = mbox_enums()
        lib = C.CDLL(str(MODEL_LIB))
        self.lib = lib
        vp, u64 = C.c_void_p, C.c_uint64
        lib.vx_new.restype = vp
        lib.vx_delete.argtypes = [vp]
        for f in (lib.vx_ram_bytes, lib.vx_cycles):
            f.argtypes, f.restype = [vp], u64
        for f in (lib.vx_done, lib.vx_pass, lib.vx_leds, lib.vx_reset):
            f.argtypes = [vp]
        lib.vx_write.argtypes = [vp, u64, C.c_char_p, u64]
        lib.vx_read.argtypes = [vp, u64, C.c_char_p, u64]
        lib.vx_run_until.argtypes = [vp, u64, C.c_uint32, u64]
        lib.vx_run_until.restype = C.c_int64
        self.h = lib.vx_new()
        sym = elf_symbols()
        self.mbox = sym[MBOX_SYMBOL] - sym[RAM_BASE_SYMBOL]
        image = IMAGE.read_bytes()
        if lib.vx_write(self.h, 0, image, len(image)) != 0:
            raise RuntimeError("firmware image larger than the SoC RAM")
        lib.vx_reset(self.h)
        # boot: crt0 clears .bss, main() fills the layout table and then sets the magic word
        if lib.vx_run_until(self.h, self.mbox, self.E["HIL_MBOX_MAGIC"], int(self.cfg["boot_max_cycles"])) < 0:
            raise RuntimeError("firmware did not boot (no mailbox magic)")
        self.boot_cycles = int(lib.vx_cycles(self.h))
        n_l = self.E["HIL_L_N"]
        self.L = struct.unpack(f"<{n_l}I", self._rd(U32, n_l * U32))
        self.cycles: dict = {}          # command name -> list of (handler cycles, post-to-ack cycles)
        self._events: list = []
        self._max_cmd = int(self.cfg["command_max_cycles"])
        self._cmd_names = {v: k[len("HIL_CMD_"):].lower() for k, v in self.E.items() if k.startswith("HIL_CMD_")}
        self.consts = bridge.Consts(self._table("HIL_CMD_CONSTS", "HIL_L_CONSTS", self.E["HIL_NCONSTS"]))
        Inputs, Outputs = bridge._make_structs(self.consts.ncyl)
        if C.sizeof(Inputs) != self._l("HIL_L_IN_SIZE") or C.sizeof(Outputs) != self._l("HIL_L_OUT_SIZE"):
            raise RuntimeError("ctypes struct layout does not match the RV64 build of ecu.h")
        self.inp, self.out = Inputs(), Outputs()
        self.knock_n_in = self._l("HIL_L_KNOCK_N_IN")
        self.reset()

    def __del__(self):
        if getattr(self, "h", None):
            self.lib.vx_delete(self.h)
            self.h = None

    # ---- mailbox primitives ----
    def _l(self, name: str) -> int:
        return self.L[self.E[name]]

    def _rd(self, off: int, n: int) -> bytes:
        buf = C.create_string_buffer(n)
        if self.lib.vx_read(self.h, self.mbox + off, buf, n) != 0:
            raise RuntimeError("mailbox read outside RAM")
        return buf.raw

    def _wr(self, off: int, data: bytes):
        if self.lib.vx_write(self.h, self.mbox + off, data, len(data)) != 0:
            raise RuntimeError("mailbox write outside RAM")

    def _command(self, name: str, *args: int) -> int:
        """Post a command, clock the core to the acknowledge, collect emitted events."""
        cmd = self.E[name]
        if args:
            self._wr(self._l("HIL_L_ARG"), struct.pack(f"<{len(args)}I", *(a & 0xFFFFFFFF for a in args)))
        self._wr(self._l("HIL_L_CMD"), struct.pack("<I", cmd))
        total = self.lib.vx_run_until(self.h, self.mbox + self._l("HIL_L_CMD"), self.E["HIL_CMD_NONE"], self._max_cmd)
        if total < 0 or self.lib.vx_done(self.h):
            raise RuntimeError(f"firmware did not acknowledge {name}")
        handler = struct.unpack("<Q", self._rd(self._l("HIL_L_CYCLES"), U64))[0]
        self.cycles.setdefault(self._cmd_names[cmd], []).append((handler, int(total)))
        n = struct.unpack("<I", self._rd(self._l("HIL_L_N_EVENTS"), U32))[0]
        if n:
            size = self._l("HIL_L_EVENT_SIZE")
            raw = self._rd(self._l("HIL_L_EVENTS"), n * size)
            ot, oa = self._l("HIL_L_EVENT_TICK"), self._l("HIL_L_EVENT_ARG")
            oy, oc = self._l("HIL_L_EVENT_TYPE"), self._l("HIL_L_EVENT_CYL")
            for i in range(n):
                e = raw[i * size:(i + 1) * size]
                self._events.append((struct.unpack_from("<I", e, ot)[0], e[oy], e[oc],
                                     struct.unpack_from("<I", e, oa)[0]))
        return handler

    def _table(self, cmd: str, field: str, n: int):
        self._command(cmd)
        return struct.unpack(f"<{n}I", self._rd(self._l(field), n * U32))

    # ---- HostEcu interface ----
    def reset(self):
        self._command("HIL_CMD_RESET")
        self._events.clear()
        self.cycles.clear()

    def crank(self, tick: int):
        self._command("HIL_CMD_CRANK", tick)

    def cam(self, tick: int):
        self._command("HIL_CMD_CAM", tick)

    def knock(self, cyl: int, counts: int):
        self._command("HIL_CMD_KNOCK", cyl, counts)

    def sched_set(self, cyl: int, adv_u: int, pw_ticks: int):
        self._command("HIL_CMD_SCHED_SET", cyl, adv_u, pw_ticks)

    def drain(self):
        ev, self._events = self._events, []
        return ev

    def slow(self, now_tick, map_pa, iat_k, clt_k, tps1, tps2, pedal, lam, vbatt):
        i = self.inp
        i.now_tick = now_tick & 0xFFFFFFFF
        i.map_pa, i.iat_k, i.clt_k, i.tps1_deg, i.tps2_deg = map_pa, iat_k, clt_k, tps1, tps2
        i.pedal, i.lambda_, i.vbatt = pedal, lam, vbatt
        self._wr(self._l("HIL_L_IN"), bytes(i))
        self._command("HIL_CMD_SLOW")
        C.memmove(C.byref(self.out), self._rd(self._l("HIL_L_OUT"), C.sizeof(self.out)), C.sizeof(self.out))
        return self.out

    def status(self) -> dict:
        vals = self._table("HIL_CMD_STATUS", "HIL_L_STATUS", self.E["HIL_NSTATUS"])
        return dict(zip(bridge.STATUS_KEYS, (int(v) for v in vals)))

    # ---- M9 extras ----
    def knock_infer(self, x_q) -> tuple:
        """Run the int8 knock MLP on the Xmlacc unit. Returns (logit, handler cycles)."""
        data = bytes(int(v) & 0xFF for v in x_q)
        if len(data) != self.knock_n_in:
            raise ValueError(f"knock window must be {self.knock_n_in} samples")
        self._wr(self._l("HIL_L_KNOCK_X"), data)
        cyc = self._command("HIL_CMD_KNOCK_INFER")
        return struct.unpack("<i", self._rd(self._l("HIL_L_KNOCK_LOGIT"), U32))[0], cyc

    def exit(self) -> bool:
        """Ask the firmware to return from main (crt0 stores PASS to tohost). True if the SoC passed."""
        self._wr(self._l("HIL_L_CMD"), struct.pack("<I", self.E["HIL_CMD_EXIT"]))
        self.lib.vx_run_until(self.h, self.mbox + self._l("HIL_L_CMD"), self.E["HIL_CMD_NONE"], self._max_cmd)
        return bool(self.lib.vx_done(self.h)) and bool(self.lib.vx_pass(self.h))

    def total_cycles(self) -> int:
        return int(self.lib.vx_cycles(self.h))


class LockstepEcu:
    """Drive `primary` (its results go to the plant) and `shadows` with the same call stream and
    compare bit for bit: every drained event (tick, type, cyl, arg) and the raw bytes of every
    slow-loop output struct. Differences are counted per shadow, the first one is kept."""

    def __init__(self, primary, shadows: dict):
        self.primary = primary
        self.shadows = shadows
        self.consts = primary.consts
        self.n_events = 0
        self.n_slow = 0
        self.diff = {name: {"events": 0, "slow": 0, "first": None} for name in shadows}
        self.inputs: list = []          # recorded call stream, for replay

    def _all(self):
        return [self.primary, *self.shadows.values()]

    def _note(self, name, kind, detail):
        d = self.diff[name]
        d[kind] += 1
        if d["first"] is None:
            d["first"] = detail

    def reset(self):
        for e in self._all():
            e.reset()

    def crank(self, tick):
        self.inputs.append(("crank", tick))
        for e in self._all():
            e.crank(tick)

    def cam(self, tick):
        self.inputs.append(("cam", tick))
        for e in self._all():
            e.cam(tick)

    def knock(self, cyl, counts):
        self.inputs.append(("knock", cyl, counts))
        for e in self._all():
            e.knock(cyl, counts)

    def drain(self):
        ev = self.primary.drain()
        self.n_events += len(ev)
        for name, s in self.shadows.items():
            other = s.drain()
            if other != ev:
                n = sum(a != b for a, b in zip(ev, other)) + abs(len(ev) - len(other))
                self.diff[name]["events"] += n - 1
                self._note(name, "events", {"call": len(self.inputs), "primary": ev, "shadow": other})
        return ev

    def slow(self, *args):
        self.inputs.append(("slow", *args))
        out = self.primary.slow(*args)
        self.n_slow += 1
        ref = bytes(out)
        for name, s in self.shadows.items():
            o = s.slow(*args)
            if bytes(o) != ref:
                fields = [f for f, _ in type(out)._fields_ if bytes(_as_bytes(getattr(out, f))) != bytes(_as_bytes(getattr(o, f)))]
                self._note(name, "slow", {"call": len(self.inputs), "fields": fields})
        return out

    def status(self):
        st = self.primary.status()
        for name, s in self.shadows.items():
            other = s.status()
            if other != st:
                self._note(name, "events", {"status": st, "shadow": other})
        return st


def _as_bytes(v):
    return bytes(v) if isinstance(v, C.Array) else struct.pack("<f", v)
