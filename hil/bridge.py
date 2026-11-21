"""ctypes bridge to the host build of the ECU core (ecu/build/libecu_host.so).

The library is (re)built with `make -C ecu lib` when missing or older than any ECU source. The
input/output structs mirror ecu/include/ecu/ecu.h; their sizes are checked against the library."""
from __future__ import annotations

import ctypes as C
import subprocess

from redline.spec import ROOT

ECU_DIR = ROOT / "ecu"
LIB = ECU_DIR / "build" / "libecu_host.so"
NSTATUS = 16
NCONSTS = 11
MAX_DRAIN = 64
STATUS_KEYS = ["sync", "tooth", "half", "period", "rev_ticks", "n_sync", "n_loss", "n_noise", "n_dropout",
               "n_cam_err", "cam_angle_u", "rev_cut", "n_cut", "n_knock", "overflow", "limp"]


def _stale() -> bool:
    if not LIB.exists():
        return True
    t = LIB.stat().st_mtime
    srcs = [*ECU_DIR.glob("src/*.[ch]"), *ECU_DIR.glob("include/ecu/*.h"), *ECU_DIR.glob("host/*.c"),
            ECU_DIR / "Makefile"]
    return any(p.stat().st_mtime > t for p in srcs)


def build():
    if _stale():
        subprocess.run(["make", "-s", "-C", str(ECU_DIR), "lib"], check=True)


class Consts:
    def __init__(self, arr):
        (self.timer_hz, self.angle_u_per_deg, self.ncyl, self.teeth_total, self.teeth_missing,
         self.tooth0_u, self.cam_edge_u, self.knock_threshold, self.ev_inj, self.ev_dwell,
         self.ev_spark) = [int(x) for x in arr]


def _make_structs(ncyl: int):
    class Inputs(C.Structure):
        _fields_ = [("now_tick", C.c_uint32), ("map_pa", C.c_float), ("iat_k", C.c_float), ("clt_k", C.c_float),
                    ("tps1_deg", C.c_float), ("tps2_deg", C.c_float), ("pedal", C.c_float), ("lambda_", C.c_float),
                    ("vbatt", C.c_float)]

    arr = C.c_float * ncyl

    class Outputs(C.Structure):
        _fields_ = [("etc_voltage", C.c_float), ("rpm", C.c_float), ("map_used", C.c_float),
                    ("throttle_target_deg", C.c_float), ("spark_deg", arr), ("knock_retard_deg", arr),
                    ("pw_s", arr), ("film_kg", arr), ("stft", C.c_float), ("ltft", C.c_float),
                    ("idle_air_deg", C.c_float), ("closed_loop", C.c_float), ("idle_active", C.c_float),
                    ("limp", C.c_float), ("faults", C.c_float), ("sync", C.c_float), ("rev_cut", C.c_float)]
    return Inputs, Outputs


class HostEcu:
    """The C ECU core with the host HAL (one static instance per process)."""

    def __init__(self, lib_path=None):
        """lib_path: another build of the same core + host HAL (M9: the portable-libm reference,
        ecu/targets/rv64/build/libecu_host_ref.so); default is the host library, rebuilt if stale."""
        if lib_path is None:
            build()
        lib = C.CDLL(str(lib_path or LIB))
        self.lib = lib
        u32 = C.c_uint32
        lib.hil_crank.argtypes = [u32]
        lib.hil_cam.argtypes = [u32]
        lib.hil_knock.argtypes = [u32, u32]
        lib.hil_sched_set.argtypes = [u32, C.c_int32, u32]
        lib.hil_sizeof_inputs.restype = u32
        lib.hil_sizeof_outputs.restype = u32
        consts = (u32 * NCONSTS)()
        lib.hil_consts(consts)
        self.consts = Consts(consts)
        Inputs, Outputs = _make_structs(self.consts.ncyl)
        if C.sizeof(Inputs) != lib.hil_sizeof_inputs() or C.sizeof(Outputs) != lib.hil_sizeof_outputs():
            raise RuntimeError("ctypes struct layout does not match ecu.h")
        self.inp = Inputs()
        self.out = Outputs()
        lib.hil_slow.argtypes = [C.POINTER(Inputs), C.POINTER(Outputs)]
        self._ticks = (u32 * MAX_DRAIN)()
        self._args = (u32 * MAX_DRAIN)()
        self._types = (C.c_uint8 * MAX_DRAIN)()
        self._cyls = (C.c_uint8 * MAX_DRAIN)()
        lib.hil_drain.argtypes = [C.POINTER(u32), C.POINTER(u32), C.POINTER(C.c_uint8), C.POINTER(C.c_uint8), C.c_int]
        lib.hil_drain.restype = C.c_int
        self._status = (u32 * NSTATUS)()
        self.reset()

    def reset(self):
        self.lib.hil_reset()

    def crank(self, tick: int):
        self.lib.hil_crank(tick & 0xFFFFFFFF)

    def cam(self, tick: int):
        self.lib.hil_cam(tick & 0xFFFFFFFF)

    def knock(self, cyl: int, counts: int):
        self.lib.hil_knock(cyl, counts)

    def sched_set(self, cyl: int, adv_u: int, pw_ticks: int):
        self.lib.hil_sched_set(cyl, adv_u, pw_ticks)

    def drain(self):
        n = self.lib.hil_drain(self._ticks, self._args, self._types, self._cyls, MAX_DRAIN)
        return [(self._ticks[i], self._types[i], self._cyls[i], self._args[i]) for i in range(n)]

    def slow(self, now_tick, map_pa, iat_k, clt_k, tps1, tps2, pedal, lam, vbatt):
        i = self.inp
        i.now_tick = now_tick & 0xFFFFFFFF
        i.map_pa, i.iat_k, i.clt_k, i.tps1_deg, i.tps2_deg = map_pa, iat_k, clt_k, tps1, tps2
        i.pedal, i.lambda_, i.vbatt = pedal, lam, vbatt
        self.lib.hil_slow(C.byref(self.inp), C.byref(self.out))
        return self.out

    def status(self) -> dict:
        self.lib.hil_status(self._status)
        return {k: int(self._status[i]) for i, k in enumerate(STATUS_KEYS)}
