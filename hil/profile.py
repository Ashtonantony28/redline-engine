"""Retire-PC profile of the ECU firmware on the Vortex-V RTL: `python -m hil.profile [window]`.

Every core cycle is charged to the next instruction that retires (hil/cosim/vortex_sim.cpp), and
the PCs are folded onto the ELF's function symbols. Used in M9.1 to find where the crank-tooth
handler and the soft-float slow loop spend their cycles. MODELLED."""
from __future__ import annotations

import bisect
import ctypes as C
import subprocess
import sys

from hil import cosim, rv64


def functions(elf) -> list:
    out = subprocess.run([rv64.NM, "-n", "--defined-only", str(elf)], check=True, capture_output=True, text=True).stdout
    syms = [(int(a, 16), n) for a, t, n in (ln.split() for ln in out.splitlines() if len(ln.split()) == 3)
            if t in "tTwW"]
    return sorted(syms)


class Profiler:
    def __init__(self, ecu: rv64.Rv64Ecu):
        self.ecu = ecu
        lib = ecu.lib
        lib.vx_prof.argtypes = [C.c_void_p, C.c_int]
        lib.vx_prof_dump.argtypes = [C.c_void_p, C.POINTER(C.c_uint64), C.POINTER(C.c_uint64), C.c_uint64]
        lib.vx_prof_dump.restype = C.c_uint64

    def start(self):
        self.ecu.lib.vx_prof(self.ecu.h, 1)

    def stop(self) -> dict:
        """{function: cycles}, most expensive first."""
        lib, h = self.ecu.lib, self.ecu.h
        n = lib.vx_prof_dump(h, None, None, 0)
        pcs, cnt = (C.c_uint64 * n)(), (C.c_uint64 * n)()
        lib.vx_prof_dump(h, pcs, cnt, n)
        lib.vx_prof(h, 0)
        syms = functions(self.ecu.elf)
        addrs = [a for a, _ in syms]
        acc: dict = {}
        for pc, c in zip(pcs, cnt):
            i = bisect.bisect_right(addrs, pc) - 1
            name = syms[i][1] if i >= 0 else hex(pc)
            acc[name] = acc.get(name, 0) + c
        return dict(sorted(acc.items(), key=lambda kv: -kv[1]))


def main(argv=None) -> int:
    window = (argv or sys.argv[1:] or ["redline"])[0]
    ecu = rv64.Rv64Ecu()
    prof = Profiler(ecu)
    prof.start()
    r = cosim.run_window(window, ecu, with_glibc=False)
    table = prof.stop()
    total = sum(table.values())
    print(f"window {window}: {total} cycles profiled, {r['slow_steps']} slow-loop steps, "
          f"{len(r['cycles']['crank'])} crank teeth (frozen-time command mode). MODELLED.")
    for name, c in list(table.items())[:25]:
        print(f"  {name:24s} {c:10d}  {100.0 * c / total:5.1f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main())
