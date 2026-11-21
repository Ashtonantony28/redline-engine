"""Closed-loop HIL runner: `python -m hil.run --target host --cycle tipin`.

Cycles (hil/cycles.yaml): idle, tipin, tipout, sweep. Prints a summary and writes
out/hil/<cycle>.csv (and out/hil/<cycle>.png with --plot). Target `host` runs the C ECU core as a
shared library against the Python MVEM. Target `rv64` (M9) runs the same core as firmware on the
Vortex-V SoC RTL in Verilator (hil/rv64.py): every ECU call is simulated cycle by cycle, about a
minute of wall time per simulated second, so shorten the cycle with --duration; the bit-exactness
and cycle-budget report is `make -C hil cosim` (hil/cosim.py).
All numbers are MODELLED: there is no physical engine.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys

from hil import runner
from redline.spec import ROOT

OUT_DIR = ROOT / "out" / "hil"
# chart tokens (dataviz reference palette, as thermo/run.py)
SURFACE, INK, INK_2, SERIES = "#fcfcfb", "#0b0b0b", "#52514e", ["#2a78d6", "#d6602a", "#3a9a5b"]


def write_csv(res, path):
    keys = runner.LOG_KEYS
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(len(res.log["t"])):
            w.writerow([f"{res.log[k][i]:.6g}" for k in keys])


def plot(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    L = res.log
    panels = [("rpm", "engine speed [rpm]"), ("lambda_sensor", "lambda (sensor)"), ("map", "MAP [Pa]"),
              ("thr", "throttle [deg]"), ("spark0", "spark cyl 1 [deg BTDC]"), ("torque_brake", "brake torque [Nm]")]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 11), sharex=True, facecolor=SURFACE)
    for ax, (k, lab) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        ax.plot(L["t"], L[k], color=SERIES[0], lw=1.2)
        if k == "thr":
            ax.plot(L["t"], L["thr_target"], color=SERIES[1], lw=1.0, ls="--", label="target")
            ax.legend(frameon=False, fontsize=8)
        ax.set_ylabel(lab, color=INK_2, fontsize=8)
        ax.tick_params(colors=INK_2, labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].set_xlabel("time [s]", color=INK_2)
    fig.suptitle(f"Redline M8 HIL, cycle {res.scenario.name} (MODELLED)", color=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None) -> int:
    cycles = [k for k in runner.load_cycles() if k != "defaults"]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--target", default="host", choices=["host", "rv64"])
    ap.add_argument("--cycle", default="tipin", choices=cycles)
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--duration", type=float, default=None, help="truncate the cycle to this many seconds")
    a = ap.parse_args(argv)
    over = {} if a.duration is None else {"duration": a.duration, "windows": {}}   # metric windows need the full cycle
    ecu = None
    if a.target == "rv64":
        from hil import rv64
        missing = rv64.toolchain_missing()
        if missing:
            print(f"target rv64 needs {missing} (run $VORTEX_V_ROOT/scripts/setup-toolchain.sh)", file=sys.stderr)
            return 2
        ecu = rv64.Rv64Ecu()
    res = runner.run(runner.Scenario.from_cycle(a.cycle, **over), ecu)
    s = runner.summary(res)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(res, OUT_DIR / f"{a.cycle}.csv")
    if a.plot:
        plot(res, OUT_DIR / f"{a.cycle}.png")
    print(f"Redline HIL, target {a.target}, cycle {a.cycle} (all values MODELLED)")
    print(json.dumps(s, indent=2, default=float))
    print(f"wrote {OUT_DIR / (a.cycle + '.csv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
