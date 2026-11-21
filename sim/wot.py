"""Full-load (wide-open-throttle) speed sweep on the M4 1D gas-dynamics engine model.

    python -m sim.wot --sweep 1000:7000:250

writes out/sim/wot.csv and out/sim/wot.png: torque, power, volumetric efficiency and BSFC against speed.
Brake = net indicated (720 deg) - FMEP from the Chen-Flynn correlation (gasdyn/friction.py).
Every operating point starts from the same cold state, so the result does not depend on --jobs.
All numbers are MODELLED: there is no physical engine (CLAUDE.md rule 5).
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from gasdyn.engine import EngineModel
from gasdyn.spec import gasdyn_spec
from redline.spec import ROOT, load

OUT_DIR = ROOT / "out" / "sim"
# chart tokens shared with thermo/run.py
SURFACE, INK, INK_2, SERIES_1, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6", "#e4e3df"
COLUMNS = ["rpm", "torque_nm", "power_kw", "bmep_bar", "imep_net_bar", "fmep_bar", "ve_trapped", "ve_delivered",
           "residual_fraction", "bsfc_g_per_kwh", "isfc_g_per_kwh", "eta_ind_net", "eta_brake", "air_kg_per_h",
           "fuel_kg_per_h", "p_max_bar", "T_ivc_K", "cycles_run", "converged", "imep_delta", "floor_hits"]

_MODEL: EngineModel | None = None
_MODEL_KEY = None


def sweep_speeds(start, stop, step) -> list[float]:
    """Inclusive speed list [rpm]."""
    n = int(np.floor((stop - start) / step + 1e-9)) + 1
    return [float(start + i * step) for i in range(n)]


def _point(args):
    rpm, overrides = args
    global _MODEL, _MODEL_KEY
    key = repr(sorted((overrides or {}).items()))
    if _MODEL is None or _MODEL_KEY != key:
        _MODEL, _MODEL_KEY = EngineModel(overrides=overrides), key
    return _MODEL.run(rpm).metrics


def run_sweep(rpms, overrides: dict | None = None, jobs: int = 1) -> list[dict]:
    """Converged WOT operating point (metrics dict) for each speed, cold-started independently."""
    tasks = [(float(r), overrides) for r in rpms]
    if jobs <= 1 or len(tasks) == 1:
        return [_point(t) for t in tasks]
    _warm_kernels()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(_point, tasks))


def _warm_kernels():
    """Compile (or load from the numba cache) the kernels in the parent before workers fork."""
    from gasdyn.network import Network
    net = Network()
    net.add_pipe(0.1, 0.01, 1.4, 287.0, n_cells=4)
    net.finalize()
    net.advance(1.0e-6)


def peaks(results: list[dict]) -> dict:
    """Peak torque and power of a sweep and where they occur."""
    it = max(range(len(results)), key=lambda i: results[i]["torque_nm"])
    ip = max(range(len(results)), key=lambda i: results[i]["power_kw"])
    ve = [r["ve_trapped"] for r in results]
    return {"torque_nm": results[it]["torque_nm"], "torque_rpm": results[it]["rpm"],
            "power_kw": results[ip]["power_kw"], "power_rpm": results[ip]["rpm"],
            "ve_min": min(ve), "ve_max": max(ve), "ve_peak_rpm": results[int(np.argmax(ve))]["rpm"],
            "bsfc_min": min(r["bsfc_g_per_kwh"] for r in results),
            "all_converged": all(r["converged"] for r in results)}


def write_csv(results, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for r in results:
            w.writerow([r[c] if isinstance(r[c], (bool, int, np.bool_)) else f"{r[c]:.6g}" for c in COLUMNS])


def plot(results, path, targets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rpm = np.array([r["rpm"] for r in results])
    panels = [("torque_nm", "Brake torque [Nm]", (targets["torque_rpm"], targets["torque_nm"])),
              ("power_kw", "Brake power [kW]", (targets["power_rpm"], targets["power_kw"])),
              ("ve_trapped", "Volumetric efficiency, trapped air / ambient density [-]", None),
              ("bsfc_g_per_kwh", "BSFC [g/kWh]", None)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), facecolor=SURFACE, sharex=True)
    for ax, (key, label, tgt) in zip(axes.flat, panels):
        y = np.array([r[key] for r in results])
        ax.set_facecolor(SURFACE)
        ax.plot(rpm, y, color=SERIES_1, lw=2, marker="o", ms=3.5)
        k = int(np.argmin(y) if key == "bsfc_g_per_kwh" else np.argmax(y))
        ax.annotate(f"{y[k]:.3g} at {rpm[k]:.0f} rpm", (rpm[k], y[k]), textcoords="offset points",
                    xytext=(0, 8 if key != "bsfc_g_per_kwh" else -14), ha="center", color=INK, fontsize=9)
        if tgt:
            ax.plot([tgt[0]], [tgt[1]], marker="D", ms=7, mfc="none", mec=INK_2, mew=1.5, ls="none")
            ax.annotate("spec target", tgt, textcoords="offset points", xytext=(8, -12), color=INK_2, fontsize=9)
        ax.set_title(label, color=INK, loc="left", fontsize=10)
        ax.grid(True, color=GRID, lw=0.6)
        ax.margins(y=0.15)
        for s in ax.spines.values():
            s.set_color(GRID)
        ax.tick_params(colors=INK_2)
    for ax in axes[1]:
        ax.set_xlabel("Engine speed [rpm]", color=INK_2)
    fig.suptitle("Redline I4-2.0 NA: full-load curve from the 1D gas-dynamics model (MODELLED, no physical engine)",
                 color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    g = gasdyn_spec()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sweep", default=None, help="start:stop:step in rpm (default: spec/gasdyn.yaml wot.sweep_rpm)")
    ap.add_argument("--coarse", action="store_true", help="use the coarse sweep of spec/gasdyn.yaml (as the tests do)")
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1), help="worker processes")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override a spec/gasdyn.yaml value for a design study, e.g. intake.runner.length=0.45")
    ap.add_argument("--out", default=str(OUT_DIR), help="output directory")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args(argv)
    if a.sweep:
        sw = [float(v) for v in a.sweep.split(":")]
    else:
        sw = g["wot"]["coarse_sweep_rpm" if a.coarse else "sweep_rpm"]
    overrides = {}
    for kv in a.set:
        k, v = kv.split("=", 1)
        try:
            overrides[k] = float(v)
        except ValueError:
            overrides[k] = v
    results = run_sweep(sweep_speeds(*sw), overrides or None, a.jobs)
    out = Path(a.out)
    write_csv(results, out / "wot.csv")
    targets = load()["targets"]
    if not a.no_plot:
        plot(results, out / "wot.png", targets)
    print("MODELLED full-load curve (no physical engine)" + (f", overrides {overrides}" if overrides else ""))
    print(f"{'rpm':>6} {'T[Nm]':>7} {'P[kW]':>7} {'BMEP':>6} {'IMEP':>6} {'FMEP':>5} {'VE':>6} {'xr':>6} {'BSFC':>6} {'cyc':>4}")
    for r in results:
        print(f"{r['rpm']:6.0f} {r['torque_nm']:7.1f} {r['power_kw']:7.1f} {r['bmep_bar']:6.2f} {r['imep_net_bar']:6.2f} "
              f"{r['fmep_bar']:5.2f} {r['ve_trapped']:6.3f} {r['residual_fraction']:6.3f} {r['bsfc_g_per_kwh']:6.1f} "
              f"{r['cycles_run']:4d}{'' if r['converged'] else '  NOT CONVERGED'}")
    pk = peaks(results)
    print(f"peak torque {pk['torque_nm']:.1f} Nm at {pk['torque_rpm']:.0f} rpm (target {targets['torque_nm']} at "
          f"{targets['torque_rpm']}); peak power {pk['power_kw']:.1f} kW at {pk['power_rpm']:.0f} rpm (target "
          f"{targets['power_kw']} at {targets['power_rpm']}); VE {pk['ve_min']:.3f}..{pk['ve_max']:.3f}; "
          f"min BSFC {pk['bsfc_min']:.1f} g/kWh")
    print(f"wrote {out / 'wot.csv'}" + ("" if a.no_plot else f" and {out / 'wot.png'}"))
    return results


if __name__ == "__main__":
    main()
