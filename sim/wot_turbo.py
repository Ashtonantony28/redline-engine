"""Full-load curve and tip-in turbo lag of the turbocharged variant Redline I4-2.0T (M10 stretch).

    python -m sim.wot_turbo --sweep 1500:6500:500          # WOT sweep + 2000 rpm tip-in
    python -m sim.wot_turbo --sweep 1500:6500:500 --no-tipin

Each speed is a cold-started steady point of turbo.system.TurboEngine: M4 1D gas dynamics with the
synthetic compressor / turbine maps, intercooler, wastegate PI boost control (Python side) and spark
knock-limited by thermo.knock on the converged trapped state. The tip-in holds the engine at
spec/turbo.yaml tipin.rpm, steps the throttle from part load to WOT with the real rotor inertia, and reports
the time to 90 % of the MAP rise to the boost target, with spark held at the knock-limited value of the
steady WOT point at that speed.

Writes out/sim/turbo/wot_turbo.csv, tipin.csv, wot_turbo.png, tipin.png, compressor_map.png.
All numbers are MODELLED: no physical engine, SYNTHETIC turbo maps (CLAUDE.md rule 5).
"""
from __future__ import annotations

import argparse
import csv
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from redline.spec import ROOT
from sim.wot import _warm_kernels, sweep_speeds
from turbo.control import rpm_to_rad
from turbo.spec import turbo_spec

OUT_DIR = ROOT / "out" / "sim" / "turbo"
SURFACE, INK, INK_2, SERIES_1, SERIES_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6", "#d6602a", "#e4e3df"
COLUMNS = ["rpm", "torque_nm", "power_kw", "bmep_bar", "bsfc_g_per_kwh", "map_bar", "target_bar", "boost_error",
           "p_charge_bar", "turbo_krpm", "wg_opening", "pr_c", "mdot_c", "mdot_c_corr", "eta_c", "m_surge",
           "surge_margin", "surge_margin_inst_min", "surge_time", "P_C_kw", "P_T_kw", "T2_K", "T_ic_K", "T03_K",
           "p03_bar", "pr_t", "eta_t", "bsr", "mdot_t", "mdot_wg", "ve_trapped", "ve_plenum", "residual_fraction",
           "p_max_bar", "spark_deg", "mbt_deg", "klsa_deg", "klsa_status", "knock_limited", "knock_margin",
           "ca50_deg", "shaft_balance", "cycles_run", "converged", "floor_hits"]
KNOCK_THRESHOLD = 1.0      # Livengood-Wu integral at autoignition (spec/combustion_ext.yaml knock.integral_threshold)

_TE = None


def _engine():
    global _TE
    if _TE is None:
        from turbo.system import TurboEngine
        _TE = TurboEngine()
    return _TE


def _point(rpm):
    return _engine().steady_point(float(rpm))


def run_sweep(rpms, jobs=1):
    if jobs <= 1 or len(rpms) == 1:
        return [_point(r) for r in rpms]
    _warm_kernels()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(_point, rpms))


def run_tipin(spark_deg):
    return _engine().tipin(spark_deg)


def fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{float(v):.6g}"


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([fmt(r[c]) for c in COLUMNS])


def plateau(rows, band):
    """Operating points whose torque lies inside the design-intent band."""
    lo, hi = band
    return [r for r in rows if lo <= r["torque_nm"] <= hi]


def _style(ax, title):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, loc="left", fontsize=10)
    ax.grid(True, color=GRID, lw=0.6)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=INK_2)


def plot_wot(rows, path, ts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rpm = np.array([r["rpm"] for r in rows])
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), facecolor=SURFACE, sharex=True)
    tb, pb = ts["variant"]["targets"]["torque_plateau_nm"], ts["variant"]["targets"]["power_kw"]
    panels = [("torque_nm", "Brake torque [Nm]", tb), ("power_kw", "Brake power [kW]", pb),
              ("map_bar", "Intake plenum pressure (MAP) [bar abs]", None), ("bsfc_g_per_kwh", "BSFC [g/kWh]", None)]
    for ax, (k, label, band) in zip(axes.flat, panels):
        _style(ax, label)
        y = np.array([r[k] for r in rows])
        ax.plot(rpm, y, color=SERIES_1, lw=2, marker="o", ms=3.5, label="modelled")
        if band:
            ax.axhspan(band[0], band[1], color=GRID, alpha=0.6, lw=0, label="design-intent band")
        if k == "map_bar":
            ax.plot(rpm, [r["target_bar"] for r in rows], color=INK_2, ls="--", lw=1.2, label="boost target")
        ax.legend(fontsize=8, frameon=False)
        ax.margins(y=0.15)
    for ax in axes[1]:
        ax.set_xlabel("Engine speed [rpm]", color=INK_2)
    fig.suptitle("Redline I4-2.0T: knock-limited full-load curve, 1D gas dynamics + synthetic turbo maps "
                 "(MODELLED, no physical engine)", color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def plot_compressor(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    te = _engine()
    cm = te.cmap
    p1, T1 = te.eng.p_amb, te.eng.T_amb
    fig, ax = plt.subplots(figsize=(8, 6), facecolor=SURFACE)
    _style(ax, "SYNTHETIC compressor map (turbo/maps.py) with the WOT operating line (MODELLED)")
    surge_m, surge_pr = [], []
    for krpm in range(40, 181, 20):
        w = rpm_to_rad(krpm * 1.0e3)
        m, pr, eta = cm.speed_line(w, p1, T1)
        ax.plot(m, pr, color=INK_2, lw=0.8)
        ax.annotate(f"{krpm}k", (m[0], pr[0]), textcoords="offset points", xytext=(-14, 2), fontsize=7, color=INK_2)
        surge_m.append(m[0])
        surge_pr.append(pr[0])
    ax.plot(surge_m, surge_pr, color=SERIES_2, lw=1.5, label="surge line (dPR/dm = 0)")
    mg = np.linspace(1.0e-3, 0.3, 120)
    wg = np.linspace(rpm_to_rad(30.0e3), rpm_to_rad(190.0e3), 120)
    M, W = np.meshgrid(mg, wg)
    PR = np.full(M.shape, np.nan)
    ETA = np.full(M.shape, np.nan)
    for i in range(M.shape[0]):
        ms, _, mch = cm.limits(W[i, 0], p1, T1)
        for j in range(M.shape[1]):
            if ms <= M[i, j] < mch:
                pt = cm.line(W[i, j], p1, T1, M[i, j])
                PR[i, j], ETA[i, j] = pt["pr"], pt["eta"]
    cs = ax.contour(M, PR, ETA, levels=[0.6, 0.65, 0.7, 0.72, 0.74, 0.75], colors=SERIES_1, linewidths=0.7)
    ax.clabel(cs, fontsize=7, fmt="%.2f")
    ax.plot([r["mdot_c_corr"] for r in rows], [r["pr_c"] for r in rows], color=INK, marker="o", ms=4, lw=1.5,
            label="WOT operating line (cycle mean)")
    for r in rows:
        ax.annotate(f"{r['rpm']:.0f}", (r["mdot_c_corr"], r["pr_c"]), textcoords="offset points", xytext=(4, -10),
                    fontsize=7, color=INK)
    ax.set_xlabel("Corrected mass flow [kg/s]", color=INK_2)
    ax.set_ylabel("Total pressure ratio [-]", color=INK_2)
    ax.set_xlim(0.0, 0.3)
    ax.set_ylim(1.0, 2.6)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def plot_tipin(tp, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = tp["log"]
    t = log[:, 0] - tp["t_step"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), facecolor=SURFACE, sharex=True)
    _style(axes[0], "Intake plenum pressure [bar abs]")
    axes[0].plot(t, log[:, 1] / 1e5, color=GRID, lw=0.6, label="instantaneous")
    axes[0].plot(tp["t"], tp["map"] / 1e5, color=SERIES_1, lw=2, marker="o", ms=2.5, label="cycle mean")
    axes[0].axhline(tp["target"] / 1e5, color=INK_2, ls="--", lw=1, label="target")
    if np.isfinite(tp["t90_rise"]):
        axes[0].axvline(tp["t90_rise"], color=SERIES_2, lw=1, label=f"90 % of rise: {tp['t90_rise']:.2f} s")
    axes[0].legend(fontsize=8, frameon=False)
    _style(axes[1], "Turbo speed [krpm] and wastegate opening [-] x 100")
    axes[1].plot(t, log[:, 3] * 60.0 / (2.0 * math.pi) / 1e3, color=SERIES_1, lw=1.5, label="turbo speed")
    axes[1].plot(t, log[:, 4] * 100.0, color=SERIES_2, lw=1.2, label="wastegate opening x 100")
    axes[1].legend(fontsize=8, frameon=False)
    _style(axes[2], "Brake torque, cycle mean [Nm]")
    axes[2].plot(tp["t"], tp["torque"], color=SERIES_1, lw=2, marker="o", ms=2.5)
    axes[2].set_xlabel(f"Time from tip-in [s] (engine held at {tp['rpm']:.0f} rpm)", color=INK_2)
    fig.suptitle("Redline I4-2.0T tip-in, part load to WOT: turbo lag (MODELLED, synthetic maps)", color=INK,
                 fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    ts = turbo_spec()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sweep", default=None, help="start:stop:step in rpm (default spec/turbo.yaml sweep.default_rpm)")
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--no-tipin", action="store_true")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--out", default=str(OUT_DIR))
    a = ap.parse_args(argv)
    sw = [float(v) for v in a.sweep.split(":")] if a.sweep else ts["sweep"]["default_rpm"]
    rpms = sweep_speeds(*sw)
    rows = run_sweep(rpms, a.jobs)
    out = Path(a.out)
    write_csv(rows, out / "wot_turbo.csv")
    print("MODELLED full-load curve of the Redline I4-2.0T (no physical engine; SYNTHETIC compressor/turbine maps)")
    print(f"{'rpm':>5} {'T[Nm]':>6} {'P[kW]':>6} {'MAP':>5} {'tgt':>5} {'err%':>6} {'Nt[k]':>6} {'WG':>5} {'PRc':>5} "
          f"{'etaC':>5} {'SM':>5} {'p03':>5} {'T03':>5} {'etaT':>5} {'spark':>6} {'KLSA':>5} {'MBT':>5} {'BSFC':>5}")
    for r in rows:
        print(f"{r['rpm']:5.0f} {r['torque_nm']:6.1f} {r['power_kw']:6.1f} {r['map_bar']:5.2f} {r['target_bar']:5.2f} "
              f"{100 * r['boost_error']:6.1f} {r['turbo_krpm']:6.1f} {r['wg_opening']:5.2f} {r['pr_c']:5.2f} "
              f"{r['eta_c']:5.2f} {r['surge_margin']:5.2f} {r['p03_bar']:5.2f} {r['T03_K']:5.0f} {r['eta_t']:5.2f} "
              f"{r['spark_deg']:6.1f} {r['klsa_deg']:5.1f} {r['mbt_deg']:5.1f} {r['bsfc_g_per_kwh']:5.0f}"
              f"{'' if r['converged'] else '  NOT CONVERGED'}")
    it = max(range(len(rows)), key=lambda i: rows[i]["torque_nm"])
    ip = max(range(len(rows)), key=lambda i: rows[i]["power_kw"])
    tb = ts["variant"]["targets"]["torque_plateau_nm"]
    pl = plateau(rows, tb)
    print(f"peak torque {rows[it]['torque_nm']:.1f} Nm at {rows[it]['rpm']:.0f} rpm; peak power "
          f"{rows[ip]['power_kw']:.1f} kW at {rows[ip]['rpm']:.0f} rpm; min BSFC "
          f"{min(r['bsfc_g_per_kwh'] for r in rows):.0f} g/kWh")
    if pl:
        print(f"torque inside the {tb[0]:.0f}-{tb[1]:.0f} Nm band at {[int(r['rpm']) for r in pl]} rpm "
              f"({min(r['torque_nm'] for r in pl):.0f}-{max(r['torque_nm'] for r in pl):.0f} Nm)")
    else:
        print(f"no speed inside the {tb[0]:.0f}-{tb[1]:.0f} Nm design-intent band")
    held = [r for r in rows if abs(r["boost_error"]) <= float(ts["boost_control"]["hold_tol"])]
    print(f"boost target held within +/-{100 * float(ts['boost_control']['hold_tol']):.0f} % at "
          f"{[int(r['rpm']) for r in held]} rpm; wastegate closed below target at "
          f"{[int(r['rpm']) for r in rows if r['boost_error'] < -float(ts['boost_control']['hold_tol'])]} rpm")
    print(f"min compressor surge margin (cycle mean, m/m_surge - 1) {min(r['surge_margin'] for r in rows):.2f}; "
          f"knock-limited at {[int(r['rpm']) for r in rows if r['knock_limited']]} rpm; max knock integral at "
          f"x_b limit {max(r['knock_margin'] for r in rows):.3f} (threshold {KNOCK_THRESHOLD})")
    if not a.no_plot:
        plot_wot(rows, out / "wot_turbo.png", ts)
        plot_compressor(rows, out / "compressor_map.png")
    if not a.no_tipin:
        rpm_t = float(ts["tipin"]["rpm"])
        ref = min(rows, key=lambda r: abs(r["rpm"] - rpm_t))
        spark = ref["spark_deg"] if abs(ref["rpm"] - rpm_t) < 1.0 else _point(rpm_t)["spark_deg"]
        tp = run_tipin(spark)
        with open(out / "tipin.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["t_s", "map_bar", "torque_nm"])
            for t, m, q in zip(tp["t"], tp["map"], tp["torque"]):
                w.writerow([fmt(t), fmt(m / 1e5), fmt(q)])
        print(f"tip-in at {tp['rpm']:.0f} rpm (spark {spark:.1f} deg BTDC): MAP {tp['map_pre'] / 1e5:.2f} -> "
              f"{tp['map_final'] / 1e5:.2f} bar (target {tp['target'] / 1e5:.2f}); time to 90 % of the rise "
              f"{tp['t90_rise']:.2f} s, to 90 % of target MAP {tp['t90_abs']:.2f} s; torque {tp['torque_pre']:.0f} -> "
              f"{tp['torque_final']:.0f} Nm; turbo {tp['speed_final_krpm']:.0f} krpm after {tp['window_s']:.2f} s; min cycle-mean "
              f"surge margin {tp['surge_margin_cycle_min']:.2f}, instantaneous surge excursions {tp['surge_time']:.3f} s")
        if not a.no_plot:
            plot_tipin(tp, out / "tipin.png")
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    main()
