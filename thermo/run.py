"""M2 runner: `python -m thermo.run --rpm 3000 --load full`.

Runs the converged single-cylinder cycle, prints a summary and writes out/thermo/:
  cycle_<rpm>_<load>.csv     crank-angle traces (0..720 deg, ADR 0001)
  summary_<rpm>_<load>.csv   headline metrics
  energy_<rpm>_<load>.csv    closed-cycle energy balance
  pv_<rpm>_<load>.png        P-V diagram (linear and log-log)
All numbers are MODELLED: there is no physical engine.
"""
from __future__ import annotations

import argparse
import csv

import numpy as np

from redline.spec import ROOT
from thermo.cycle import CYCLE_DEG, energy_balance, run_cycle
from thermo.props import build_tables, thermo_spec

OUT_DIR = ROOT / "out" / "thermo"
# chart tokens (dataviz reference palette: surface, ink, categorical slot 1)
SURFACE, INK, INK_2, SERIES_1, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#2a78d6", "#e4e3df"
PA_PER_BAR, M3_PER_CC = 1.0e5, 1.0e-6


def write_outputs(res, eb, tag):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    th = np.mod(res.theta, CYCLE_DEG)
    order = np.argsort(th, kind="stable")
    cols = {
        "theta_deg": th, "V_m3": res.V, "p_Pa": res.p, "T_u_K": res.T_u, "T_b_K": res.T_b,
        "T_bulk_K": res.T_bulk, "x_b": res.x_b, "hrr_chem_J_per_deg": res.hrr_chem,
        "hrr_apparent_J_per_deg": res.hrr_apparent,
    }
    with open(OUT_DIR / f"cycle_{tag}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(cols))
        for i in order:
            w.writerow([f"{cols[k][i]:.10g}" for k in cols])
    with open(OUT_DIR / f"summary_{tag}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in res.metrics.items():
            w.writerow([k, v])
    with open(OUT_DIR / f"energy_{tag}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["term", "J"])
        for k, v in eb.items():
            w.writerow([k, v])


def plot_pv(res, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    V = res.V / M3_PER_CC
    p = res.p / PA_PER_BAR
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), facecolor=SURFACE)
    for ax, log in zip(axes, (False, True)):
        ax.set_facecolor(SURFACE)
        ax.plot(V, p, color=SERIES_1, lw=2)
        if log:
            ax.set_xscale("log")
            ax.set_yscale("log")
        ax.set_xlabel("Cylinder volume [cm³]", color=INK_2)
        ax.set_ylabel("Pressure [bar]", color=INK_2)
        ax.set_title("log-log" if log else "linear", color=INK, loc="left", fontsize=10)
        ax.grid(True, color=GRID, lw=0.6, which="both")
        for s in ax.spines.values():
            s.set_color(GRID)
        ax.tick_params(colors=INK_2)
    m = res.metrics
    fig.suptitle(f"Redline I4 single cylinder, {res.rpm:.0f} rpm, {res.load} load (MODELLED): "
                 f"IMEPn {m['imep_net_bar']:.2f} bar, pmax {m['p_max_bar']:.1f} bar",
                 color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"pv_{tag}.png", dpi=130, facecolor=SURFACE)
    plt.close(fig)


def summary_text(res, eb):
    m = res.metrics
    lines = [
        f"Redline M2 two-zone cycle (MODELLED, no physical engine): {res.rpm:.0f} rpm, load={res.load}, "
        f"step {res.step_deg} deg, {m['cycles']} cycles, converged={m['converged']}",
        f"  IMEP gross / net        {m['imep_gross_bar']:.3f} / {m['imep_net_bar']:.3f} bar  (PMEP {m['pmep_bar']:.3f})",
        f"  indicated eff. gross/net {m['eta_ind_gross']:.4f} / {m['eta_ind_net']:.4f}  (LHV {m['lhv_spec']:.4g} J/kg)",
        f"  peak pressure           {m['p_max_bar']:.2f} bar at {m['theta_p_max_deg']:.1f} deg ATDC",
        f"  CA10 / CA50 / CA90      {m['ca10_deg']:.2f} / {m['ca50_deg']:.2f} / {m['ca90_deg']:.2f} deg",
        f"  T_b max / T_u max       {m['T_b_max_K']:.0f} / {m['T_u_max_K']:.0f} K",
        f"  residual gas fraction   {m['residual_fraction']:.4f}",
        f"  volumetric efficiency   {m['volumetric_efficiency']:.4f}",
        f"  trapped mass / fuel     {m['m_trapped_kg']*1e6:.2f} / {m['m_fuel_kg']*1e6:.3f} mg",
        f"  IVC state               {m['p_ivc_bar']:.3f} bar, {m['T_ivc_K']:.1f} K",
        f"  wall heat closed / open {m['Q_wall_closed_J']:.1f} / {m['Q_wall_open_J']:.1f} J",
        f"  energy balance IVC-EVO  fuel {eb['E_fuel']:.1f} J = W {eb['W']:.1f} + Q {eb['Q_wall']:.1f} "
        f"+ dU_s {eb['dU_sensible']:.1f} + incompl. {eb['E_incomplete']:.2f}; residual {eb['rel_residual']:.2e}",
        f"  mass through cycle      fresh in {m['fresh_in_kg']*1e6:.4f} mg, exhaust out {m['exhaust_out_net_kg']*1e6:.4f} mg",
    ]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rpm", type=float, default=3000.0)
    ap.add_argument("--load", default="full", choices=sorted(thermo_spec()["loads"]))
    ap.add_argument("--step", type=float, default=None, help="crank-angle step [deg]")
    ap.add_argument("--no-heat-transfer", action="store_true")
    a = ap.parse_args(argv)
    tables = build_tables()
    res = run_cycle(a.rpm, a.load, a.step, heat_transfer=not a.no_heat_transfer, tables=tables)
    eb = energy_balance(res.closed, tables)
    tag = f"{a.rpm:.0f}_{a.load}"
    write_outputs(res, eb, tag)
    plot_pv(res, tag)
    print(summary_text(res, eb))
    print(f"  wrote {OUT_DIR}/{{cycle,summary,energy}}_{tag}.csv and pv_{tag}.png")
    return res


if __name__ == "__main__":
    main()
