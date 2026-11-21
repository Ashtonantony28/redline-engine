"""WLTP Class 3b fuel consumption and CO2 of the generic compact sedan with the Redline engine.

    python -m vehicle.wltp

prints L/100 km and g CO2/km (overall and per phase) and writes out/vehicle/wltp.json and
out/vehicle/wltp.png. Backward simulation on the 1 Hz GTR 15 trace (vehicle.model.simulate_cycle),
WLTP test mass from GTR 15, economy shift strategy (not the GTR 15 Annex 2 gear-shift algorithm),
fuel from the Willans fit to the M7 maps, overrun fuel cut and idle at standstill. Hot start: no cold-
start enrichment and no catalyst heating, which a real type-approval test includes.
ALL NUMBERS ARE MODELLED (no physical car or engine).
"""
from __future__ import annotations

import argparse
import json
from collections import Counter

import numpy as np

from redline.spec import ROOT
from vehicle import params as vp
from vehicle.model import Powertrain, simulate_cycle
from vehicle.wltc import load_trace, phase_distances

OUT_DIR = ROOT / "out" / "vehicle"
M_PER_KM = 1.0e3
L_PER_100KM_SCALE = 100.0
G_PER_KG = 1.0e3
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"


def co2_per_kg_fuel(vs: dict | None = None) -> float:
    """kg CO2 per kg fuel from complete combustion of C_x H_y (all carbon to CO2)."""
    f = (vs or vp.vehicle_spec())["fuel"]
    m_fuel = f["carbon_atoms"] * f["atomic_mass_C"] + f["hydrogen_atoms"] * f["atomic_mass_H"]
    m_co2 = f["carbon_atoms"] * (f["atomic_mass_C"] + 2.0 * f["atomic_mass_O"])
    return m_co2 / m_fuel


def consumption(fuel_kg: float, dist_m: float, vs: dict | None = None) -> dict:
    vs = vs or vp.vehicle_spec()
    f = vs["fuel"]
    km = dist_m / M_PER_KM
    lhv = float(vp.engine_spec()["combustion"]["lhv"])
    gas_l = fuel_kg * lhv / float(f["gasoline_lhv"]) / f["gasoline_density"]
    return {
        "fuel_kg": fuel_kg,
        "l_per_100km": fuel_kg / f["density"] / km * L_PER_100KM_SCALE,
        "g_co2_per_km": fuel_kg * co2_per_kg_fuel(vs) * G_PER_KG / km,
        "gasoline_equiv_l_per_100km": gas_l / km * L_PER_100KM_SCALE,
        "gasoline_equiv_g_co2_per_km": gas_l * f["gasoline_co2_per_litre"] * G_PER_KG / km,
    }


def run(vs: dict | None = None) -> dict:
    vs = vs or vp.vehicle_spec()
    t, v = load_trace(vs)
    pt = Powertrain(vs=vs)
    mass = vp.wltp_test_mass(vs)
    res = simulate_cycle(t, v, pt, mass)
    out = {"modelled": True, "test_mass_kg": mass, "distance_km": res.distance / M_PER_KM,
           "duration_s": float(t[-1] - t[0]), **consumption(res.fuel_mass, res.distance, vs)}
    phases = {}
    pd = phase_distances(t, v, vs)
    for name, ph in vs["simulation"]["wltc_phases"].items():
        a, b = ph["span_s"]
        m = (res.t > a) & (res.t < b)
        phases[name] = {"distance_km": pd[name] / M_PER_KM,
                        **consumption(float(np.sum(res.fuel_rate[m]) * res.dt), pd[name], vs)}
    out["phases"] = phases
    modes = Counter(res.mode)
    p_eng = res.torque * res.rpm * vp.RAD_S_PER_RPM
    pos = p_eng > 0.0
    e_eng = float(np.sum(p_eng[pos]) * res.dt)
    out.update({
        "mode_seconds": dict(modes),
        "gear_seconds": {str(g): int(c) for g, c in sorted(Counter(res.gear.tolist()).items())},
        "idle_fuel_share": float(np.sum(res.fuel_rate[[m == "idle" for m in res.mode]]) * res.dt / res.fuel_mass),
        "positive_wheel_energy_kwh": float(np.sum(np.clip(res.p_wheel, 0.0, None)) * res.dt / 3.6e6),
        "positive_engine_energy_kwh": e_eng / 3.6e6,
        "cycle_mean_bsfc_g_per_kwh": float(np.sum(res.fuel_rate[pos]) * res.dt / e_eng * 3.6e6 * G_PER_KG),
        "max_trace_shortfall_N": float(res.shortfall.max()),
        "seconds_with_shortfall": int(np.sum(res.shortfall > 0.0)),
        "max_engine_rpm": float(res.rpm.max()),
        "fuel_share_inside_m7_map_coverage": float(np.sum(res.fuel_rate[pt.emap.in_map_coverage(res.rpm, res.torque)])
                                                   / np.sum(res.fuel_rate)),
        "willans_fit": {k: v for k, v in pt.emap.fit_report().items() if "err" in k},
    })
    return out, res


def plot(res, summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from vehicle.engine_map import engine_map

    em = engine_map()
    fig = plt.figure(figsize=(12, 8), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.6, 1.0])
    axs = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    ax_map = fig.add_subplot(gs[:, 1])
    axs[0].plot(res.t, res.v * vp.KMH_PER_MS, color=C1, lw=1.2)
    axs[0].set_title("WLTC Class 3b speed [km/h]", loc="left", fontsize=10, color=INK)
    axs[1].step(res.t, res.gear, color=C2, lw=1.0, where="mid")
    axs[1].set_title("gear (0 = clutch open, idle)", loc="left", fontsize=10, color=INK)
    axs[2].plot(res.t, res.fuel_rate * 3.6e6, color=C3, lw=1.0)
    axs[2].set_title("fuel flow [g/h]", loc="left", fontsize=10, color=INK)
    axs[2].set_xlabel("time [s]", color=INK_2)
    for ph in vp.vehicle_spec()["simulation"]["wltc_phases"].values():
        for ax in axs:
            ax.axvline(ph["span_s"][1], color=GRID, lw=1.0)
    n = np.linspace(em.idle_rpm, em.redline_rpm, 120)
    tq = np.linspace(1.0, float(em.wot_torque.max()), 100)
    N, TQ = np.meshgrid(n, tq)
    B = em.bsfc(N, TQ)
    B = np.where(TQ <= em.max_torque(N), B, np.nan)
    cs = ax_map.contour(N, TQ, B, levels=[250, 260, 280, 300, 350, 400, 500, 700], colors=INK_2, linewidths=0.7)
    ax_map.clabel(cs, fontsize=7, fmt="%.0f")
    ax_map.plot(em.wot_rpm, em.wot_torque, color=INK, lw=1.6, label="WOT (M4 sim.wot fixture)")
    drive = res.torque > 0
    ax_map.scatter(res.rpm[drive], res.torque[drive], s=4, color=C1, alpha=0.35, label="WLTP operating points")
    ax_map.set_title("BSFC [g/kWh], Willans fit to M7 maps", loc="left", fontsize=10, color=INK)
    ax_map.set_xlabel("engine speed [rpm]", color=INK_2)
    ax_map.set_ylabel("brake torque [N m]", color=INK_2)
    ax_map.legend(frameon=False, fontsize=8, loc="upper left")
    for ax in axs + [ax_map]:
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.5)
        ax.tick_params(colors=INK_2, labelsize=8)
    fig.suptitle(f"Redline compact sedan, WLTP Class 3b (MODELLED): {summary['l_per_100km']:.2f} L/100 km, "
                 f"{summary['g_co2_per_km']:.0f} g CO2/km", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=120, facecolor=SURFACE)
    plt.close(fig)


def summary_text(s: dict) -> str:
    L = ["Redline M10 WLTP Class 3b. ALL NUMBERS MODELLED (generic compact sedan, Redline I4-2.0 NA, 6MT).",
         f"  test mass {s['test_mass_kg']:.0f} kg, distance {s['distance_km']:.3f} km, duration {s['duration_s']:.0f} s",
         f"  fuel {s['fuel_kg']:.3f} kg iso-octane -> {s['l_per_100km']:.2f} L/100 km, {s['g_co2_per_km']:.1f} g CO2/km",
         f"  energy-equivalent E5 gasoline: {s['gasoline_equiv_l_per_100km']:.2f} L/100 km, "
         f"{s['gasoline_equiv_g_co2_per_km']:.1f} g CO2/km", "  phases:"]
    for k, p in s["phases"].items():
        L.append(f"    {k:10s} {p['distance_km']:6.3f} km  {p['l_per_100km']:5.2f} L/100 km  {p['g_co2_per_km']:6.1f} g/km")
    L += [f"  positive wheel energy {s['positive_wheel_energy_kwh']:.3f} kWh, engine {s['positive_engine_energy_kwh']:.3f} kWh,"
          f" cycle-mean BSFC {s['cycle_mean_bsfc_g_per_kwh']:.0f} g/kWh, idle fuel share {s['idle_fuel_share'] * 100:.1f} %",
          f"  modes [s]: {s['mode_seconds']}", f"  gears [s]: {s['gear_seconds']}",
          f"  trace shortfall: {s['seconds_with_shortfall']} s (max {s['max_trace_shortfall_N']:.0f} N);"
          f" max engine speed {s['max_engine_rpm']:.0f} rpm",
          f"  Willans fit to M7 maps: throttled rms {s['willans_fit']['throttled_rms_rel_err'] * 100:.2f} %"
          f" (max {s['willans_fit']['throttled_max_abs_rel_err'] * 100:.2f} %), WOT max"
          f" {s['willans_fit']['wot_max_abs_rel_err'] * 100:.2f} %",
          f"  fuel burned inside the speed/BMEP range sampled by the M7 maps: "
          f"{s['fuel_share_inside_m7_map_coverage'] * 100:.0f} % (rest: Willans extrapolation to light load and idle)",
          "  hot start, no cold-start enrichment; economy shift strategy, not GTR 15 Annex 2"]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args(argv)
    s, res = run()
    print(summary_text(s))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "wltp.json").write_text(json.dumps(s, indent=2))
    if not a.no_plot:
        plot(res, s, OUT_DIR / "wltp.png")
        print(f"  -> {OUT_DIR / 'wltp.json'}, {OUT_DIR / 'wltp.png'}")


if __name__ == "__main__":
    main()
