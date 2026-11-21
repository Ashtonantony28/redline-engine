"""Full-throttle acceleration: 0-100 km/h and 80-120 km/h in gear.

    python -m vehicle.accel

Forward simulation (vehicle.model.simulate_wot) at mass in running order (driver included), WOT
torque from the M4 sim.wot fixture, clutch-slip launch, upshifts at the accel shift speed with a
torque interruption, front-axle traction limit with load transfer. Writes out/vehicle/accel.json and
out/vehicle/accel.png. ALL NUMBERS ARE MODELLED (no physical car or engine).
"""
from __future__ import annotations

import argparse
import json

from redline.spec import ROOT
from vehicle import params as vp
from vehicle.model import Powertrain, simulate_wot

OUT_DIR = ROOT / "out" / "vehicle"
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2 = "#2a78d6", "#eb6834"


def run(vs: dict | None = None):
    vs = vs or vp.vehicle_spec()
    tg = vs["simulation"]["accel_targets_kmh"]
    pt = Powertrain(vs=vs)
    mass = vs["vehicle"]["mass_running_order"]
    v100 = tg["zero_to"] / vp.KMH_PER_MS
    v80, v120 = tg["in_gear_from"] / vp.KMH_PER_MS, tg["in_gear_to"] / vp.KMH_PER_MS
    standing = simulate_wot(0.0, v100, 1, pt, mass, allow_shifts=True, launch=True)
    g = vs["simulation"]["in_gear_gear"]
    elastic = simulate_wot(v80, v120, g, pt, mass, allow_shifts=False, launch=False)
    peak_kw = float(pt.emap.max_power(pt.emap.wot_rpm).max()) / 1.0e3
    s = {"modelled": True, "mass_kg": mass, "power_to_mass_w_per_kg": peak_kw * 1.0e3 / mass,
         "peak_power_kw": peak_kw,
         "zero_to_100_s": standing.time_between(0.0, v100),
         "shifts_0_100": [{"t_s": t, "to_gear": gg} for t, gg in standing.shifts],
         "in_gear": g, "in_gear_80_120_s": elastic.time_between(v80, v120),
         "rpm_at_80_in_gear": pt.rpm(v80, g), "rpm_at_120_in_gear": pt.rpm(v120, g),
         "traction_limit_N": pt.traction_limit(mass)}
    return s, standing, elastic


def plot(standing, elastic, s, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), facecolor=SURFACE)
    axs[0].plot(standing.t, standing.v * vp.KMH_PER_MS, color=C1, lw=1.8)
    for t, g in standing.shifts:
        axs[0].axvline(t, color=GRID, lw=1.0)
        axs[0].annotate(f"{g}", (t, 5), color=INK_2, fontsize=8)
    axs[0].set_title(f"0-100 km/h: {s['zero_to_100_s']:.2f} s (MODELLED)", loc="left", fontsize=10, color=INK)
    axs[1].plot(elastic.t, elastic.v * vp.KMH_PER_MS, color=C2, lw=1.8)
    axs[1].set_title(f"80-120 km/h in gear {s['in_gear']}: {s['in_gear_80_120_s']:.2f} s (MODELLED)", loc="left",
                     fontsize=10, color=INK)
    for ax in axs:
        ax.set_facecolor(SURFACE)
        ax.set_xlabel("time [s]", color=INK_2)
        ax.set_ylabel("speed [km/h]", color=INK_2)
        ax.grid(True, color=GRID, lw=0.5)
        ax.tick_params(colors=INK_2, labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120, facecolor=SURFACE)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args(argv)
    s, standing, elastic = run()
    print("Redline M10 acceleration. ALL NUMBERS MODELLED (generic compact sedan, Redline I4-2.0 NA, 6MT).")
    print(f"  mass {s['mass_kg']:.0f} kg, peak {s['peak_power_kw']:.1f} kW ({s['power_to_mass_w_per_kg']:.1f} W/kg),"
          f" traction limit {s['traction_limit_N']:.0f} N")
    print(f"  0-100 km/h: {s['zero_to_100_s']:.2f} s  (shifts at "
          + ", ".join(f"{x['t_s']:.2f} s -> {x['to_gear']}" for x in s["shifts_0_100"]) + ")")
    print(f"  80-120 km/h in gear {s['in_gear']}: {s['in_gear_80_120_s']:.2f} s"
          f" ({s['rpm_at_80_in_gear']:.0f} -> {s['rpm_at_120_in_gear']:.0f} rpm)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "accel.json").write_text(json.dumps(s, indent=2))
    if not a.no_plot:
        plot(standing, elastic, s, OUT_DIR / "accel.png")
        print(f"  -> {OUT_DIR / 'accel.json'}, {OUT_DIR / 'accel.png'}")


if __name__ == "__main__":
    main()
