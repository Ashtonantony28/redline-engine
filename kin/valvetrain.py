"""M3 valvetrain analysis: cams, VVT, piston-to-valve clearance, dynamics, Hertz stress.

``python -m kin.valvetrain`` prints the report and writes ``out/valvetrain/*.csv|png``.

Every number is MODELLED: there is no physical engine (CLAUDE.md rule 5). SI units, crank
angle 0 = firing TDC of cylinder 1 (ADR 0001). Engine-level parameters come from
``spec/engine.yaml``; valvetrain-only parameters from ``spec/valvetrain.yaml``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
import math
from pathlib import Path

import numpy as np

from kin import slider_crank as sc
from kin import valve_dynamics as vd
from kin.cam import CamProfile, TWO_PI, design_cams
from redline.spec import ROOT, geometry, load

VT_SPEC_PATH = ROOT / "spec" / "valvetrain.yaml"
OUT_DIR = ROOT / "out" / "valvetrain"
CYCLE = 2.0 * TWO_PI                     # one four-stroke cycle in crank radians
VALVES = ("intake", "exhaust")


def load_specs():
    return load(), load(str(VT_SPEC_PATH))


# ----------------------------------------------------------------------------- geometry checks
def required_tappet_diameter(profile: CamProfile, edge_margin: float, samples: int) -> float:
    """Flat follower: contact sits |dy/dphi| off the tappet axis; keep it inside the face."""
    phi = np.linspace(-profile.half_event, profile.half_event, samples)
    return 2.0 * (float(np.abs(profile.contact_offset(phi)).max()) + edge_margin)


def coil_bind_margin(vt: dict, max_valve_lift: float) -> float:
    sp = vt["spring"]
    return sp["installed_height"] - max_valve_lift - sp["solid_height"]


# ----------------------------------------------------------------------------- VVT + clearance
def vvt_advances(engine: dict, vt: dict) -> np.ndarray:
    """Intake phaser positions [crank rad], 0 (parked) .. full advance."""
    span = math.radians(engine["valves"]["vvt_intake_range_deg"])
    return np.linspace(0.0, span, int(vt["simulation"]["vvt_steps"]))


def valve_pocket_gap(theta, lift, piston_s, valve_angle, closed_gap, pocket):
    """Valve face to pocket floor, along the valve axis.

    The pocket floor is a plane normal to the valve axis; the piston moving down by ``s``
    along the bore axis opens the gap by ``s cos(alpha)``; the valve opening closes it by
    its lift. ``closed_gap + pocket`` is the gap at TDC with the valve shut."""
    return closed_gap + pocket + piston_s * math.cos(valve_angle) - lift


def clearance_map(profile: CamProfile, which: str, engine: dict, vt: dict, advances):
    """Gap [m] on a (len(advances), crank_samples) grid over one 720 deg cycle."""
    g = geometry(engine)
    cc = vt["combustion_chamber"]
    theta = np.linspace(0.0, CYCLE, int(vt["simulation"]["crank_samples"]))
    s = sc.displacement_from_tdc(theta, g.crank_radius, g.conrod_length)
    alpha = math.radians(cc[f"{which}_valve_angle_deg"])
    gaps = np.array([
        valve_pocket_gap(theta, profile.valve_lift_crank(theta, adv), s, alpha,
                         cc["valve_to_crown_at_tdc"], cc[f"{which}_pocket_depth"])
        for adv in np.atleast_1d(advances)
    ])
    return theta, gaps


# ----------------------------------------------------------------------------- Hertz
def reduced_modulus(e1, nu1, e2, nu2):
    """E* with 1/E* = (1 - nu1^2)/E1 + (1 - nu2^2)/E2."""
    return 1.0 / ((1.0 - nu1**2) / e1 + (1.0 - nu2**2) / e2)


def hertz_line_pressure(force, radius, width, e_star):
    """Peak Hertz pressure of a cylinder (radius R) on a flat, line contact of length w:
    p0 = sqrt(P' E* / (pi R)), P' = F / w (Johnson, Contact Mechanics, eq. 4.45-4.46)."""
    force = np.maximum(np.asarray(force, dtype=float), 0.0)
    return np.sqrt(force * e_star / (math.pi * np.asarray(radius, dtype=float) * width))


def hertz_contact_halfwidth(force, radius, width, e_star):
    """b = sqrt(4 P' R / (pi E*))."""
    return np.sqrt(4.0 * np.asarray(force) / width * radius / (math.pi * e_star))


# ----------------------------------------------------------------------------- analysis
@dataclass
class ValveResult:
    name: str
    profile: CamProfile
    model: vd.ValvetrainModel
    seat_open_deg: float
    seat_close_deg: float
    min_rho: float
    max_offset: float
    max_accel: float
    min_accel: float
    max_jerk: float
    required_tappet_diameter: float
    coil_bind_margin: float
    min_clearance: float
    min_clearance_theta_deg: float
    min_clearance_advance_deg: float
    qs_float_rpm: float
    float_rpm: float
    surge_hz: np.ndarray
    surge_continuous_hz: np.ndarray
    system_hz: np.ndarray
    peak_hertz: float
    peak_hertz_source: str
    seat_speed_redline: float
    sweep: vd.SimResult = field(repr=False)
    clearance_theta: np.ndarray = field(repr=False)
    clearance_min_per_advance: np.ndarray = field(repr=False)
    advances: np.ndarray = field(repr=False)


def analyse_valve(name, profile, engine, vt) -> ValveResult:
    sim = vt["simulation"]
    samples = int(sim["crank_samples"])
    phi = np.linspace(-profile.half_event, profile.half_event, samples)
    d = profile.derivatives(phi)
    rho = profile.radius_of_curvature(phi)
    o, c = profile.valve_open_crank()

    advances = vvt_advances(engine, vt) if name == "intake" else np.zeros(1)
    theta, gaps = clearance_map(profile, name, engine, vt, advances)
    i_adv, i_th = np.unravel_index(np.argmin(gaps), gaps.shape)

    model = vd.model_from_spec(name, vt)
    revs, safety = int(sim["revolutions"]), sim["step_safety"]
    float_rpm, sweep = vd.float_speed(model, profile, sim["sweep_rpm"], sim["refine_rpm_step"], revs, safety)

    # Hertz: quasi-static single-mass over the idle..redline band (dominant at low speed,
    # where inertia does not unload the nose) and the multi-mass sweep within the band.
    m = vt["materials"]
    e_star = reduced_modulus(m["cam"]["E"], m["cam"]["nu"], m["tappet"]["E"], m["tappet"]["nu"])
    width = vt["cam"]["lobe_width"]
    lo, hi = vt["simulation"]["hertz_rpm"]
    qs_rpm = np.linspace(lo, hi, int(sim["vvt_steps"]))
    qs_p = max(float(hertz_line_pressure(vd.quasi_static_contact_force(model, profile, phi, r), rho, width,
                                         e_star).max()) for r in qs_rpm)
    band = (sweep.rpm >= lo) & (sweep.rpm <= hi)
    dyn_p = float(np.sqrt(sweep.max_stress_ratio[band] * e_star / (math.pi * width)).max()) if band.any() else 0.0
    peak, src = (qs_p, "quasi-static") if qs_p >= dyn_p else (dyn_p, "multi-mass")
    redline = engine["engine"]["redline_rpm"]
    seat = float(np.interp(redline, sweep.rpm, sweep.max_seat_speed))

    return ValveResult(
        name=name, profile=profile, model=model,
        seat_open_deg=math.degrees(o), seat_close_deg=math.degrees(c),
        min_rho=float(rho.min()), max_offset=float(np.abs(d[1]).max()),
        max_accel=float(d[2].max()), min_accel=float(d[2].min()), max_jerk=float(np.abs(d[3]).max()),
        required_tappet_diameter=required_tappet_diameter(profile, vt["tappet"]["edge_margin"], samples),
        coil_bind_margin=coil_bind_margin(vt, profile.valve_max_lift),
        min_clearance=float(gaps[i_adv, i_th]),
        min_clearance_theta_deg=math.degrees(theta[i_th]),
        min_clearance_advance_deg=math.degrees(advances[i_adv]),
        qs_float_rpm=vd.quasi_static_float_rpm(model, profile, samples),
        float_rpm=float_rpm,
        surge_hz=model.spring_chain_frequencies(),
        surge_continuous_hz=vd.surge_frequencies_continuous(model.spring_rate, model.spring_mass),
        system_hz=model.natural_frequencies(contact=True, seat=False),
        peak_hertz=peak, peak_hertz_source=src, seat_speed_redline=seat,
        sweep=sweep, clearance_theta=theta, clearance_min_per_advance=gaps.min(axis=1),
        advances=advances,
    )


def analyse(engine=None, vt=None):
    if engine is None or vt is None:
        engine, vt = load_specs()
    cams = design_cams(engine, vt)
    return {n: analyse_valve(n, cams[n], engine, vt) for n in VALVES}


# ----------------------------------------------------------------------------- report
def _deg(x):
    return math.degrees(x)


def report(results, engine, vt) -> str:
    lim = vt["limits"]
    lines = ["Redline M3 valvetrain report (ALL VALUES MODELLED - no physical engine)", ""]
    for r in results.values():
        p = r.profile
        lines += [
            f"[{r.name}] valve lift {p.valve_max_lift * 1e3:.2f} mm, lash {p.lash * 1e3:.2f} mm, "
            f"base circle r {p.base_radius * 1e3:.1f} mm",
            f"  main half-angle {_deg(p.main_half_angle):.2f} cam deg, ramp {_deg(p.ramp_angle):.2f} cam deg "
            f"(lift at spec timing = {vt['timing']['check_lift'] * 1e3:.1f} mm)",
            f"  seat-to-seat {r.seat_open_deg:.1f} .. {r.seat_close_deg:.1f} crank deg "
            f"({r.seat_close_deg - r.seat_open_deg:.1f} deg)",
            f"  accel +{r.max_accel * 1e3:.2f} / {r.min_accel * 1e3:.2f} mm/rad^2, |jerk| max "
            f"{r.max_jerk * 1e3:.1f} mm/rad^3",
            f"  min radius of curvature {r.min_rho * 1e3:.2f} mm (limit {lim['min_radius_of_curvature'] * 1e3:.1f})",
            f"  max contact offset {r.max_offset * 1e3:.2f} mm -> tappet diameter needed "
            f"{r.required_tappet_diameter * 1e3:.2f} mm (fitted {vt['tappet']['diameter'] * 1e3:.1f})",
            f"  coil-bind margin {r.coil_bind_margin * 1e3:.2f} mm (min {vt['spring']['min_coil_bind_margin'] * 1e3:.1f})",
            f"  min piston-valve clearance {r.min_clearance * 1e3:.3f} mm at {r.min_clearance_theta_deg:.1f} crank deg, "
            f"VVT advance {r.min_clearance_advance_deg:.1f} deg (limit {lim['min_piston_valve_clearance'] * 1e3:.2f})",
            f"  spring surge (fixed-fixed, {r.model.n_spring}-mass): "
            + ", ".join(f"{f:.0f}" for f in r.surge_hz[:3])
            + " Hz; continuous: " + ", ".join(f"{f:.0f}" for f in r.surge_continuous_hz) + " Hz",
            "  system modes (contact closed, valve off seat): " + ", ".join(f"{f:.0f}" for f in r.system_hz) + " Hz",
            f"  float speed: single-mass quasi-static {r.qs_float_rpm:.0f} rpm, multi-mass dynamic "
            f"{r.float_rpm:.0f} rpm (required > {lim['float_free_rpm']})",
            f"  seating speed at redline {r.seat_speed_redline:.3f} m/s",
            f"  peak Hertz pressure {r.peak_hertz / 1e6:.0f} MPa ({r.peak_hertz_source}; limit "
            f"{lim['hertz_limit'] / 1e6:.0f} MPa)",
            "",
        ]
    return "\n".join(lines)


def write_outputs(results, engine, vt, out_dir: Path = OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = int(vt["simulation"]["crank_samples"])
    for r in results.values():
        p = r.profile
        phi = np.linspace(-p.half_event, p.half_event, samples)
        d = p.derivatives(phi)
        with open(out_dir / f"cam_{r.name}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cam_deg", "cam_lift_m", "velocity_m_per_rad", "accel_m_per_rad2", "jerk_m_per_rad3",
                        "valve_lift_m", "rho_m"])
            for row in zip(np.degrees(phi), d[0], d[1], d[2], d[3], p.valve_lift(phi), p.radius_of_curvature(phi)):
                w.writerow([f"{v:.8g}" for v in row])
        s = r.sweep
        with open(out_dir / f"sweep_{r.name}.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["engine_rpm", "min_contact_main_N", "max_contact_N", "max_seat_speed_m_s", "float"])
            for row in zip(s.rpm, s.min_contact_main, s.max_contact, s.max_seat_speed, s.floats):
                w.writerow([f"{row[0]:.0f}", f"{row[1]:.3f}", f"{row[2]:.3f}", f"{row[3]:.4f}", int(row[4])])
    ri = results["intake"]
    with open(out_dir / "clearance_intake_vvt.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["vvt_advance_crank_deg", "min_clearance_m"])
        for a, g in zip(np.degrees(ri.advances), ri.clearance_min_per_advance):
            w.writerow([f"{a:.2f}", f"{g:.6g}"])
    _plots(results, engine, vt, out_dir)


SERIES = ("#2a78d6", "#eb6834", "#1baf7a")   # dataviz reference palette, slots 1-3 (light)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=INK, fontsize=10, loc="left")
    ax.set_xlabel(xlabel, color=MUTED, fontsize=9)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    ax.grid(color=GRID, linewidth=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def _plots(results, engine, vt, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    samples = int(vt["simulation"]["crank_samples"])
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    labels = ("lift [mm]", "velocity [mm/rad]", "acceleration [mm/rad^2]", "jerk [mm/rad^3]")
    for k, ax in enumerate(axes.flat):
        for colour, r in zip(SERIES, results.values()):
            p = r.profile
            phi = np.linspace(-p.half_event, p.half_event, samples)
            ax.plot(np.degrees(phi), p.derivatives(phi)[k] * 1e3, color=colour, linewidth=1.6, label=r.name)
        _style(ax, labels[k].split(" [")[0].capitalize() + " vs cam angle (modelled)", "cam angle from nose [deg]",
               labels[k])
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(out_dir / "cam_profiles.png", dpi=130)
    plt.close(fig)

    # VVT lift curves and piston-valve gap
    ri, re = results["intake"], results["exhaust"]
    theta = ri.clearance_theta
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    ax = axes[0]
    for colour, adv in zip(SERIES, (ri.advances[0], ri.advances[len(ri.advances) // 2], ri.advances[-1])):
        ax.plot(np.degrees(theta), ri.profile.valve_lift_crank(theta, adv) * 1e3, color=colour, linewidth=1.6,
                label=f"intake, advance {math.degrees(adv):.0f} deg")
    ax.plot(np.degrees(theta), re.profile.valve_lift_crank(theta) * 1e3, color=MUTED, linewidth=1.2,
            linestyle="--", label="exhaust")
    _style(ax, "Valve lift over the VVT range (modelled)", "crank angle [deg]", "valve lift [mm]")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    ax.plot(np.degrees(ri.advances), ri.clearance_min_per_advance * 1e3, color=SERIES[0], linewidth=1.6,
            marker="o", markersize=3, label="intake")
    ax.axhline(re.min_clearance * 1e3, color=SERIES[1], linewidth=1.6, label="exhaust (not phased)")
    ax.axhline(vt["limits"]["min_piston_valve_clearance"] * 1e3, color=MUTED, linestyle=":", linewidth=1.2,
               label="limit")
    _style(ax, "Minimum piston-valve clearance (modelled)", "intake VVT advance [crank deg]", "clearance [mm]")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out_dir / "vvt_clearance.png", dpi=130)
    plt.close(fig)

    # dynamics: min main-event contact force vs speed
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    for colour, r in zip(SERIES, results.values()):
        ax.plot(r.sweep.rpm, r.sweep.min_contact_main, color=colour, linewidth=1.6, label=r.name)
    ax.axvline(vt["limits"]["float_free_rpm"], color=MUTED, linestyle=":", linewidth=1.2, label="float-free limit")
    _style(ax, "Min cam/tappet force in the main event, multi-mass (modelled)", "engine speed [rpm]",
           "contact force [N]")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out_dir / "float_sweep.png", dpi=130)
    plt.close(fig)


def main():
    engine, vt = load_specs()
    results = analyse(engine, vt)
    print(report(results, engine, vt))
    write_outputs(results, engine, vt)
    print(f"wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
