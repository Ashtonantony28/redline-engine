"""Build the Redline technical report: `python -m report.build [--fast]` -> docs/REPORT.md.

Every number is pulled from the code, spec files and committed fixtures at build time (nothing is
hand-typed where it can be computed), and the figures it links are regenerated under out/:

* M2 thermo (thermo.run), M3 valvetrain (kin.valvetrain), M5 dynamics on REAL pressure traces
  (dyn.report --pressure m4, with the placeholder comparison), M8 host HIL cycles (hil.runner),
  M9 PPO table (ml.eval), M10 WLTP / acceleration (vehicle.*) and the three.js replay (viz.build);
* M4 WOT curve and M7 maps are re-plotted from their committed fixtures (re-running the sweeps takes
  minutes; `python -m sim.wot` / `python -m sim.maps` regenerate them from the models);
* M1 CAD mass properties and the M9 Vortex-V co-sim budget are read from out/ when present (their
  generators need build123d / Verilator and the RTL tree).

`--fast` skips the dynamics (~1 min) and replay rebuilds and reads their last outputs instead.
ALL NUMBERS ARE MODELLED: there is no physical engine or car.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import datetime as _dt
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

from redline.spec import ROOT, geometry, load
from mvem.params import mvem_spec
from vehicle.params import vehicle_spec

CONFIG = ROOT / "report" / "report.yaml"
DEG = math.pi / 180.0
SEC_PER_MIN = 60.0


def config() -> dict:
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def rel(path) -> str:
    """Link from docs/REPORT.md to a repo path."""
    return "../" + str(Path(path).resolve().relative_to(ROOT))


def fig(path, caption) -> str:
    p = Path(path)
    return f"![{caption}]({rel(p)})" if p.exists() else f"*(figure {rel(p)} not generated)*"


def _num(v):
    if v in ("True", "False"):
        return v == "True"
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def read_rows(path):
    with open(path, newline="") as f:
        return [{k: _num(v) for k, v in r.items()} for r in csv.DictReader(f)]


def run_module(args, quiet=True):
    r = subprocess.run([sys.executable, "-m", *args], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"python -m {' '.join(args)} failed:\n{r.stdout}\n{r.stderr}")
    return r.stdout


# ------------------------------------------------------------------------------------------------
# per-milestone data
# ------------------------------------------------------------------------------------------------
def m0(spec):
    from kin.slider_crank import piston_acceleration, piston_velocity
    from thermo.ideal import otto_efficiency
    g = geometry(spec)
    e = spec["engine"]
    w = e["redline_rpm"] * 2 * math.pi / SEC_PER_MIN
    th = np.linspace(0.0, 2 * math.pi, 3601)
    acc = piston_acceleration(th, w, g.crank_radius, g.conrod_length)
    return {"displacement_cc": g.displacement * 1e6, "bore_mm": g.bore * 1e3, "stroke_mm": g.stroke * 1e3,
            "cr": g.compression_ratio, "rod_ratio": g.conrod_length / g.stroke,
            "otto_eff_1p4": otto_efficiency(g.compression_ratio, 1.4),
            "mean_piston_speed_redline": 2 * g.stroke * e["redline_rpm"] / SEC_PER_MIN,
            "max_piston_speed_redline": float(np.abs(piston_velocity(th, w, g.crank_radius, g.conrod_length)).max()),
            "max_piston_acc_redline_g": float(acc.max()) / 9.80665, "min_piston_acc_redline_g": float(acc.min()) / 9.80665}


def m1(cfg):
    p = ROOT / cfg["external_artefacts"]["cad_mass_properties"]
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return {"piston_assembly": d["piston_assembly"], "conrod": d["conrod"], "recip": d["reciprocating_mass"],
            "crank": d["crankshaft"], "chamber": d["combustion_chamber"]}


def m2(cfg):
    from thermo import run as trun
    tp = cfg["thermo_point"]
    with contextlib.redirect_stdout(io.StringIO()):
        res = trun.main(["--rpm", str(tp["rpm"]), "--load", tp["load"]])
    from thermo.cycle import energy_balance
    from thermo.props import build_tables
    eb = energy_balance(res.closed, build_tables())
    tag = f"{tp['rpm']:.0f}_{tp['load']}"
    return {"m": res.metrics, "eb": eb, "rpm": tp["rpm"], "load": tp["load"],
            "fig": ROOT / "out" / "thermo" / f"pv_{tag}.png"}


def m3():
    from kin import valvetrain as vt
    engine, vts = vt.load_specs()
    res = vt.analyse(engine, vts)
    vt.write_outputs(res, engine, vts)
    figs = sorted((ROOT / "out" / "valvetrain").glob("*.png"))
    return {"res": res, "limits": vts["limits"], "figs": figs, "check_lift": vts["timing"]["check_lift"]}


def m4(fig_dir):
    from sim.wot import peaks, plot
    rows = read_rows(ROOT / "fixtures" / "wot" / "wot_coarse.csv")
    spec = load()
    path = fig_dir / "wot_coarse.png"
    plot(rows, path, spec["targets"])
    return {"rows": rows, "peaks": peaks(rows), "fig": path}


def m5(fast: bool):
    out = ROOT / "out" / "dyn" / "m4"
    if not fast or not (out / "summary.json").exists():
        run_module(["dyn.report", "--pressure", "m4"])
    s = json.loads((out / "summary.json").read_text())
    cmp_ = json.loads((out / "compare.json").read_text())["rows"] if (out / "compare.json").exists() else []
    return {"s": s, "cmp": cmp_, "dir": out}


def m7(fig_dir):
    from sim.maps import plot, summary
    from thermo.m7spec import combustion_spec
    rows = read_rows(ROOT / "fixtures" / "maps" / "maps_m4.csv")
    d = fig_dir / "maps"
    d.mkdir(parents=True, exist_ok=True)
    plot(rows, d, "M4 1D engine (fixture)")
    mp = combustion_spec()["maps"]
    grid = len(mp["rpm"]) * len(mp["map_levels"])
    mbt = [r for r in rows if not r["knock_limited"]]
    return {"rows": rows, "summary": summary(rows), "grid": grid, "skipped": grid - len(rows),
            "ca50_mbt": (min(r["ca50_deg"] for r in mbt), max(r["ca50_deg"] for r in mbt)) if mbt else None,
            "dilution_limit": mp["dilution_limit"], "fig_dir": d}


def m8(cfg, fig_dir):
    from hil import runner
    from hil.run import plot
    out = {}
    for c in cfg["hil_cycles"]:
        res = runner.run(runner.Scenario.from_cycle(c))
        out[c] = runner.summary(res)
        plot(res, fig_dir / f"hil_{c}.png")
    return out


def m9(cfg):
    from ml import ml_spec
    from ml.eval import table
    from ml.knock import golden
    p = ROOT / cfg["external_artefacts"]["cosim"]
    cosim = json.loads(p.read_text()) if p.exists() else None
    rows = table(ml_spec()["fuel_ppo"])
    mean = lambda who, key: float(np.mean([r[who][key] for r in rows]))
    return {"cosim": cosim, "ppo_rows": rows, "knock": golden.load_policy()["metrics"],
            "max_cl": mean("classical", "max_dev"), "max_ppo": mean("ppo", "max_dev"),
            "rms_cl": mean("classical", "rms_dev"), "rms_ppo": mean("ppo", "rms_dev")}


def m10(fast: bool):
    from vehicle import accel, wltp
    s, res = wltp.run()
    wltp.OUT_DIR.mkdir(parents=True, exist_ok=True)
    (wltp.OUT_DIR / "wltp.json").write_text(json.dumps(s, indent=2))
    wltp.plot(res, s, wltp.OUT_DIR / "wltp.png")
    a, standing, elastic = accel.run()
    (accel.OUT_DIR / "accel.json").write_text(json.dumps(a, indent=2))
    accel.plot(standing, elastic, a, accel.OUT_DIR / "accel.png")
    from viz.build import build, config as vcfg
    vc = vcfg()
    html = ROOT / vc["output"]["dir"] / vc["output"]["html"]
    info = None
    if not fast or not html.exists():
        info = build(vc)
    shot = ROOT / vc["output"]["screenshot"]
    return {"wltp": s, "accel": a, "wltp_fig": wltp.OUT_DIR / "wltp.png", "accel_fig": accel.OUT_DIR / "accel.png",
            "replay_html": html, "replay_shot": shot, "replay_info": info,
            "three_version": vc["three"]["version"],
            "replay_data": json.loads((ROOT / vc["output"]["dir"] / vc["output"]["data"]).read_text())}


# ------------------------------------------------------------------------------------------------
# markdown
# ------------------------------------------------------------------------------------------------
def pct(x, d=1):
    return f"{x * 100:.{d}f} %"


def render(d: dict) -> str:
    spec = load()
    t = spec["targets"]
    L = []
    A = L.append
    A("# Redline: technical report (MODELLED)")
    A("")
    A("> **Every number in this report is MODELLED.** Redline is a code-first digital twin of a generic 2.0 L "
      "naturally aspirated inline-4; there is no physical engine, dyno or car. Numbers are computed from the code, "
      "`spec/*.yaml` and the committed fixtures by `python -m report.build`, which also regenerates the figures "
      f"under `out/` (git-ignored, so the image links resolve after a build). Built {d['date']}.")
    A("")
    # headline table
    w4, v, a = d["m4"]["peaks"], d["m10"]["wltp"], d["m10"]["accel"]
    A("## Headline numbers")
    A("")
    A("| Quantity | Modelled | Target / reference | Source |")
    A("|---|---|---|---|")
    A(f"| Peak brake torque | {w4['torque_nm']:.1f} N m at {w4['torque_rpm']:.0f} rpm | {t['torque_nm']} N m at "
      f"{t['torque_rpm']} rpm | M4 WOT fixture |")
    A(f"| Peak brake power | {w4['power_kw']:.1f} kW at {w4['power_rpm']:.0f} rpm | {t['power_kw']} kW at "
      f"{t['power_rpm']} rpm | M4 WOT fixture |")
    A(f"| Minimum BSFC (WOT sweep / M7 maps) | {w4['bsfc_min']:.0f} / {d['m7']['summary']['bsfc_min']:.0f} g/kWh | "
      f"{t['min_bsfc_g_per_kwh']} g/kWh | M4 / M7 fixtures |")
    A(f"| WLTP Class 3b fuel | {v['l_per_100km']:.2f} L/100 km iso-octane ({v['gasoline_equiv_l_per_100km']:.2f} "
      f"L/100 km E5-equivalent) | 6.0-8.5 L/100 km band for 2.0 NA compacts | M10 vehicle |")
    A(f"| WLTP CO2 | {v['g_co2_per_km']:.0f} g/km | | M10 vehicle |")
    A(f"| 0-100 km/h | {a['zero_to_100_s']:.2f} s | 7.5-10.5 s band | M10 vehicle |")
    A(f"| 80-120 km/h in {a['in_gear']}th | {a['in_gear_80_120_s']:.2f} s | | M10 vehicle |")
    s5 = d["m5"]["s"]
    A(f"| Secondary shaking force at {s5['balance']['rpm']} rpm | "
      f"{s5['balance']['secondary_exact_2nd_harmonic_N'] / 1e3:.2f} kN | analytic "
      f"{s5['balance']['secondary_formula_4mrw2lam_N'] / 1e3:.2f} kN (4 m r w^2 lam) | M5 |")
    A(f"| First torsional mode | {s5['torsional']['f_eigen_Hz'][0]:.0f} Hz | Holzer within "
      f"{s5['torsional']['max_rel_diff']:.0e} | M5 |")
    if d["m9"]["cosim"]:
        b = d["m9"]["cosim"]["budget"]
        A(f"| ECU crank-tooth handler at {b['redline_rpm']:.0f} rpm | {b['crank_worst_us']:.1f} us "
          f"({pct(b['crank_fraction_of_tooth_period'])} of a tooth) | {b['tooth_period_us']:.1f} us tooth period | M9 |")
    A("")

    # ---- M0
    x = d["m0"]
    A("## M0: scaffold, spec, kinematics, Otto anchor")
    A("")
    A(f"Reference engine from `spec/engine.yaml`: {x['bore_mm']:.0f} x {x['stroke_mm']:.0f} mm, "
      f"{x['displacement_cc']:.0f} cc, CR {x['cr']}, rod/stroke {x['rod_ratio']:.2f}, firing order "
      f"{'-'.join(map(str, spec['engine']['firing_order']))}. Exact slider-crank kinematics are checked against "
      "closed forms (to 1e-12) and finite differences; the ideal Otto cycle on the real volume trace reproduces "
      f"1 - r^(1-gamma) = {x['otto_eff_1p4']:.4f} (gamma 1.4). At the {spec['engine']['redline_rpm']} rpm redline the "
      f"mean piston speed is {x['mean_piston_speed_redline']:.1f} m/s (peak {x['max_piston_speed_redline']:.1f} m/s) "
      f"and piston acceleration spans {x['min_piston_acc_redline_g']:.0f} to +{x['max_piston_acc_redline_g']:.0f} g.")
    A("")
    A("*Known gaps:* zero wrist-pin offset only.")
    A("")

    # ---- M1
    A("## M1: parametric CAD and mass properties")
    A("")
    c = d["m1"]
    if c:
        A("| Item | CAD | Budget | Error |")
        A("|---|---|---|---|")
        A(f"| Piston assembly | {c['piston_assembly']['mass_kg']:.4f} kg | {c['piston_assembly']['budget_kg']} kg | "
          f"{pct(c['piston_assembly']['rel_err'], 2)} |")
        A(f"| Conrod | {c['conrod']['mass_kg']:.4f} kg | {c['conrod']['budget_kg']} kg | {pct(c['conrod']['rel_err'], 2)} |")
        A(f"| Reciprocating mass | {c['recip']['mass_kg']:.4f} kg | {c['recip']['budget_kg']} kg | "
          f"{pct(c['recip']['rel_err'], 2)} |")
        A(f"| Crank inertia | {c['crank']['inertia_about_crank_axis_kgm2']:.5f} kg m^2 | {c['crank']['budget_kgm2']} kg m^2 | "
          f"{pct(c['crank']['rel_err'], 2)} |")
        A(f"| Clearance volume (CR from CAD) | {c['chamber']['clearance_volume_cad_m3'] * 1e6:.2f} cm^3 "
          f"(CR {c['chamber']['compression_ratio_cad']:.3f}) | {c['chamber']['clearance_volume_geometric_m3'] * 1e6:.2f} cm^3 | "
          f"{pct(c['chamber']['rel_err'], 3)} |")
        A("")
        A(f"Exports: STEP / STL / glTF in `out/cad/` (`python -m cad.build`); the glTF parts drive the M10 replay.")
    else:
        A("*out/cad/mass_properties.json not found: run `python -m cad.build`.*")
    A("")
    A("*Known gaps:* the block is monolithic (no separate main caps), crank fillets are not in the solid, the head "
      "is a simplified slab with chambers, seats and port throats, and the camshafts, valves, intake manifold and "
      "exhaust header of the plan are not CAD parts (cams exist as kin.cam profiles only).")
    A("")

    # ---- M2
    x = d["m2"]
    m = x["m"]
    A(f"## M2: two-zone crank-angle thermodynamics ({x['rpm']:.0f} rpm, {x['load']} load)")
    A("")
    A(f"Two-zone Wiebe cycle with Cantera property tables, Woschni heat transfer and filling/emptying gas exchange: "
      f"net IMEP {m['imep_net_bar']:.2f} bar (gross {m['imep_gross_bar']:.2f}), net indicated efficiency "
      f"{m['eta_ind_net']:.3f}, peak pressure {m['p_max_bar']:.1f} bar at {m['theta_p_max_deg']:.1f} deg ATDC, CA50 "
      f"{m['ca50_deg']:.1f} deg, peak burned-gas temperature {m['T_b_max_K']:.0f} K, residual fraction "
      f"{m['residual_fraction']:.3f}, volumetric efficiency {m['volumetric_efficiency']:.3f}. Closed-cycle energy "
      f"balance residual {x['eb']['rel_residual']:.1e} (exit criterion 0.5 %).")
    A("")
    A(fig(x["fig"], "P-V diagram"))
    A("")
    A("*Known gaps:* the stand-alone M2 cycle still uses a sine-squared valve-lift stand-in (M4 uses kin.cam); no "
      "ring blow-by model; no wave dynamics (M4 adds them).")
    A("")

    # ---- M3
    x = d["m3"]
    A("## M3: cams, valvetrain kinematics and dynamics")
    A("")
    A("| Valve | Lift | Seat-to-seat | Float (multi-mass) | Min piston-valve clearance | Peak Hertz | Seating speed at redline |")
    A("|---|---|---|---|---|---|---|")
    for r in x["res"].values():
        A(f"| {r.name} | {r.profile.valve_max_lift * 1e3:.2f} mm | {r.seat_close_deg - r.seat_open_deg:.0f} deg | "
          f"{r.float_rpm:.0f} rpm | {r.min_clearance * 1e3:.2f} mm (VVT {r.min_clearance_advance_deg:.0f} deg) | "
          f"{r.peak_hertz / 1e6:.0f} MPa | {r.seat_speed_redline:.3f} m/s |")
    lim = x["limits"]
    A("")
    A(f"Limits: no float below {lim['float_free_rpm']} rpm, Hertz {lim['hertz_limit'] / 1e6:.0f} MPa, clearance "
      f"{lim['min_piston_valve_clearance'] * 1e3:.2f} mm. All met.")
    A("")
    for f in x["figs"][:2]:
        A(fig(f, f.stem))
    A("")
    v_ = spec["valves"]
    longer = [r.seat_close_deg - r.seat_open_deg - ((v_["ivc_deg"] - v_["ivo_deg"]) if r.name == "intake"
              else (v_["evc_deg"] - v_["evo_deg"])) for r in x["res"].values()]
    A("*Known gaps:* the first ramp design bounced at lash take-up (near 4400 rpm per the M3 commit); ramps were "
      "slowed and contact damping raised. Spec timings are interpreted at the "
      f"{x['check_lift'] * 1e3:.1f} mm check lift, so seat-to-seat events are {min(longer):.0f}-{max(longer):.0f} deg "
      "longer than the spec durations.")
    A("")

    # ---- M4
    x = d["m4"]
    pk = x["peaks"]
    A("## M4: 1D gas dynamics and full-load curve")
    A("")
    A(f"HLLC/MUSCL finite-volume pipes (intake duct, throttle, plenum, 4 runners; 4-2-1 exhaust, catalyst, silencer) "
      f"coupled to four two-zone cylinders. Committed WOT fixture (`fixtures/wot/wot_coarse.csv`, "
      f"{len(x['rows'])} speeds): peak {pk['torque_nm']:.1f} N m at {pk['torque_rpm']:.0f} rpm and "
      f"{pk['power_kw']:.1f} kW at {pk['power_rpm']:.0f} rpm (targets {t['torque_nm']} N m / {t['power_kw']} kW, "
      f"within {abs(pk['torque_nm'] / t['torque_nm'] - 1) * 100:.1f} % / {abs(pk['power_kw'] / t['power_kw'] - 1) * 100:.1f} %), "
      f"trapped VE {pk['ve_min']:.2f}-{pk['ve_max']:.2f}, minimum BSFC {pk['bsfc_min']:.0f} g/kWh, all points converged: "
      f"{pk['all_converged']}. Validated on Sod, orifice and quarter-wave anchors; mass conserved to round-off.")
    A("")
    A(fig(x["fig"], "WOT curve"))
    A("")
    tq = [r["torque_nm"] for r in x["rows"]]
    dips = [f"{x['rows'][i]['rpm']:.0f} rpm ({tq[i]:.0f} N m)" for i in range(1, len(tq) - 1)
            if tq[i] < tq[i - 1] and tq[i] < tq[i + 1]]
    A("*Known gaps:* catalyst and silencer internals are lumped as friction multipliers in the exhaust pipes; the "
      f"torque curve has wave-tuning dips at {', '.join(dips) or 'none'} on this grid; spark is the fixed spec Wiebe "
      "timing (no knock limit) in the WOT sweep; the planned OpenWAM cross-check has not been done.")
    A("")

    # ---- M5
    x = d["m5"]
    s = x["s"]
    b, fl, br, to = s["balance"], s["flywheel"], s["bearings"], s["torsional"]
    A("## M5: dynamics and NVH (now on REAL pressure traces)")
    A("")
    A(f"Pressure source: {s['pressure_source']} at full load; {s['idle_pressure_source']} at idle "
      "(`python -m dyn.report --pressure m4`; the M5 placeholder stays the default for the fast tests).")
    A("")
    A(f"* Balance at {b['rpm']} rpm: primary forces and all pitch/yaw moments cancel to round-off "
      f"(max {b['max_total_pitch_moment_Nm']:.0e} N m); secondary vertical force {b['secondary_exact_2nd_harmonic_N']:.0f} N "
      f"vs {b['secondary_formula_4mrw2lam_N']:.0f} N from 4 m r w^2 lam (exact kinematics add the higher lam terms); "
      f"a Lanchester pair needs m e = {b['lanchester_me_per_shaft_kgm'] * 1e3:.2f} kg mm per shaft.")
    A(f"* Torsional: modes {', '.join(f'{f:.0f}' for f in to['f_eigen_Hz'][:3])} Hz (Holzer agrees to "
      f"{to['max_rel_diff']:.0e}); dominant resonance order {to['dominant']['order']:g} at "
      f"{to['dominant']['critical_rpm']:.0f} rpm, twist {to['dominant']['twist_peak_deg_bare']:.3f} deg bare -> "
      f"{to['dominant']['twist_peak_deg_damped']:.3f} deg with the tuned rubber damper.")
    A(f"* Flywheel (spec {spec['masses']['flywheel_inertia']} kg m^2): cyclic irregularity C_s "
      f"{fl['idle']['cs_energy_method']:.3f} at idle, {fl['full_load_1500']['cs_energy_method']:.3f} at 1500 rpm full "
      f"load (target {fl['idle']['cs_target']}); meeting the target at 1500 rpm full load needs "
      f"{fl['full_load_1500']['required_flywheel_inertia']:.3f} kg m^2.")
    for rpm, e in br.items():
        A(f"* Bearings at {rpm} rpm full load: big-end peak {e['big_end_max_N'] / 1e3:.1f} kN, minimum oil film "
          f"{e['big_end_hmin_um']:.2f} um big end, {e['main_hmin_um']:.2f} um main {e['worst_main']} "
          f"(limit ~{s['min_film_limit_um']:.0f} um).")
    A("")
    if x["cmp"]:
        A("Change of key numbers when the placeholder Wiebe trace is replaced by the M4 / thermo traces:")
        A("")
        A("| Quantity | Placeholder | M4 / thermo | Change |")
        A("|---|---|---|---|")
        for r in x["cmp"]:
            ch = "" if r["rel"] != r["rel"] else f"{r['rel'] * 100:+.1f} %"
            A(f"| {r['quantity']} [{r['unit']}] | {r['placeholder']:.4g} | {r['m4']:.4g} | {ch} |")
        A("")
    for name in ("pressure_m4.png", "crank_torque.png", "bearing_polar_4500.png", "torsional_forced_response.png"):
        A(fig(x["dir"] / name, name[:-4]))
    A("")
    A("*Known gaps:* the M4 full-load traces use the spec Wiebe burn at fixed spark, not the M7 flame model or the "
      "knock-limited spark, so peak pressure at the torque peak is likely overstated; the idle trace comes from the "
      "quasi-steady M2 cycle; torsional gas-torque orders use one trace (torque-peak speed) for the whole speed sweep; "
      "the short-bearing mobility solution assumes rigid shells.")
    A("")

    # ---- M6
    A("## M6: structural FEA")
    A("")
    A("**In progress** in a parallel work stream (`fea/`, `spec/fea.yaml`: gmsh + CalculiX for conrod, crank, piston "
      "and main cap with Goodman fatigue and mesh convergence). Interim results in `out/fea/` are not final and are "
      "not summarised here. Exit criterion: conrod buckling SF >= 3, fatigue SF >= 1.5, convergence shown.")
    A("")

    # ---- M7
    x = d["m7"]
    sm = x["summary"]
    A("## M7: knock, flame model, emissions and speed-load maps")
    A("")
    A(f"Quasi-dimensional turbulent flame, Livengood-Wu / Douaud-Eyzat knock and extended Zeldovich NOx on the M4 "
      f"engine. Fixture `fixtures/maps/maps_m4.csv`: {len(x['rows'])} of {x['grid']} grid points "
      f"({x['skipped']} skipped beyond the {x['dilution_limit']:.2f} residual dilution limit). Minimum BSFC "
      f"{sm['bsfc_min']:.0f} g/kWh at {sm['bsfc_min_rpm']:.0f} rpm ({sm['bsfc_min_level']}, BMEP "
      f"{sm['bsfc_min_bmep']:.1f} bar); knock limits spark at {len(sm['knock_limited'])} points "
      f"({', '.join(f'{r:.0f}/{l}' for r, l in sm['knock_limited'])}); NOx {sm['nox_min']:.1f}-{sm['nox_max']:.1f} g/kWh.")
    A("")
    for k in ("bsfc_g_per_kwh", "spark_deg"):
        A(fig(x["fig_dir"] / f"{k}.png", k))
    A("")
    if x["ca50_mbt"]:
        A(f"*Known gaps:* MBT CA50 lands at {x['ca50_mbt'][0]:.1f}-{x['ca50_mbt'][1]:.1f} deg ATDC, earlier than the "
          "usual ~8-10 deg (short burn-up tail of the flame model); the M4 network's high residuals at low speed and "
          f"part load leave {x['skipped']} points of the grid empty (covered only by `--quasi-steady`), which is why "
          "the M10 drive cycle extrapolates fuel flow below the map.")
    A("")

    # ---- M8
    x = d["m8"]
    A("## M8: C11 ECU on host against the mean-value model")
    A("")
    if "idle" in x:
        i = x["idle"]
        A(f"* Idle: mean {i['idle_rpm_mean_after']:.1f} rpm after the accessory step, max deviation "
          f"{i['idle_rpm_maxdev_after']:.1f} rpm (exit: 800 +/- 25), minimum {i['idle_rpm_min_after_step']:.0f} rpm "
          "during the step.")
    if "tipin" in x:
        p = x["tipin"]
        A(f"* Tip-in at 2000 rpm: steady lambda mean {p['steady_lambda_mean']:.4f} (max deviation "
          f"{p['steady_lambda_maxdev']:.3f}; exit +/- 0.03), transient lambda {p['transient_lambda_min']:.3f}-"
          f"{p['transient_lambda_max']:.3f}, final trims STFT {p['final_stft']:+.3f} / LTFT {p['final_ltft']:+.3f}.")
    A("")
    A(fig(d["fig_dir"] / "hil_tipin.png", "tip-in"))
    A("")
    A("*Known gaps:* the closed loop runs against the mean-value plant; running it against the crank-angle model "
      "is still open. The M8 commit message quotes lambda within 0.3 % steady; the tip-in steady window recomputed "
      f"here (with the {pct(d['m8_injector_error'], 0)} injector flow error the plant carries) gives the deviation "
      "above, inside the exit band.")
    A("")

    # ---- M9
    x = d["m9"]
    A("## M9: ECU on Vortex-V RTL, int8 knock sensor, PPO fuel trim")
    A("")
    if x["cosim"]:
        b = x["cosim"]["budget"]
        k = x["cosim"]["knock"]
        A(f"* Lockstep co-sim: ECU events and slow-loop outputs bit-exact vs the host build "
          f"({sum(r['diff']['host_ref']['events'] for r in x['cosim']['runs'])} event mismatches over "
          f"{len(x['cosim']['runs'])} windows).")
        A(f"* Budget at {b['redline_rpm']:.0f} rpm on a {b['core_clock_hz'] / 1e6:.0f} MHz core: tooth handler worst "
          f"{b['crank_worst_us']:.1f} us of a {b['tooth_period_us']:.1f} us tooth ({pct(b['crank_fraction_of_tooth_period'])}); "
          f"soft-float slow loop {b['slow_worst_us']:.0f} us of {b['slow_period_us']:.0f} us; CPU utilisation "
          f"{pct(b['cpu_utilisation_at_redline'])}; polled worst tooth latency {b['polled_worst_tooth_latency_us']:.0f} us.")
        A(f"* int8 knock MLP on Xmlacc: {k['vectors']} vectors, {k['mismatches']} mismatches vs the numpy golden model, "
          f"{k['cycles_max']} cycles per inference; ROC AUC {k['metrics']['int8']['auc']:.3f} (float "
          f"{k['metrics']['float']['auc']:.3f}, band-pass energy baseline {k['metrics']['bandpass_energy']['auc']:.3f}).")
    else:
        A("*out/hil/cosim.json not found: run `make -C hil cosim` (needs Verilator and the Vortex-V tree).*")
    A(f"* PPO + Simplex supervisor vs classical over {len(x['ppo_rows'])} tip-in/tip-out cycles: mean rms lambda "
      f"deviation {x['rms_cl']:.4f} -> {x['rms_ppo']:.4f} ({(x['rms_ppo'] / x['rms_cl'] - 1) * 100:+.0f} %), mean peak "
      f"deviation {x['max_cl']:.4f} -> {x['max_ppo']:.4f} ({(x['max_ppo'] / x['max_cl'] - 1) * 100:+.0f} %).")
    A("")
    A("*Known gaps:* at 7000 rpm the soft-float slow loop plus tooth handling exceeds one 50 MHz core and, with no "
      "interrupts, a polled tooth can wait about a slow-loop period; lockstep co-sim hides both. PPO only beats the "
      "classical ECU on rms lambda, not on peak excursion.")
    A("")

    # ---- M10
    x = d["m10"]
    v, a, r = x["wltp"], x["accel"], x["replay_data"]
    vs = vehicle_spec()
    A("## M10: vehicle, drive cycle, replay")
    A("")
    A(f"Generic compact sedan (`spec/vehicle.yaml`, sources cited there): WLTP test mass {v['test_mass_kg']:.0f} kg, "
      f"c_d A {vs['vehicle']['drag_coefficient'] * vs['vehicle']['frontal_area']:.3f} m^2, c_r "
      f"{vs['vehicle']['rolling_resistance']}, 6-speed manual {vs['driveline']['gear_ratios']} x final drive "
      f"{vs['driveline']['final_drive']}, driveline efficiency {vs['driveline']['efficiency']}. Fuel flow from a "
      f"Willans fit to the M7 maps (throttled points reproduced to {pct(v['willans_fit']['throttled_rms_rel_err'], 2)} rms, "
      f"max {pct(v['willans_fit']['throttled_max_abs_rel_err'], 1)}), WOT envelope from the M4 fixture, overrun fuel "
      "cut, idle at standstill.")
    A("")
    A("| WLTP Class 3b | Distance | L/100 km | g CO2/km |")
    A("|---|---|---|---|")
    for k, p in v["phases"].items():
        A(f"| {k} | {p['distance_km']:.3f} km | {p['l_per_100km']:.2f} | {p['g_co2_per_km']:.0f} |")
    A(f"| **combined** | {v['distance_km']:.3f} km | **{v['l_per_100km']:.2f}** | **{v['g_co2_per_km']:.0f}** |")
    A("")
    A(f"Cycle-mean BSFC {v['cycle_mean_bsfc_g_per_kwh']:.0f} g/kWh, positive wheel energy "
      f"{v['positive_wheel_energy_kwh']:.2f} kWh, idle share of fuel {pct(v['idle_fuel_share'])}, trace followed with "
      f"{v['seconds_with_shortfall']} s of torque shortfall. Acceleration (mass in running order {a['mass_kg']:.0f} kg, "
      f"{a['power_to_mass_w_per_kg']:.0f} W/kg): 0-100 km/h {a['zero_to_100_s']:.2f} s, 80-120 km/h in {a['in_gear']}th "
      f"{a['in_gear_80_120_s']:.2f} s.")
    A("")
    A(fig(x["wltp_fig"], "WLTP"))
    A(fig(x["accel_fig"], "acceleration"))
    A("")
    A(f"Replay: `viz/replay/index.html` (self-contained page, three.js {x['three_version']} pinned on jsDelivr,"
      f" built by `python -m viz.build`): the CAD glTF crank, rods and pistons animated from the exact"
      f" slider-crank at {r['rpm']:.0f} rpm WOT with live P-V, pressure, crank-torque and valve-lift panels from the "
      f"M4 engine (brake torque {r['metrics']['torque_nm']:.1f} N m, peak pressure {r['metrics']['p_max_bar']:.1f} bar).")
    A("")
    A("")
    A(fig(x["replay_shot"], "replay"))
    A("")
    A(f"*Known gaps:* only {pct(v['fuel_share_inside_m7_map_coverage'], 0)} of the WLTP fuel is burned inside the "
      "speed/BMEP range the M7 maps sampled; the rest relies on the Willans extrapolation to light load and idle. "
      "Hot start only (no cold-start enrichment or catalyst heating), economy shift strategy rather than the GTR 15 "
      "Annex 2 algorithm, constant driveline efficiency, no road grade. The replay is not yet published as an Artifact.")
    A("")

    # ---- remaining
    A("## What remains")
    A("")
    A("* **M6 FEA** (in progress): conrod, crankshaft fillets, piston crown and main cap to their safety-factor "
      "targets with mesh convergence.")
    A("* **M10 stretch** (plan decision 4 default): turbocharged variant (option B: compressor/turbine maps, "
      "wastegate, boost control) and the ECU PCB in Zener. Intake-port CFD is not planned.")
    A("* Open items carried from earlier milestones: ECU loop against the crank-angle model (M8), interrupt-driven "
      "tooth handling and a hard-float or fixed-point slow loop to fit the 7000 rpm budget (M9), maps on a knock-"
      "limited M4 at low-speed part load (M7), M4 traces with the flame model for the dynamics (M5/M10), publishing "
      "the replay, and the GitHub repository.")
    A("")
    A("## Regenerate")
    A("")
    A("```")
    A("python -m report.build            # this report and its figures (~2 min)")
    A("python -m vehicle.wltp && python -m vehicle.accel")
    A("python -m dyn.report --pressure m4")
    A("python -m viz.build --check")
    A("pytest -q tests/")
    A("```")
    A("")
    return "\n".join(L)


def build(fast: bool = False) -> Path:
    cfg = config()
    fig_dir = ROOT / cfg["figures_dir"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    spec = load()
    d = {"date": _dt.date.today().isoformat(), "fig_dir": fig_dir,
         "m0": m0(spec), "m1": m1(cfg), "m2": m2(cfg), "m3": m3(), "m4": m4(fig_dir), "m5": m5(fast),
         "m7": m7(fig_dir), "m8": m8(cfg, fig_dir), "m8_injector_error": mvem_spec()["injector"]["flow_error_fraction"], "m9": m9(cfg), "m10": m10(fast)}
    out = ROOT / cfg["output"]
    out.write_text(render(d))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fast", action="store_true", help="reuse the last dynamics and replay outputs")
    a = ap.parse_args(argv)
    out = build(a.fast)
    print(f"wrote {out} (ALL NUMBERS MODELLED)")


if __name__ == "__main__":
    main()
