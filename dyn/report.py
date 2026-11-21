"""M5 dynamics/NVH report: `python -m dyn.report [--pressure placeholder|m4]` prints a summary and writes
plots to out/dyn/ (placeholder) or out/dyn/m4/ (real traces).

All numbers are MODELLED (no physical engine). `--pressure placeholder` (default, fast, used by the M5
tests) takes cylinder pressure from the single-zone Wiebe stand-in in dyn/pressure.py; `--pressure m4`
(M10) uses REAL traces: the M4 1D engine at full load at each speed the report needs, and the thermo/
cycle solved for zero brake torque at idle (dyn.pressure.PressureProvider). With m4 the report also
prints and writes the change of the key numbers against the placeholder (out/dyn/m4/compare.json).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from dyn import balance, bearings, flywheel, forces, params, torsional
from dyn.pressure import PressureProvider

OUT_DIR = params.rspec.ROOT / "out" / "dyn"
GRID = 2880                    # points per 720-deg cycle (0.25 deg) for torque/load work
REV_GRID = 4096                # points per revolution for balance harmonics
DEG = 180.0 / math.pi

# reference categorical palette (dataviz skill, light mode), fixed order
C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, MUTED = "#0b0b0b", "#52514e"


def compute(spec=None, dspec=None, with_orbits: bool = True, pressure: str | PressureProvider = "placeholder") -> dict:
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    eng = spec["engine"]
    g = params.geometry(spec)
    m = params.masses(spec)
    th = forces.cycle_grid(GRID)
    prov = pressure if isinstance(pressure, PressureProvider) else PressureProvider(pressure, spec, dspec)
    torque_rpm = spec["targets"]["torque_rpm"]
    # the full-load trace at the torque-peak speed stands for "full load" where one trace is needed
    # (torsional gas-torque orders); with the placeholder every speed shares the same trace
    traces = {"full_load": prov("full_load", torque_rpm), "idle": prov("idle", eng["idle_rpm"])}
    res: dict = {"modelled": True, "pressure": prov.kind, "pressure_source": traces["full_load"].source,
                 "idle_pressure_source": traces["idle"].source}

    # (1) forces and torque
    tr = traces["full_load"]
    torque = {}
    used = {}
    for label, op, rpm in [("idle", "idle", eng["idle_rpm"]),
                           ("full_load_torque_rpm", "full_load", torque_rpm),
                           ("full_load_redline", "full_load", eng["redline_rpm"])]:
        w = params.omega_from_rpm(rpm)
        trace = prov(op, rpm)
        used[label] = trace
        T_cyl, T = forces.multi_cylinder_torque(th, trace, w, spec, dspec)
        f1 = forces.cylinder_forces(th, trace, w, spec, dspec)
        torque[label] = {
            "rpm": rpm, "mean_torque_Nm": float(T.mean()),
            "imep_identity_Nm": forces.mean_torque_identity(trace, spec),
            "imep_bar": forces.imep(trace, spec) / 1e5, "peak_pressure_bar": trace.peak / 1e5,
            "peak_pressure_deg": trace.peak_angle_deg,
            "peak_Nm": float(T.max()), "min_Nm": float(T.min()),
            "inertia_torque_mean_Nm": float(f1.T_inertia.mean()),
            "peak_rod_force_N": float(f1.F_rod.max()), "peak_side_thrust_N": float(np.abs(f1.F_side).max()),
            "_T": T, "_T1": T_cyl[0], "_f1": f1,
        }
    res["torque"] = torque
    res["peak_pressure_bar"] = {op: t.peak / 1e5 for op, t in traces.items()}
    res["_traces"] = used

    # (2) balance
    thr = np.linspace(0.0, params.REV, REV_GRID, endpoint=False)
    w7 = params.omega_from_rpm(eng["redline_rpm"])
    sh = balance.shaking(thr, w7, spec, dspec)
    sec_an = balance.secondary_force_analytic(thr, w7, spec)
    pair = balance.lanchester_sizing(spec)
    fz_b, fy_b = balance.lanchester_force(thr, w7, pair)
    res["balance"] = {
        "rpm": eng["redline_rpm"],
        "m_rec_kg": m.reciprocating,
        "secondary_exact_2nd_harmonic_N": abs(balance.harmonic(sh.Fz, thr, 2)),
        "secondary_formula_4mrw2lam_N": abs(balance.harmonic(sec_an, thr, 2)),
        "fourth_order_N": abs(balance.harmonic(sh.Fz, thr, 4)),
        "max_primary_force_N": float(np.abs(sh.Fz_primary).max()),
        "max_primary_moment_Nm": float(np.abs(sh.My_primary).max()),
        "max_secondary_moment_Nm": float(np.abs(sh.My_secondary).max()),
        "max_total_pitch_moment_Nm": float(np.abs(sh.My).max()),
        "max_total_yaw_moment_Nm": float(np.abs(sh.Mz).max()),
        "max_rotating_force_N": float(np.hypot(sh.Fz_rot, sh.Fy).max()),
        "lanchester_me_per_shaft_kgm": pair.me_per_shaft,
        "residual_2nd_with_shafts_N": abs(balance.harmonic(sh.Fz + fz_b, thr, 2)),
        "_th": thr, "_sh": sh, "_sec": sec_an, "_fzb": fz_b,
    }

    # (3) flywheel
    fw = dspec["flywheel"]
    cases = [("idle", "idle", eng["idle_rpm"], fw["cs_target_idle"])] + [
        (f"full_load_{r}", "full_load", r, fw["cs_target_full_load"]) for r in fw["full_load_rpm"]]
    sweep = np.linspace(fw["inertia_sweep"]["min"], fw["inertia_sweep"]["max"], fw["inertia_sweep"]["points"])
    fly = {}
    for label, op, rpm, target in cases:
        w = params.omega_from_rpm(rpm)
        _, T = forces.multi_cylinder_torque(th, prov(op, rpm), w, spec, dspec)
        dE = flywheel.fluctuation_energy(th, T)
        I = flywheel.rotating_inertia(None, spec, dspec)
        _, ws = flywheel.simulate_speed(th, T, I, w)
        fly[label] = {
            "rpm": rpm, "dE_J": dE, "system_inertia": I,
            "cs_energy_method": flywheel.cyclic_irregularity(dE, I, w),
            "cs_time_integration": float((ws.max() - ws.min()) / ws.mean()),
            "cs_target": target,
            "required_flywheel_inertia": flywheel.required_flywheel_inertia(dE, w, target, spec, dspec),
            "_sweep": [flywheel.cyclic_irregularity(dE, flywheel.rotating_inertia(j, spec, dspec), w) for j in sweep],
        }
    res["flywheel"] = fly
    res["flywheel_sweep"] = sweep

    # (4) bearings
    bs = dspec["bearings"]
    film = bearings.film_for(dspec)
    bres = {}
    for rpm in bs["rpm_cases"]:
        w = params.omega_from_rpm(rpm)
        tr_b = prov("full_load", rpm)
        be = bearings.big_end_loads(th, tr_b, w, spec, dspec)
        ml = bearings.main_bearing_loads(th, tr_b, w, spec, dspec)
        mag_be = np.hypot(*be.F_pin)
        mag_m = np.hypot(ml[:, 0], ml[:, 1])
        worst = int(np.argmax(mag_m.max(1)))
        entry = {
            "big_end_max_N": float(mag_be.max()), "big_end_mean_N": float(mag_be.mean()),
            "main_max_N": mag_m.max(1).tolist(), "main_mean_N": mag_m.mean(1).tolist(),
            "worst_main": worst + 1, "_be": be, "_ml": ml,
        }
        if with_orbits:
            ob = bearings.journal_orbit(th, be.F_pin, w, be.omega_rod, bearings.bearing("big_end", dspec), w,
                                        film, bs["orbit_cycles"], bs["orbit_steps_per_cycle"])
            om = bearings.journal_orbit(th, -ml[worst], w, 0.0, bearings.bearing("main", dspec), w,
                                        film, bs["orbit_cycles"], bs["orbit_steps_per_cycle"])
            entry.update({"big_end_hmin_um": ob.min_film * 1e6, "big_end_hmin_deg": ob.min_film_angle_deg,
                          "main_hmin_um": om.min_film * 1e6, "main_hmin_deg": om.min_film_angle_deg,
                          "_ob": ob, "_om": om})
        bres[rpm] = entry
    res["bearings"] = bres
    res["min_film_limit_um"] = bs["min_film_limit"] * 1e6

    # (5) torsional
    ts = dspec["torsional"]
    model = torsional.build_model(spec, dspec)
    f_eig, modes = torsional.eigen_frequencies(model.J, model.k, return_modes=True)
    f_hol = torsional.holzer_frequencies(model.J, model.k, f_eig[-1] * 1.05)
    orders, Cg, Ci = torsional.torque_orders(tr, ts["max_order"], spec=spec, dspec=dspec)
    rpm = np.linspace(ts["rpm_sweep"]["min"], ts["rpm_sweep"]["max"], ts["rpm_sweep"]["points"])
    damp = torsional.den_hartog_damper(model, ts["damper"]["ring_inertia"], ts["damper"]["tuned_mode"])
    model_d = torsional.with_damper(model, damp)
    f_eig_d = torsional.eigen_frequencies(model_d.J, model_d.k)
    f_hol_d = torsional.holzer_frequencies(model_d.J, model_d.k, f_eig_d[-1] * 1.05)
    r0 = torsional.forced_response(model, rpm, orders, Cg, Ci, spec)
    r1 = torsional.forced_response(model_d, rpm, orders, Cg, Ci, spec)
    tw0 = np.abs(r0[:, :, model.hub] - r0[:, :, -1])
    tw1 = np.abs(r1[:, :, model_d.hub] - r1[:, :, -1])
    f1 = f_eig[1]
    w_ref = params.omega_from_rpm(ts["rpm_sweep"]["max"])
    crit = []
    for i, q in enumerate(orders):
        n_c = 60.0 * f1 / q
        if ts["rpm_sweep"]["min"] <= n_c <= ts["rpm_sweep"]["max"]:
            wc = params.omega_from_rpm(n_c)
            Tq = abs(Cg[i] + wc**2 * Ci[i])
            crit.append({"order": float(q), "critical_rpm": n_c, "torque_harmonic_Nm": Tq,
                         "phase_sum": torsional.phase_sum(modes[:, 1], model.throws, q, spec),
                         "twist_peak_deg_bare": float(tw0[i].max() * DEG),
                         "twist_peak_deg_damped": float(tw1[i].max() * DEG)})
    crit.sort(key=lambda c: -c["twist_peak_deg_bare"])
    dominant = crit[0]
    res["torsional"] = {
        "J": model.J.tolist(), "k": model.k.tolist(),
        "f_eigen_Hz": f_eig[1:].tolist(), "f_holzer_Hz": f_hol[1:len(f_eig)].tolist(),
        "max_rel_diff": float(np.max(np.abs(f_hol[1:len(f_eig)] / f_eig[1:] - 1))),
        "f_eigen_damped_Hz": f_eig_d[1:].tolist(), "f_holzer_damped_Hz": f_hol_d[1:len(f_eig_d)].tolist(),
        "critical_orders": crit[:8], "dominant": dominant,
        "damper": {k: float(v) for k, v in damp.__dict__.items()},
        "twist_sum_peak_deg_bare": float(tw0.sum(0).max() * DEG),
        "twist_sum_peak_deg_damped": float(tw1.sum(0).max() * DEG),
        "_rpm": rpm, "_orders": orders, "_Cg": Cg, "_Ci": Ci, "_tw0": tw0, "_tw1": tw1, "_modes": modes,
        "_labels": model.labels, "_w_ref": w_ref,
    }
    return res


# --------------------------------------------------------------------------------------------

def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=INK, fontsize=10, loc="left")
    ax.set_xlabel(xlabel, color=MUTED, fontsize=9)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    ax.grid(True, color="#e6e5e0", linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)


def plots(res: dict, out: Path, spec=None, dspec=None) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    out.mkdir(parents=True, exist_ok=True)
    files = []
    lw = 1.6
    real = res.get("pressure") == "m4"
    ptag = "M4/thermo pressure" if real else "placeholder pressure"
    tag = f"  (modelled; {ptag})"

    def save(fig, name):
        fig.tight_layout()
        p = out / name
        fig.savefig(p, dpi=130)
        plt.close(fig)
        files.append(str(p))

    th_deg = forces.cycle_grid(GRID) * DEG
    # pressure
    fig, ax = plt.subplots(figsize=(7, 3.6))
    seen = set()
    for (label, t), c in zip(res["_traces"].items(), (C2, C1, C3)):
        if id(t) in seen:
            continue
        seen.add(id(t))
        ax.plot(t.theta * DEG, t.p / 1e5, color=c, lw=lw,
                label=f"{label.replace('_', ' ')} ({res['torque'][label]['rpm']} rpm)" if real else label.split("_")[0])
    title = ("REAL cylinder pressure: M4 1D engine (full load), thermo/ zero-torque cycle (idle)" if real
             else "PLACEHOLDER cylinder pressure (single-zone Wiebe)")
    _style(ax, title + tag, "crank angle [deg]", "p [bar]")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "pressure_m4.png" if real else "pressure_placeholder.png")

    # single-cylinder forces
    tq = res["torque"]["full_load_torque_rpm"]
    f1 = tq["_f1"]
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for arr, c, lab in ((f1.F_gas, C1, "gas"), (f1.F_inertia, C2, "inertia"), (f1.F_rod, C3, "conrod"),
                        (f1.F_side, C4, "side thrust"), (f1.F_tangential, C5, "tangential")):
        ax.plot(th_deg, arr / 1e3, color=c, lw=lw, label=lab)
    _style(ax, f"Cylinder 1 forces, full load {tq['rpm']} rpm" + tag, "crank angle [deg]", "force [kN]")
    ax.legend(frameon=False, fontsize=8, ncol=3)
    save(fig, "forces_single_cylinder.png")

    # torque
    fig, axs = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True)
    for ax, key in zip(axs, ("idle", "full_load_torque_rpm", "full_load_redline")):
        d = res["torque"][key]
        ax.plot(th_deg, d["_T1"], color=C2, lw=1.2, label="cylinder 1")
        ax.plot(th_deg, d["_T"], color=C1, lw=lw, label="engine (1-3-4-2)")
        ax.axhline(d["mean_torque_Nm"], color=MUTED, lw=1, ls="--", label="mean")
        _style(ax, f"Crank torque, {key.replace('_', ' ')} ({d['rpm']} rpm)", "", "T [N m]")
    axs[-1].set_xlabel("crank angle [deg]", color=MUTED, fontsize=9)
    axs[0].legend(frameon=False, fontsize=8, ncol=3)
    save(fig, "crank_torque.png")

    # shaking
    b = res["balance"]
    fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    td = b["_th"] * DEG
    axs[0].plot(td, b["_sh"].Fz / 1e3, color=C1, lw=lw, label="exact vertical force")
    axs[0].plot(td, b["_sec"] / 1e3, color=C2, lw=1.2, ls="--", label="4 m r w^2 lam cos 2t")
    axs[0].plot(td, (b["_sh"].Fz + b["_fzb"]) / 1e3, color=C3, lw=1.2, label="with Lanchester pair")
    _style(axs[0], f"Shaking force at {b['rpm']} rpm" + tag, "", "F_z [kN]")
    axs[0].legend(frameon=False, fontsize=8)
    axs[1].plot(td, b["_sh"].My_primary, color=C1, lw=lw, label="primary pitch moment")
    axs[1].plot(td, b["_sh"].My_secondary, color=C2, lw=lw, label="secondary pitch moment")
    axs[1].plot(td, b["_sh"].My, color=C3, lw=1.2, label="total pitch moment")
    # scale: primary moment of one end cylinder alone, so round-off (~1e-13) reads as zero
    g = params.geometry(spec)
    w7 = params.omega_from_rpm(b["rpm"])
    scale = b["m_rec_kg"] * g.crank_radius * w7**2 * np.abs(params.cylinder_positions(spec, dspec)).max()
    axs[1].set_ylim(-scale, scale)
    _style(axs[1], f"Pitch moments (cancel; y-range = one end cylinder's primary moment, max |sum| "
                   f"{b['max_total_pitch_moment_Nm']:.0e} N m)", "crank angle [deg]", "M_y [N m]")
    axs[1].legend(frameon=False, fontsize=8)
    save(fig, "shaking.png")

    # flywheel
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for (k, d), c in zip(res["flywheel"].items(), (C1, C2, C3, C4)):
        ax.plot(res["flywheel_sweep"], d["_sweep"], color=c, lw=lw, label=f"{k.replace('_', ' ')} ({d['rpm']} rpm)")
    ax.axhline(dspec["flywheel"]["cs_target_idle"], color=MUTED, ls="--", lw=1, label="target")
    ax.axvline(spec["masses"]["flywheel_inertia"], color=MUTED, ls=":", lw=1, label="spec flywheel")
    ax.set_yscale("log")
    _style(ax, "Cyclic speed variation vs flywheel inertia" + tag, "flywheel inertia [kg m^2]", "C_s [-]")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "flywheel_cs.png")

    # bearings: polar diagrams
    for rpm, d in res["bearings"].items():
        be, ml = d["_be"], d["_ml"]
        fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))
        axs[0].plot(be.pin_frame[1] / 1e3, be.pin_frame[0] / 1e3, color=C1, lw=lw)
        _style(axs[0], "Big end on crankpin (crankpin frame)", "tangential [kN]", "radial outward [kN]")
        axs[1].plot(be.rod_frame[1] / 1e3, be.rod_frame[0] / 1e3, color=C2, lw=lw)
        _style(axs[1], "Big-end shell load (rod frame)", "transverse [kN]", "axial to small end [kN]")
        w = d["worst_main"] - 1
        axs[2].plot(ml[w, 1] / 1e3, ml[w, 0] / 1e3, color=C3, lw=lw)
        _style(axs[2], f"Main bearing {w + 1} load (engine frame)", "y [kN]", "z [kN]")
        for ax in axs:
            ax.set_aspect("equal", adjustable="datalim")
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.axvline(0, color=MUTED, lw=0.6)
        fig.suptitle(f"Bearing load polars, full load {rpm} rpm (modelled; {ptag})", fontsize=10)
        save(fig, f"bearing_polar_{rpm}.png")
        if "_ob" in d:
            fig, ax = plt.subplots(figsize=(7, 3.6))
            ax.plot(d["_ob"].theta * DEG, d["_ob"].h_min * 1e6, color=C1, lw=lw, label="big end")
            ax.plot(d["_om"].theta * DEG, d["_om"].h_min * 1e6, color=C3, lw=lw, label=f"main {d['worst_main']}")
            ax.axhline(res["min_film_limit_um"], color=MUTED, ls="--", lw=1, label="~1 um limit")
            _style(ax, f"Minimum oil film, short-bearing mobility solution, {rpm} rpm" + tag,
                   "crank angle [deg]", "h_min [um]")
            ax.legend(frameon=False, fontsize=8)
            save(fig, f"oil_film_{rpm}.png")

    # torsional
    t = res["torsional"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for i, c in zip(range(1, 4), (C1, C2, C3)):
        ph = t["_modes"][:, i] / np.max(np.abs(t["_modes"][:, i]))
        axs[0].plot(range(len(ph)), ph, color=c, lw=lw, marker="o", ms=4, label=f"mode {i}: {t['f_eigen_Hz'][i - 1]:.0f} Hz")
    axs[0].set_xticks(range(len(t["_labels"])), t["_labels"], fontsize=8)
    _style(axs[0], "Torsional mode shapes (bare crank)", "", "relative angle")
    axs[0].legend(frameon=False, fontsize=8)
    wr = t["_w_ref"]
    axs[1].bar(t["_orders"] - 0.1, np.abs(t["_Cg"]), width=0.2, color=C1, label="gas (full load)")
    axs[1].bar(t["_orders"] + 0.1, np.abs(t["_Ci"]) * wr**2, width=0.2, color=C2, label="inertia (redline)")
    _style(axs[1], "Single-cylinder torque harmonics by order", "order", "|T_q| [N m]")
    axs[1].legend(frameon=False, fontsize=8)
    save(fig, "torsional_modes_orders.png")

    fig, ax = plt.subplots(figsize=(7, 3.8))
    q = t["dominant"]["order"]
    i = int(round(2 * q)) - 1
    ax.plot(t["_rpm"], t["_tw0"][i] * DEG, color=C1, lw=lw, label=f"order {q:g}, bare")
    ax.plot(t["_rpm"], t["_tw1"][i] * DEG, color=C2, lw=lw, label=f"order {q:g}, tuned damper")
    ax.plot(t["_rpm"], t["_tw0"].sum(0) * DEG, color=C3, lw=1.2, ls="--", label="all orders summed, bare")
    ax.plot(t["_rpm"], t["_tw1"].sum(0) * DEG, color=C4, lw=1.2, ls="--", label="all orders summed, damper")
    _style(ax, "Crank twist amplitude hub to flywheel" + tag, "engine speed [rpm]", "twist [deg]")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "torsional_forced_response.png")
    return files


def _public(obj):
    if isinstance(obj, dict):
        return {str(k): _public(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [_public(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def summary_text(res: dict) -> str:
    L = ["Redline M5 dynamics/NVH report. ALL NUMBERS MODELLED (no physical engine).",
         f"Cylinder pressure ({res.get('pressure', 'placeholder')}): full load {res['pressure_source']};"
         f" idle {res.get('idle_pressure_source', res['pressure_source'])}", ""]
    L.append("(1) Forces and crank torque")
    for k, d in res["torque"].items():
        L.append(f"  {k:22s} {d['rpm']:5d} rpm  IMEP {d['imep_bar']:5.2f} bar  mean T {d['mean_torque_Nm']:7.2f} N m"
                 f" (IMEP*Vd/4pi {d['imep_identity_Nm']:7.2f})  peak {d['peak_Nm']:7.1f}  min {d['min_Nm']:7.1f}"
                 f"  inertia-torque mean {d['inertia_torque_mean_Nm']:.1e}  peak rod {d['peak_rod_force_N'] / 1e3:5.1f} kN"
                 f"  side {d['peak_side_thrust_N'] / 1e3:4.2f} kN")
    L.append(f"  peak pressure: full load {res['peak_pressure_bar']['full_load']:.1f} bar,"
             f" idle {res['peak_pressure_bar']['idle']:.1f} bar")
    b = res["balance"]
    L += ["", f"(2) Shaking at {b['rpm']} rpm (m_rec = {b['m_rec_kg']:.3f} kg)",
          f"  secondary vertical force: exact 2nd harmonic {b['secondary_exact_2nd_harmonic_N']:.0f} N,"
          f" 4 m r w^2 lam {b['secondary_formula_4mrw2lam_N']:.0f} N"
          f" (ratio {b['secondary_exact_2nd_harmonic_N'] / b['secondary_formula_4mrw2lam_N']:.4f}); 4th order"
          f" {b['fourth_order_N']:.0f} N",
          f"  max primary force {b['max_primary_force_N']:.1e} N, primary moment {b['max_primary_moment_Nm']:.1e} N m,"
          f" secondary moment {b['max_secondary_moment_Nm']:.1e} N m, total pitch {b['max_total_pitch_moment_Nm']:.1e},"
          f" yaw {b['max_total_yaw_moment_Nm']:.1e}, rotating force {b['max_rotating_force_N']:.1e} N",
          f"  Lanchester pair at 2x: m e = {b['lanchester_me_per_shaft_kgm'] * 1e3:.3f} kg mm per shaft;"
          f" residual 2nd order {b['residual_2nd_with_shafts_N']:.1e} N"]
    L += ["", "(3) Flywheel (energy-fluctuation method)"]
    for k, d in res["flywheel"].items():
        L.append(f"  {k:16s} {d['rpm']:5d} rpm  dE {d['dE_J']:6.1f} J  I_sys {d['system_inertia']:.4f}  C_s"
                 f" {d['cs_energy_method']:.4f} (time integration {d['cs_time_integration']:.4f})  target"
                 f" {d['cs_target']:.3f} -> flywheel needed {d['required_flywheel_inertia']:.3f} kg m^2")
    L += ["", f"(4) Bearings (short-bearing mobility solution; film limit ~{res['min_film_limit_um']:.1f} um)"]
    for rpm, d in res["bearings"].items():
        s = (f"  {rpm} rpm: big end max {d['big_end_max_N'] / 1e3:.1f} kN mean {d['big_end_mean_N'] / 1e3:.1f} kN;"
             f" mains max {', '.join(f'{x / 1e3:.1f}' for x in d['main_max_N'])} kN"
             f" mean {', '.join(f'{x / 1e3:.1f}' for x in d['main_mean_N'])} kN")
        L.append(s)
        if "big_end_hmin_um" in d:
            L.append(f"      h_min big end {d['big_end_hmin_um']:.2f} um at {d['big_end_hmin_deg']:.0f} deg;"
                     f" main {d['worst_main']} {d['main_hmin_um']:.2f} um at {d['main_hmin_deg']:.0f} deg")
    t = res["torsional"]
    L += ["", "(5) Torsional (bare crank: hub, 4 throws, flywheel)",
          "  eigen  [Hz]: " + ", ".join(f"{f:.1f}" for f in t["f_eigen_Hz"]),
          "  Holzer [Hz]: " + ", ".join(f"{f:.1f}" for f in t["f_holzer_Hz"])
          + f"   max rel diff {t['max_rel_diff']:.1e}",
          "  with damper eigen [Hz]: " + ", ".join(f"{f:.1f}" for f in t["f_eigen_damped_Hz"]),
          "  critical orders of mode 1 in the speed range (sorted by bare resonant twist):"]
    for c in t["critical_orders"]:
        L.append(f"    order {c['order']:4.1f} at {c['critical_rpm']:6.0f} rpm  |T_q| {c['torque_harmonic_Nm']:6.1f} N m"
                 f"  phase sum {c['phase_sum']:.2f}  twist {c['twist_peak_deg_bare']:.3f} deg bare,"
                 f" {c['twist_peak_deg_damped']:.3f} deg damped")
    d = t["damper"]
    L.append(f"  tuned rubber damper: ring {d['ring_inertia']:.4f} kg m^2, mu {d['mass_ratio']:.3f}, tuned to"
             f" {d['tuning_ratio'] * d['target_hz']:.0f} Hz, k {d['stiffness']:.0f} N m/rad,"
             f" c {d['damping']:.2f} N m s/rad (zeta {d['damping_ratio']:.3f})")
    dom = t["dominant"]
    L.append(f"  dominant resonance order {dom['order']:g}: {dom['twist_peak_deg_bare']:.3f} -> "
             f"{dom['twist_peak_deg_damped']:.3f} deg"
             f" ({dom['twist_peak_deg_bare'] / dom['twist_peak_deg_damped']:.1f}x lower);"
             f" all orders summed peak {t['twist_sum_peak_deg_bare']:.3f} -> {t['twist_sum_peak_deg_damped']:.3f} deg")
    return "\n".join(L)


# key numbers compared between pressure sources: (label, path into the public summary, unit)
KEY_NUMBERS = [
    ("peak pressure, full load @ torque rpm", ("torque", "full_load_torque_rpm", "peak_pressure_bar"), "bar"),
    ("peak pressure, full load @ redline", ("torque", "full_load_redline", "peak_pressure_bar"), "bar"),
    ("peak pressure, idle", ("torque", "idle", "peak_pressure_bar"), "bar"),
    ("IMEP, full load @ torque rpm", ("torque", "full_load_torque_rpm", "imep_bar"), "bar"),
    ("IMEP, full load @ redline", ("torque", "full_load_redline", "imep_bar"), "bar"),
    ("IMEP, idle", ("torque", "idle", "imep_bar"), "bar"),
    ("peak crank torque @ torque rpm", ("torque", "full_load_torque_rpm", "peak_Nm"), "N m"),
    ("peak crank torque @ redline", ("torque", "full_load_redline", "peak_Nm"), "N m"),
    ("peak rod force @ torque rpm", ("torque", "full_load_torque_rpm", "peak_rod_force_N"), "N"),
    ("flywheel dE, idle", ("flywheel", "idle", "dE_J"), "J"),
    ("C_s idle", ("flywheel", "idle", "cs_energy_method"), "-"),
    ("C_s full load 1500", ("flywheel", "full_load_1500", "cs_energy_method"), "-"),
    ("flywheel needed, idle", ("flywheel", "idle", "required_flywheel_inertia"), "kg m^2"),
    ("flywheel needed, full load 1500", ("flywheel", "full_load_1500", "required_flywheel_inertia"), "kg m^2"),
    ("big-end max load 4500", ("bearings", "4500", "big_end_max_N"), "N"),
    ("big-end max load 7000", ("bearings", "7000", "big_end_max_N"), "N"),
    ("big-end h_min 4500", ("bearings", "4500", "big_end_hmin_um"), "um"),
    ("big-end h_min 7000", ("bearings", "7000", "big_end_hmin_um"), "um"),
    ("main h_min 4500 (worst main)", ("bearings", "4500", "main_hmin_um"), "um"),
    ("main h_min 7000 (worst main)", ("bearings", "7000", "main_hmin_um"), "um"),
    ("dominant torsional order twist, bare", ("torsional", "dominant", "twist_peak_deg_bare"), "deg"),
    ("dominant torsional order twist, damped", ("torsional", "dominant", "twist_peak_deg_damped"), "deg"),
    ("all-order twist peak, bare", ("torsional", "twist_sum_peak_deg_bare"), "deg"),
    ("all-order twist peak, damped", ("torsional", "twist_sum_peak_deg_damped"), "deg"),
]


def _dig(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def compare(new: dict, ref: dict) -> list[dict]:
    """Key numbers of `new` (public summary) against `ref` with absolute and relative change."""
    rows = []
    for label, path, unit in KEY_NUMBERS:
        a, b = _dig(ref, path), _dig(new, path)
        if a is None or b is None:
            continue
        rows.append({"quantity": label, "unit": unit, "placeholder": a, "m4": b, "delta": b - a,
                     "rel": (b / a - 1.0) if a else float("nan")})
    dom_ref, dom_new = _dig(ref, ("torsional", "dominant", "order")), _dig(new, ("torsional", "dominant", "order"))
    rows.append({"quantity": "dominant torsional order", "unit": "-", "placeholder": dom_ref, "m4": dom_new,
                 "delta": (dom_new or 0) - (dom_ref or 0), "rel": float("nan")})
    return rows


def compare_text(rows: list[dict]) -> str:
    L = ["Key numbers, REAL (M4 / thermo) vs PLACEHOLDER pressure (all MODELLED):",
         f"  {'quantity':40s} {'placeholder':>12s} {'m4':>12s} {'change':>9s}"]
    for r in rows:
        rel = "" if r["rel"] != r["rel"] else f"{r['rel'] * 100:+8.1f}%"
        L.append(f"  {r['quantity'] + ' [' + r['unit'] + ']':40s} {r['placeholder']:12.4g} {r['m4']:12.4g} {rel:>9s}")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pressure", choices=["placeholder", "m4"], default="placeholder",
                    help="cylinder-pressure source: placeholder Wiebe stand-in or REAL M4/thermo traces")
    ap.add_argument("--no-orbits", action="store_true", help="skip the journal-orbit integration")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--no-compare", action="store_true", help="with --pressure m4: skip the placeholder comparison")
    args = ap.parse_args(argv)
    res = compute(with_orbits=not args.no_orbits, pressure=args.pressure)
    print(summary_text(res))
    out = OUT_DIR if args.pressure == "placeholder" else OUT_DIR / args.pressure
    out.mkdir(parents=True, exist_ok=True)
    pub = _public(res)
    (out / "summary.json").write_text(json.dumps(pub, indent=2))
    if args.pressure != "placeholder" and not args.no_compare:
        ref = _public(compute(with_orbits=not args.no_orbits, pressure="placeholder"))
        rows = compare(pub, ref)
        print("\n" + compare_text(rows))
        (out / "compare.json").write_text(json.dumps({"modelled": True, "rows": rows}, indent=2))
    if not args.no_plots:
        files = plots(res, out)
        print("\nplots: " + ", ".join(Path(f).name for f in files) + f"  -> {out}")


if __name__ == "__main__":
    main()
