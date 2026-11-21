"""Markdown report and JSON results for `python -m fea.run --all` (out/fea/). All numbers MODELLED."""
from __future__ import annotations

import datetime
import json
from pathlib import Path

MPA = 1.0e-6


def _f(v, fmt="{:.3g}"):
    return "n/a" if v is None else fmt.format(v)


def part_json(st) -> dict:
    return {"sf": st.sf, "targets": st.targets, "passed": st.passed, "convergence": st.convergence,
            "iterations": st.iterations, "final_params": st.final_params, "changed": st.changed,
            "mass_final": getattr(st, "mass_final", None), "mass_ok": getattr(st, "mass_ok", None),
            "data": st.finest.data if st.finest else None,
            "levels": [{"level": r.level, "nodes": r.n_nodes, "elements": r.n_elems, "peak_Pa": r.peak, "sf": r.sf}
                       for r in st.levels],
            "figures": getattr(st, "figures", []), "seconds": st.seconds,
            "stop_reason": getattr(st, "stop_reason", ""), "volume_check": getattr(st, "volume_check", None),
            "probes": getattr(st, "probes", [])}


def write_results(studies, anchors, spec_hash, out: Path) -> dict:
    res = {"note": "MODELLED results: there is no physical engine", "spec_hash": spec_hash,
           "generated": datetime.datetime.now().isoformat(timespec="seconds"),
           "anchors": anchors, "parts": {s.name: part_json(s) for s in studies}}
    res["all_passed"] = all(s.passed for s in studies)
    (out / "results.json").write_text(json.dumps(res, indent=2, default=_np))
    return res


def _np(o):
    import numpy as np
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _conv_table(st) -> list[str]:
    c = st.convergence
    lines = ["| level | nodes | elements | peak von Mises [MPa] | SFs |", "|---|---|---|---|---|"]
    for r in st.levels:
        lines.append(f"| {r.level} | {r.n_nodes} | {r.n_elems} | {r.peak * MPA:.1f} | "
                     + ", ".join(f"{k} {v:.2f}" for k, v in r.sf.items()) + " |")
    if "change_finest_two" in c:
        lines.append("")
        lines.append(f"Change of the peak between the two finest meshes: **{c['change_finest_two']:.2%}** "
                     f"(criterion < 5%: {'met' if c['ok'] else 'NOT met'}). Richardson extrapolation "
                     f"(refinement ratio {_f(c['refinement_ratio'], '{:.2f}')}): observed order "
                     f"{_f(c['richardson_order'], '{:.2f}')}, extrapolated peak "
                     f"{_f(c['richardson_extrapolated_Pa'] and c['richardson_extrapolated_Pa'] * MPA, '{:.1f}')} MPa, "
                     f"GCI {_f(c['gci_fine'], '{:.2%}')}" + (" (non-monotone: no extrapolation)" if c["richardson_order"] is None else "") + ".")
    else:
        lines.append(c.get("note", ""))
    return lines


def _iter_table(st) -> list[str]:
    keys = list(st.iterations[0]["params"])
    lines = ["| iteration | " + " | ".join(f"{k} [mm]" for k in keys) + " | SFs (medium mesh) | piston asm / conrod mass vs budget |",
             "|---" * (len(keys) + 3) + "|"]
    for it in st.iterations:
        m = it["mass"]
        lines.append(f"| {it['iter']} | " + " | ".join(f"{it['params'][k] * 1e3:.2f}" for k in keys) + " | "
                     + ", ".join(f"{k} {v:.2f}" for k, v in it["sf"].items())
                     + f" | {m['piston_rel_err']:+.2%} / {m['conrod_rel_err']:+.2%} |")
    return lines


def write_report(studies, anchors, res, out: Path, history: list | None = None):
    L = ["# Redline M6: structural FEA report", "",
         "> **All results are MODELLED.** There is no physical engine; loads come from the M2 thermo and M5 "
         "dynamics models, geometry from the M1 CAD parameters (spec/cad.yaml).", "",
         f"Generated {res['generated']}, spec hash `{res['spec_hash']}`. Toolchain: gmsh (2nd-order tets, C3D10) "
         "+ CalculiX ccx. Verify: `python -m fea.run --all` and `pytest -q tests/test_m6_fea.py`.", "",
         "## Summary", "",
         "| part | safety factors (finest mesh) | targets | convergence (finest two) | result |",
         "|---|---|---|---|---|"]
    for st in studies:
        c = st.convergence.get("change_finest_two")
        L.append(f"| {st.name} | " + ", ".join(f"{k} **{v:.2f}**" for k, v in st.sf.items()) + " | "
                 + ", ".join(f"{k} >= {v}" for k, v in st.targets.items()) + f" | {_f(c, '{:.2%}')} | "
                 + ("PASS" if st.passed else "**FAIL**") + " |")
    L += ["", f"Overall: **{'all parts meet their targets and the convergence criterion' if res['all_passed'] else 'NOT all targets met'}**.", ""]
    L += ["## Analytic anchors (tests/test_m6_fea.py)", "",
          "| case | FEA | reference | error |", "|---|---|---|---|"]
    a = anchors
    if a:
        c = a["cantilever"]
        L.append(f"| cantilever tip deflection | {c['tip_fea_m'] * 1e3:.4f} mm | {c['tip_eb_m'] * 1e3:.4f} mm (Euler-Bernoulli) | {c['tip_rel_err']:+.2%} |")
        L.append(f"| cantilever bending stress at L/4 | {c['stress_fea_Pa'] * MPA:.2f} MPa | {c['stress_eb_Pa'] * MPA:.2f} MPa | {c['stress_rel_err']:+.2%} |")
        c = a["lame"]
        L.append(f"| Lame cylinder hoop / radial (max error over all nodes, / inner hoop) | {c['hoop_inner_fea_Pa'] * MPA:.1f} MPa inner hoop | {c['hoop_inner_lame_Pa'] * MPA:.1f} MPa | {c['hoop_max_err_rel']:.2%} / {c['radial_max_err_rel']:.2%} |")
        c = a["plate_hole"]
        L.append(f"| plate with hole Kt (gross, d/W = 0.1) | {c['kt_fea']:.3f} | {c['kt_ref']:.3f} (Heywood fit; 3 for an infinite plate) | {c['rel_err']:+.2%} |")
        c = a["slab"]
        L.append(f"| 1-D conduction slab, end temperature | {c['T_end_fea_K']:.3f} K | {c['T_end_exact_K']:.3f} K | max nodal error {c['max_err_K']:.2e} K |")
    L.append("")
    for st in studies:
        L += _part_section(st)
    if history:
        L += ["## Design history (spec/cad.yaml changes applied from earlier runs)", "",
              "| date | part | parameter changes | SFs before -> after (medium mesh) |", "|---|---|---|---|"]
        for h in history:
            if "changes" in h:      # loop-applied parameter changes
                ch = ", ".join(f"{k}: {v[0] * 1e3:.2f} -> {v[1] * 1e3:.2f} mm" for k, v in h["changes"].items())
                L.append(f"| {h['date']} | {h['part']} | {ch} | {h.get('sf_before')} -> {h.get('sf_after')} |")
            else:                   # manual architecture iterations (see the entry in fea/design_history.json)
                its = h.get("iterations", [])
                sf0, sf1 = (its[0]["sf"], its[-1]["sf"]) if its else (None, None)
                L.append(f"| {h['date']} | {h['part']} | {h.get('kind', 'manual change')}: {len(its)} iterations, "
                         f"start for the loop {h.get('accepted_start_for_loop')} | {sf0} -> {sf1} |")
        L.append("")
    L += ["## Method notes and limits", "",
          "- Geometry: each part is reconstructed for meshing from the same spec/cad.yaml parameters as the M1 CAD "
          "(fea/geometry.py), adding the blend radii the CAD omits (sharp re-entrant corners make elastic peak "
          "stress mesh dependent) and the defeaturing listed in spec/fea.yaml. Volumes are compared with the CAD "
          "solids in the part sections.",
          "- Fatigue: modified Goodman with the spec material S_e / S_ut (taken as component-corrected values), "
          "equivalent stresses = von Mises of the alternating / mean tensors (mean signed by the first invariant), "
          "compressive mean stress not credited.",
          "- Loads: thermo WOT sweep (spec/fea.yaml loads.rpm_sweep), dyn forces with CAD masses; worst speed picked "
          "per part. Load cases combine the worst firing case and the redline overrun case where stated (envelope, "
          "conservative).",
          "- Not modelled: bearing-shell contact and clearance (cosine pressure instead), bolt/cap joint of the rod, "
          "oil holes, residual stresses and surface treatments (fillet rolling / nitriding would add margin), "
          "thermal fatigue (low-cycle) of the piston, ring-land gas pressure, piston side thrust.", ""]
    (out / "report.md").write_text("\n".join(L))


def _part_section(st) -> list[str]:
    d = st.finest.data if st.finest else {}
    L = [f"## {st.name}", ""]
    if st.name == "conrod":
        ld, bk = d["loads"], d["buckling"]
        cc = ld["compression_case"]
        L += [f"Loads (CAD masses: piston assembly {ld['masses']['piston_assembly']:.4f} kg, rod {ld['masses']['conrod']:.4f} kg, "
              f"small-end share {ld['masses']['conrod_small_end_fraction']:.3f}, reciprocating {ld['masses']['reciprocating']:.4f} kg):",
              f"- compression {ld['compression_N'] / 1e3:.2f} kN at {cc['rpm']:.0f} rpm, {cc['theta_deg']:.0f} deg ATDC "
              f"(gas {cc['F_gas_N'] / 1e3:.2f} kN, inertia {cc['F_inertia_N'] / 1e3:+.2f} kN), worst WOT speed of the sweep;",
              f"- tension {ld['tension_N'] / 1e3:.2f} kN at {ld['redline_rpm']:.0f} rpm overrun, gas-exchange TDC "
              f"(WOT at redline gives {-(ld['wot_redline_min_F_rod_N'] or 0) / 1e3:.2f} kN).",
              f"- Fatigue cycle: tension <-> compression envelope (conservative).", "",
              "Buckling of the shank (smallest section, at the small-end pocket end, over the centre distance):", "",
              "| plane | K | slenderness KL/rho | C_c | regime | P_cr [kN] | SF |", "|---|---|---|---|---|---|---|"]
        for pl in ("in_plane", "out_of_plane"):
            b = bk[pl]
            L.append(f"| {pl.replace('_', ' ')} | {b['K']} | {b['slenderness']:.1f} | {b['C_c']:.1f} | {b['regime']} | {b['P_cr_N'] / 1e3:.1f} | {b['sf']:.2f} |")
        fc = d["fatigue_critical"]
        L += ["", f"FEA: peak von Mises {d['peak_vm_compression_Pa'] * MPA:.1f} MPa in compression (nominal shank "
              f"{d['nominal_shank_stress_Pa'] * MPA:.1f} MPa), {d['peak_vm_tension_Pa'] * MPA:.1f} MPa in tension; static "
              f"yield SF {d['yield_sf']:.2f}. Fatigue-critical point at {_mm(fc['at_m'])} mm: sigma_a "
              f"{fc['sigma_a_Pa'] * MPA:.1f} MPa, sigma_m {fc['sigma_m_Pa'] * MPA:.1f} MPa (S_e {fc['se_Pa'] * MPA:.0f}, "
              f"S_ut {fc['sut_Pa'] * MPA:.0f} MPa)."]
    elif st.name == "crankshaft":
        fc, ws, kt = d["fatigue_critical"], d["worst_state"], d["kt"]
        L += [f"Load histories: dyn.bearings crankpin forces and dyn.forces torques over the 720-deg cycle for every "
              f"WOT speed and all four throws (front torque = sum of the cylinders ahead), plus the throw's own "
              f"centrifugal load.",
              f"- Worst combined state: throw {ws['throw']} ({ws['side']} half) at {ws['rpm']:.0f} rpm, "
              f"{ws['theta_deg']:.0f} deg: F_radial {ws['Fr_N'] / 1e3:.2f} kN, F_tangential {ws['Ft_N'] / 1e3:.2f} kN, "
              f"T_in {ws['T_in_Nm']:.0f} N m; peak von Mises {ws['vm'] * MPA:.1f} MPa at {_mm(ws['at_m'])} mm.",
              f"- Fatigue-critical: throw {fc['throw']} ({fc['side']}) at {fc['rpm']:.0f} rpm, sigma_a "
              f"{fc['sigma_a_Pa'] * MPA:.1f} MPa, sigma_m {fc['sigma_m_Pa'] * MPA:.1f} MPa at {_mm(fc['at_m'])} mm "
              f"(S_e {d['se_Pa'] * MPA:.0f}, S_ut {d['sut_Pa'] * MPA:.0f} MPa).",
              f"- Fillet stress concentration (unit radial pin load, pin fillet): FEA Kt = {kt['fea']:.2f} on the "
              f"nominal pin bending stress at the pin/web junction; Peterson stepped round bar (D = {kt['D_m'] * 1e3:.0f} mm "
              f"web envelope, d = {kt['d_m'] * 1e3:.0f} mm, r = {kt['r_m'] * 1e3:.1f} mm) gives {kt['peterson']:.2f} "
              f"(ratio {kt['ratio']:.2f}). Same order; the web is not a round shoulder and bends with the pin, so the "
              f"crank fillet sees a higher Kt than the stepped-bar idealisation (consistent with crank-specific SCF "
              f"correlations, which typically give 2.5 to 4 in bending)."]
    elif st.name == "piston":
        ld, fc = d["loads"], d["fatigue_critical"]
        th = d["thermal"]
        jt, hq = th["jets"], th["heat_W_quarter"]
        L += [f"Loads: peak pressure {ld['p_max_Pa'] / 1e5:.1f} bar at {ld['p_max_rpm']:.0f} rpm (largest of the WOT sweep), "
              f"gas force {d['F_gas_full_N'] / 1e3:.2f} kN; crown heat flux {ld['q_mean_Wm2'] / 1e3:.0f} kW/m^2 cycle-mean "
              f"Woschni at {ld['q_rpm']:.0f} rpm and T_wall = {ld['T_piston_spec_K']:.0f} K (h_mean {ld['h_mean_Wm2K']:.0f} "
              f"W/m^2K), x{ld['piston_area_factor']} area factor -> {d['q_crown_applied_Wm2'] / 1e3:.0f} kW/m^2 on the flat crown.",
              f"- Temperatures: max {d['T_max_K']:.0f} K, crown centre {d['T_crown_centre_K']:.0f} K, min {d['T_min_K']:.0f} K. "
              f"Underside splash film {th['h_splash_Wm2K']:.0f} W/m^2K (frozen M6 baseline) plus piston cooling jets: "
              + (f"jets Re {jt['Re']:.0f}, Pr {jt['Pr']:.0f} ({jt['orifice_d_m'] * 1e3:.1f} mm, {jt['u_ms']:.1f} m/s, "
                 f"{jt['q_engine_Lmin']:.1f} L/min engine total), Stevens & Webb (1991) stagnation h "
                 f"{jt['h_stagnation_Wm2K']:.0f} W/m^2K, aimed {jt['target_offset_m'] * 1e3:.0f} mm off axis on the thrust "
                 f"axis ({jt['small_end_clearance_m'] * 1e3:.1f} mm clear of the small end); jets remove "
                 f"{4 * hq['jet_excess_over_splash']:.0f} W above the splash film of {4 * hq['crown_in']:.0f} W into the crown. "
                 f"Cost: oil-pump power {jt['pump_cost']['pump_power_W']:.0f} W above the jet check-valve speed, added to "
                 f"FMEP (gasdyn/friction.py fmep_total): "
                 + ", ".join(f"{v:.3f} bar at {k} rpm" for k, v in jt["pump_cost"]["fmep_bar_at_rpm"].items()) + ". "
                 if jt.get("flowing_at_rpm") else "jets off. ")
              + f"Area-mean crown temperature {th['crown_mean_T_K']:.1f} K vs thermo.yaml T_piston "
              f"{th['T_piston_spec_K']:.0f} K (the wall temperature behind the Woschni flux): "
              f"{th['consistency_delta_K']:+.1f} K, tolerance {th['consistency_tol_K']:.0f} K -> "
              + ("consistent." if th["consistent"] else "NOT consistent: update thermo.yaml T_piston and re-run."),
              "", "| region | peak von Mises [MPa] | T [K] | S_ut(T) [MPa] | static SF | fatigue SF | at [mm] |",
              "|---|---|---|---|---|---|---|"]
        for reg in ("crown", "pin_boss"):
            r = d[reg]
            L.append(f"| {reg.replace('_', ' ')} | {r['peak_vm_Pa'] * MPA:.1f} | {r['T_K']:.0f} | {r['sut_derated_Pa'] * MPA:.0f} | "
                     f"{r['static_sf']:.2f} | {r['fatigue_sf']:.2f} | {_mm(r['fatigue_at_m'])} |")
        L += ["", f"Fatigue-critical point: {_mm(fc['at_m'])} mm at {fc['T_K']:.0f} K, sigma_a {fc['sigma_a_Pa'] * MPA:.1f} MPa, "
              f"sigma_m {fc['sigma_m_Pa'] * MPA:.1f} MPa vs derated S_e {fc['se_derated_Pa'] * MPA:.1f} / S_ut "
              f"{fc['sut_derated_Pa'] * MPA:.1f} MPa. Peak von Mises mechanical only {d['peak_vm_mech_only_Pa'] * MPA:.1f} MPa, "
              f"thermal only {d['peak_vm_thermal_only_Pa'] * MPA:.1f} MPa.",
              "- Strength at temperature: S_ut(T) and S_e(T) of the piston material interpolated from its cited "
              "spec/engine.yaml tables at each node's temperature (forged 2618-T61; sources and the ASSUMED parts of "
              "the tables are listed there); E held at its room-temperature value. Static SF = S_ut(T) / von Mises.",
              "- Not modelled: side thrust on the skirt, reciprocating inertia at the firing TDC (would relieve the pin "
              "load), inertia-driven pin-load reversal at gas-exchange TDC, pin bending / ovalisation (the pin reaction is "
              "a rigid cosine pressure along the boss), ring-land gas pressure, crown valve pockets."]
    else:
        ld, b, fc = d["loads"], d["bolt"], d["fatigue_critical"]
        L += [f"Load: peak main-bearing load into the cap {ld['F_mag_N'] / 1e3:.2f} kN (main {ld['main'] + 1}, {ld['rpm']:.0f} rpm, "
              f"{ld['theta_deg']:.0f} deg; lateral component {ld['Fy_N'] / 1e3:.2f} kN folded into the magnitude).",
              f"- Bolts (2 x M10x1.25, class 10.9): preload {b['preload_N'] / 1e3:.1f} kN each, joint factor C = {b['C']:.3f}, "
              f"bolt fatigue SF {b['fatigue_sf']:.2f}, joint separation SF {b['separation_sf']:.2f}, "
              f"proof-load factor {b['proof_load_factor']:.2f}.",
              f"- Cap: peak von Mises {d['peak_vm_Pa'] * MPA:.1f} MPa at {_mm(d['peak_at_m'])} mm (preload only "
              f"{d['peak_vm_preload_only_Pa'] * MPA:.1f}, bearing only {d['peak_vm_bearing_only_Pa'] * MPA:.1f}); yield SF "
              f"{d['yield_sf']:.2f}; fatigue-critical sigma_a {fc['sigma_a_Pa'] * MPA:.1f} MPa, sigma_m {fc['sigma_m_Pa'] * MPA:.1f} MPa."]
    L += ["", "Mesh convergence:", ""] + _conv_table(st)
    L += ["", "Design iterations (spec/cad.yaml parameters):", ""] + _iter_table(st)
    if getattr(st, "stop_reason", ""):
        L += ["", f"**Loop stopped: {st.stop_reason}.**"]
    if st.changed:
        L += ["", "**Parameters changed by the loop:** " + ", ".join(f"{k} = {v}" for k, v in st.changed.items())
              + (" (apply with `python -m fea.run --all --apply`)." if st.passed else
                 " (not applied to spec/cad.yaml: the target is still missed).")]
    probes = getattr(st, "probes", [])
    if probes:
        L += ["", "Budget-free probes (medium mesh, M1 mass budget ignored; spec/fea.yaml probes):", "",
              "| parameters [mm] | SFs | peak von Mises [MPa] | mass [kg] |", "|---|---|---|---|"]
        for pr in probes:
            L.append("| " + ", ".join(f"{k} {v * 1e3:.1f}" for k, v in pr["params"].items()) + " | "
                     + ", ".join(f"{k} {v:.2f}" for k, v in pr["sf"].items()) + f" | {pr['peak_Pa'] * MPA:.0f} | {pr['mass_kg']:.4f} |")
    for fig in getattr(st, "figures", []):
        L += ["", f"![{fig}]({fig})"]
    vc = getattr(st, "volume_check", None)
    if vc:
        L += ["", f"Reconstruction vs CAD volume: {vc}"]
    L.append("")
    return L


def _mm(xyz):
    return "(" + ", ".join(f"{v * 1e3:.1f}" for v in xyz) + ")"
