"""Piston: steady thermal field with the Woschni crown heat flux, thermal stress, gas load at peak
firing pressure; crown and pin-boss stresses against the piston alloy's cited strength at the local
temperature (spec/engine.yaml tables); Goodman fatigue (gas load cycles 0 -> peak on top of the steady
thermal stress) and a static check S_ut(T) / von Mises.

Model: quarter of the defeatured piston (fea.geometry.piston_solid), symmetry planes x = 0 (pin-axis
normal plane... thrust plane, ux = 0) and y = 0 (pin plane, uy = 0); adiabatic on both for the thermal
problem. Gas load (p_max - p_crankcase) on the crown; the pin reaction is a cosine pressure over the
upper half of the pin bore along the boss, equal to the gas force (self-equilibrated); one node at the
boss bottom fixes uz. Reciprocating inertia at TDC (which would relieve the pin load) is ignored
(conservative); ring-land gas pressure and side thrust are not modelled. Thermal boundary conditions:
crown = cycle-mean Woschni flux (thermo/) at the spec wall temperature times piston_area_factor; ring
belt film to the liner temperature; underside splash film plus the piston cooling jets (Stevens & Webb
impingement h, thermo/oiljet.py) and pin bore film to oil (spec/fea.yaml, spec/thermo.yaml). The area-mean
crown temperature is an output and must match thermo.yaml T_piston within consistency_tol_K (ADR 0006).
E held constant with temperature; S_ut(T), S_e(T) from the material tables (fea.config.strength_at).
"""
from __future__ import annotations

import numpy as np

from fea import ccx
from fea.common import (LevelResult, bore_faces, cosine_pressure, mesh_options, plane_nodes, roi_mask, solve,
                        workdir)
from fea.config import material, strength_at
from fea.fatigue import goodman_sf, signed_von_mises
from fea.geometry import export, piston_dims, piston_solid
from fea.loads import piston_loads
from fea.mesh import mesh_step

NAME = "piston"
SECTION = "piston"
FACE_TOL_FRAC = 0.02
NORMAL_TOL = 0.99         # |n . axis| above this counts as a face normal to that axis


def build(p, fspec):
    return piston_solid(p, fspec[SECTION], quarter=True)


def classify_faces(mesh, d, k):
    """Indices of boundary faces: crown top, ring belt outer, pin bore, symmetry planes, underside."""
    _e, _f, _cn, _six, c, n, _a = mesh.boundary_faces()
    top, rg, pin_r = d["top"] / k, d["R"] / k, d["pin_r"] / k
    tol = FACE_TOL_FRAC * pin_r
    rxy = np.hypot(c[:, 0], c[:, 1])
    crown = (np.abs(c[:, 2] - top) < tol) & (n[:, 2] > NORMAL_TOL)
    belt = (np.abs(rxy - rg) < tol) & (np.abs(n[:, 0] * c[:, 0] + n[:, 1] * c[:, 1]) / np.maximum(rxy, 1e-30) > NORMAL_TOL)
    land = (rxy > rg + tol)                                            # (none with filled grooves)
    pin = np.abs(np.hypot(c[:, 1], c[:, 2]) - pin_r) < tol
    symx = (np.abs(c[:, 0]) < tol) & (np.abs(n[:, 0]) > NORMAL_TOL)
    symy = (np.abs(c[:, 1]) < tol) & (np.abs(n[:, 1]) > NORMAL_TOL)
    under = ~(crown | belt | land | pin | symx | symy)
    return {k2: np.nonzero(v)[0] for k2, v in
            {"crown": crown, "belt": belt, "land": land, "pin": pin, "under": under}.items()}


def underside_films(mesh, under, th, p, rpm):
    """Film coefficient on each underside face: splash everywhere, max(splash, jet) on faces that face the
    piston cooling jets (thermo/oiljet.py, Stevens & Webb 1991). Jets sit on the x = 0 (thrust-axis) plane at
    y = +-target_offset; the quarter model (x, y >= 0) sees the +y jet, its mirror images are the same jet
    (x mirror) and the twin orifice (y mirror). Returns (h per face, jet summary)."""
    from thermo.oiljet import jet_spec, jet_state, jets_flowing, local_h, stagnation_nusselt
    h = np.full(len(under), float(th["underside_h_splash"]))
    js = jet_spec()
    info = {"enabled": bool(js["enabled"])}
    info["flowing_at_rpm"] = bool(jets_flowing(rpm))
    if not info["flowing_at_rpm"]:
        return h, info
    jt = th["jets"]
    st = jet_state()
    y0 = float(js["target_offset"])
    # geometric check: the jet must clear the conrod small end (outer radius = pin radius + small-end wall)
    r_se = p.section("wrist_pin")["outer_diameter"] / 2 + p.section("conrod")["small_end_wall"]
    clear = y0 - st.d / 2 - r_se
    if clear < float(js["min_small_end_clearance"]):
        raise ValueError(f"cooling jet clears the small end by {clear * 1e3:.2f} mm only")
    _e, _f, _cn, _six, c, n, a = mesh.boundary_faces()
    cu, nu_ = c[under], n[under]                          # mesh coordinates are in m
    r = np.hypot(cu[:, 0], np.abs(cu[:, 1]) - y0)          # distance from the nearest jet axis
    z = float(js["nozzle_to_crown_mean"])
    hj = local_h(r, st, z, jt["correlation"])
    facing = nu_[:, 2] < float(jt["facing_nz"])
    h = np.where(facing, np.maximum(h, hj), h)
    h0 = stagnation_nusselt(st, z, jt["correlation"]) * st.k / st.d
    from thermo.oiljet import fmep_jets, pump_power
    info["pump_cost"] = {"pump_power_W": pump_power(rpm), "fmep_bar_at_rpm": {
        str(int(n_)): fmep_jets(float(n_)) / 1.0e5 for n_ in (js["on_above_rpm"], 4500.0, rpm)}}
    info.update({"orifice_d_m": st.d, "u_ms": st.u, "Re": st.re, "Pr": st.pr, "q_orifice_m3s": st.q_orifice,
                 "q_engine_Lmin": st.q_total * 6.0e4, "target_offset_m": y0, "small_end_clearance_m": clear,
                 "h_stagnation_Wm2K": float(h0), "faces_above_splash": int((h > th["underside_h_splash"]).sum()),
                 "area_above_splash_m2_quarter": float(a[under][h > th["underside_h_splash"]].sum())})
    return h, info


def analyse(p, fspec, level: int, loads: dict | None = None, tag: str = "") -> LevelResult:
    fd = fspec[SECTION]
    th = fd["thermal"]
    loads = loads or piston_loads(fspec)
    step = export(build(p, fspec), workdir(fspec, NAME) / f"{NAME}{tag}.step")
    h, hmin, cp = mesh_options(fd["mesh"], level)
    mesh = mesh_step(step, 1.0 / p.k, h, hmin, cp)
    d = piston_dims(p)
    k = p.k
    key = p.cad["materials"]["piston"]
    m = material(p.engine, key)
    tm = fspec["materials"]["piston_thermal"][key]
    mat = ccx.Material(key, m["E"], fspec["materials"]["poisson"][key], m["rho"], tm["alpha"], tm["k"])
    faces = classify_faces(mesh, d, k)
    e, f, *_ = mesh.boundary_faces()
    # ---- thermal ----
    q_crown = loads["q_mean_Wm2"] * loads["piston_area_factor"]
    crown_faces = faces["crown"]
    bf = mesh.boundary_faces()
    crown_six, crown_area = bf[3][crown_faces], bf[6][crown_faces]
    h_under, jet = underside_films(mesh, faces["under"], th, p, loads["q_rpm"])

    heat = ccx.Deck(mesh, mat, t_ref=tm["t_ref"])
    films = [(int(e[i]), int(f[i]), loads["T_liner_K"], th["ring_belt_h"]) for i in faces["belt"]]
    films += [(int(e[i]), int(f[i]), th["T_oil"], float(hu)) for i, hu in zip(faces["under"], h_under)]
    films += [(int(e[i]), int(f[i]), th["T_oil"], th["pin_bore_h"]) for i in faces["pin"]]
    dfl = [(int(e[i]), int(f[i]), q_crown) for i in crown_faces]
    heat.steps.append(ccx.Step(kind="heat", films=films, dfluxes=dfl))
    T = solve(heat, fspec, NAME, f"{NAME}{tag}_heat_L{level}")[0]["NDTEMP"][: mesh.n_nodes]
    t_crown = float((T[crown_six].mean(axis=1) * crown_area).sum() / crown_area.sum())
    # heat balance: what enters the crown leaves through the films (steady state)
    a_all, c_all = bf[6], bf[4]
    six = bf[3]
    def film_heat(idx, sink, hh):
        tf = T[six[idx]].mean(axis=1)
        return float((np.asarray(hh) * a_all[idx] * (tf - sink)).sum())
    q_in = float(q_crown * crown_area.sum())
    q_belt = film_heat(faces["belt"], loads["T_liner_K"], th["ring_belt_h"])
    q_under = film_heat(faces["under"], th["T_oil"], h_under)
    q_pin = film_heat(faces["pin"], th["T_oil"], th["pin_bore_h"])
    jet_mask = h_under > th["underside_h_splash"]
    q_jet_excess = film_heat(faces["under"][jet_mask], th["T_oil"], h_under[jet_mask] - th["underside_h_splash"])
    tol = th["consistency_tol_K"]
    thermal = {"crown_mean_T_K": t_crown, "T_piston_spec_K": loads["T_piston_spec_K"],
               "consistency_delta_K": t_crown - loads["T_piston_spec_K"], "consistency_tol_K": tol,
               "consistent": abs(t_crown - loads["T_piston_spec_K"]) <= tol,
               "h_splash_Wm2K": th["underside_h_splash"], "jets": jet,
               "heat_W_quarter": {"crown_in": q_in, "belt_out": q_belt, "underside_out": q_under,
                                  "pin_bore_out": q_pin, "jet_excess_over_splash": q_jet_excess,
                                  "balance_rel": (q_in - q_belt - q_under - q_pin) / q_in}}
    # ---- mechanical ----
    deck = ccx.Deck(mesh, mat, t_ref=tm["t_ref"])
    for n_ in plane_nodes(mesh, 0):
        deck.boundaries.append((int(n_), 1, 1, 0.0))
    for n_ in plane_nodes(mesh, 1):
        deck.boundaries.append((int(n_), 2, 2, 0.0))
    anchor_pt = np.array([(d["boss_gap"] + d["flat_x"]) / 2 / k, 0.0, -d["boss_r"] / k])
    anchor = int(np.argmin(np.linalg.norm(mesh.nodes - anchor_pt, axis=1)))
    deck.boundaries.append((anchor, 3, 3, 0.0))
    dp = loads["p_max_Pa"] - loads["crankcase_pressure_Pa"]
    crown_dl = [(int(e[i]), int(f[i]), dp) for i in faces["crown"]]
    a, nrm = mesh.boundary_faces()[6], mesh.boundary_faces()[5]
    F_gas = float((dp * a[faces["crown"]] * nrm[faces["crown"], 2]).sum())    # downward resultant magnitude (> 0)
    sel, cosw = bore_faces(mesh, (0.0, 0.0), d["pin_r"] / k, (0.0, 1.0), FACE_TOL_FRAC * d["pin_r"] / k)
    pin_dl, _ = cosine_pressure(mesh, sel, cosw, (0, 0, 1), F_gas)
    deck.steps.append(ccx.Step(dloads=crown_dl + pin_dl, temperatures=np.full(mesh.n_nodes, tm["t_ref"])))
    deck.steps.append(ccx.Step(temperatures=T))
    res = solve(deck, fspec, NAME, f"{NAME}{tag}_L{level}")
    S_m = res[0]["STRESS"][: mesh.n_nodes]
    S_t = res[1]["STRESS"][: mesh.n_nodes]
    roi = roi_mask(mesh, [anchor_pt], fd["roi_exclusion"])
    S_peak = S_m + S_t
    vm = ccx.von_mises(S_peak)
    sut_T, se_T = strength_at(T, m)
    sa = 0.5 * ccx.von_mises(S_m)
    sm = signed_von_mises(S_t + 0.5 * S_m)
    sf_n = np.where(roi, goodman_sf(sa, sm, se_T, sut_T), np.inf)
    static_n = np.where(roi, sut_T / np.maximum(vm, 1.0), np.inf)
    i_f = int(np.argmin(sf_n))
    i_p = int(np.argmax(np.where(roi, vm, 0.0)))
    # regions for reporting: crown (above the underside) and pin boss (around the pin axis)
    z = mesh.nodes[:, 2]
    crown_r = roi & (z > d["crown_under"] / k - p.L(fd["defeature"]["crown_under_fillet"]) / k)
    boss_r = roi & (np.hypot(mesh.nodes[:, 1], z) < (d["boss_r"] / k) * 1.05) & (mesh.nodes[:, 0] > d["boss_gap"] / k * 0.999)
    def region(mask):
        j = int(np.argmax(np.where(mask, vm, 0.0)))
        jf = int(np.argmin(np.where(mask, sf_n, np.inf)))
        return {"peak_vm_Pa": float(vm[j]), "at_m": mesh.nodes[j].tolist(), "T_K": float(T[j]),
                "sut_derated_Pa": float(sut_T[j]), "static_sf": float(sut_T[j] / vm[j]),
                "fatigue_sf": float(sf_n[jf]), "fatigue_at_m": mesh.nodes[jf].tolist(), "fatigue_T_K": float(T[jf]),
                "se_derated_Pa": float(se_T[jf]), "sigma_a_Pa": float(sa[jf]), "sigma_m_Pa": float(sm[jf])}
    centre = int(np.argmin(np.linalg.norm(mesh.nodes - np.array([0.0, 0.0, d["top"] / k]), axis=1)))
    data = {"loads": {k2: v for k2, v in loads.items() if k2 != "sweep"}, "q_crown_applied_Wm2": q_crown,
            "thermal": thermal,
            "F_gas_full_N": 4 * F_gas,
            "T_max_K": float(T.max()), "T_crown_centre_K": float(T[centre]), "T_min_K": float(T.min()),
            "crown": region(crown_r), "pin_boss": region(boss_r),
            "fatigue_critical": {"sf": float(sf_n[i_f]), "at_m": mesh.nodes[i_f].tolist(), "T_K": float(T[i_f]),
                                 "sigma_a_Pa": float(sa[i_f]), "sigma_m_Pa": float(sm[i_f]),
                                 "se_derated_Pa": float(se_T[i_f]), "sut_derated_Pa": float(sut_T[i_f]),
                                 "S_mech_Pa": S_m[i_f].tolist(), "S_thermal_Pa": S_t[i_f].tolist()},
            "static_sf_min": float(static_n.min()),
            "peak_vm_Pa": float(vm[i_p]), "peak_at_m": mesh.nodes[i_p].tolist(),
            "peak_vm_mech_only_Pa": float(ccx.von_mises(S_m)[roi].max()),
            "peak_vm_thermal_only_Pa": float(ccx.von_mises(S_t)[roi].max())}
    return LevelResult(level, mesh.n_nodes, mesh.n_elems, float(vm[i_p]), {"fatigue": float(sf_n[i_f]), "static": float(static_n.min())}, data, mesh,
                       {"T": T, "vm": vm, "goodman_sf": sf_n, "S_m": S_m, "S_t": S_t})


def plots(r: LevelResult, outdir, p):
    from fea.common import plot_field
    plot_field(r.mesh, r.fields["T"], outdir / "piston_temperature.png",
               "Piston steady temperature (quarter model, mirrored)", "T [K]", view=(25, -50), mirror=(0, 1))
    plot_field(r.mesh, r.fields["vm"] / 1.0e6, outdir / "piston_vm.png",
               "Piston von Mises, peak pressure + thermal", "von Mises [MPa]", view=(-25, -50), mirror=(0, 1))
    return ["piston_temperature.png", "piston_vm.png"]
