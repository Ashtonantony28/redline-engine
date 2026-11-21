"""Main bearing cap: bolt preload + peak main-bearing load from dyn.bearings (worst main, speed and
crank angle at WOT), Goodman fatigue of the cap, analytic bolt fatigue and joint separation.

Model: quarter cap (x >= 0 mid plane, y >= 0 symmetry plane), joint face z = 0 as a frictionless
normal support (valid while the joint stays closed, checked analytically). Preload = pressure on the
bolt-head bearing annulus; bearing load = cosine pressure over the lower half bore, directed straight
down with the magnitude of the peak load vector (its small lateral component, reported, is folded into
the magnitude). Bolt stiffness and joint stiffness factor C from Shigley (eq. 8-17, 8-20 frusta,
Shigley "Mechanical Engineering Design" 10th ed.); bolt Goodman with preload, Shigley eq. 8-38.
"""
from __future__ import annotations

import math

import numpy as np

from fea import ccx
from fea.common import (LevelResult, bore_faces, cosine_pressure, mesh_options, plane_nodes, roi_mask, solve,
                        workdir)
from fea.config import material
from fea.fatigue import goodman_sf, signed_von_mises
from fea.geometry import export, main_cap_dims, main_cap_solid, _halfspace
from fea.loads import main_bearing_peak
from fea.mesh import mesh_step

NAME = "main_cap"
SECTION = "main_cap"
FACE_TOL_FRAC = 0.02
NORMAL_TOL = 0.99
CONE_TAN = math.tan(math.radians(30.0))   # Shigley pressure-cone half angle


def build(p, fspec):
    return main_cap_solid(p, half=True) & _halfspace(x=True)


def frustum_k(E, d, D, t):
    """Stiffness of one 30-deg pressure frustum (Shigley eq. 8-20)."""
    a = 2 * CONE_TAN
    return math.pi * E * d * CONE_TAN / math.log((a * t + D - d) * (D + d) / ((a * t + D + d) * (D - d)))


def bolt_check(p, fspec, F_bearing):
    fd = fspec[SECTION]
    b = fspec["materials"]["bolt"]
    cap = fspec["materials"]["main_cap"]
    d = main_cap_dims(p)
    k = p.k
    As, dn = fd["bolt_stress_area"], fd["bolt_nominal_d"]
    grip = (d["depth"]) / k            # bolt clamps the cap from its bottom face to the joint
    Fi = fd["preload_fraction"] * b["sp"] * As
    kb = As * b["E"] / grip
    E_block = float(p.engine["materials"][p.cad["materials"]["block"]]["E"])
    D = d["head_r"] * 2 / k
    k1 = frustum_k(cap["E"], dn, D, grip / 2)
    k2 = frustum_k(E_block, dn, D, grip / 2)
    km = 1.0 / (1.0 / k1 + 1.0 / k2)
    C = kb / (kb + km)
    P = F_bearing / 2                  # per bolt
    sig_i = Fi / As
    sig_a = C * P / (2 * As)
    sig_m = sig_i + sig_a
    nf = b["se"] * (b["sut"] - sig_i) / (b["sut"] * sig_a + b["se"] * (sig_m - sig_i))
    n_sep = Fi / (P * (1 - C))
    n_proof = b["sp"] * As / (C * P + Fi)
    return {"preload_N": Fi, "kb_N_per_m": kb, "km_N_per_m": km, "C": C, "P_per_bolt_N": P, "sigma_i_Pa": sig_i,
            "sigma_a_Pa": sig_a, "sigma_m_Pa": sig_m, "fatigue_sf": nf, "separation_sf": n_sep,
            "proof_load_factor": n_proof, "grip_m": grip}


def analyse(p, fspec, level: int, loads: dict | None = None, tag: str = "") -> LevelResult:
    fd = fspec[SECTION]
    loads = loads or main_bearing_peak(p, fspec)
    step = export(build(p, fspec), workdir(fspec, NAME) / f"{NAME}{tag}.step")
    h, hmin, cp = mesh_options(fd["mesh"], level)
    mesh = mesh_step(step, 1.0 / p.k, h, hmin, cp)
    d = main_cap_dims(p)
    k = p.k
    cm = fspec["materials"]["main_cap"]
    deck = ccx.Deck(mesh, ccx.Material("main_cap", cm["E"], cm["nu"], cm["rho"]))
    for n_ in plane_nodes(mesh, 0):
        deck.boundaries.append((int(n_), 1, 1, 0.0))
    for n_ in plane_nodes(mesh, 1):
        deck.boundaries.append((int(n_), 2, 2, 0.0))
    joint = plane_nodes(mesh, 2, 0.0)
    for n_ in joint:
        deck.boundaries.append((int(n_), 3, 3, 0.0))
    bolt = bolt_check(p, fspec, loads["F_mag_N"])
    e, f, _cn, _six, c, n, a = mesh.boundary_faces()
    tol = FACE_TOL_FRAC * d["hole_r"] / k
    rb = np.hypot(c[:, 0], c[:, 1] - d["bolt_y"] / k)
    seat = np.nonzero((np.abs(c[:, 2] + d["depth"] / k) < tol) & (n[:, 2] < -NORMAL_TOL) &
                      (rb > d["hole_r"] / k) & (rb < d["head_r"] / k))[0]
    half_bolt = 0.5   # the quarter holds half of one bolt's head annulus
    pre_dl, _ = cosine_pressure(mesh, seat, np.ones(len(seat)), (0, 0, 1), half_bolt * bolt["preload_N"])
    sel, cosw = bore_faces(mesh, (0.0, 0.0), d["bore_r"] / k, (0.0, -1.0), FACE_TOL_FRAC * d["bore_r"] / k)
    quarter = 0.25
    brg_dl, _ = cosine_pressure(mesh, sel, cosw, (0, 0, -1), quarter * loads["F_mag_N"])
    deck.steps.append(ccx.Step(dloads=pre_dl))
    deck.steps.append(ccx.Step(dloads=brg_dl))
    res = solve(deck, fspec, NAME, f"{NAME}{tag}_L{level}")
    S_pre = res[0]["STRESS"][: mesh.n_nodes]
    S_b = res[1]["STRESS"][: mesh.n_nodes]
    roi = np.ones(mesh.n_nodes, dtype=bool)
    S_max = S_pre + S_b
    vm = ccx.von_mises(S_max)
    sa = 0.5 * ccx.von_mises(S_b)
    sm = signed_von_mises(S_pre + 0.5 * S_b)
    sf_n = np.where(roi, goodman_sf(sa, sm, cm["se"], cm["sut"]), np.inf)
    i_f = int(np.argmin(sf_n))
    i_p = int(np.argmax(vm))
    # joint stays closed in the FE model if the joint-face reaction stays compressive everywhere is not
    # checked node-wise (linear support); the analytic separation SF covers it
    data = {"loads": loads, "bolt": bolt,
            "peak_vm_Pa": float(vm[i_p]), "peak_at_m": mesh.nodes[i_p].tolist(),
            "peak_vm_preload_only_Pa": float(ccx.von_mises(S_pre).max()),
            "peak_vm_bearing_only_Pa": float(ccx.von_mises(S_b).max()),
            "fatigue_critical": {"sf": float(sf_n[i_f]), "at_m": mesh.nodes[i_f].tolist(), "sigma_a_Pa": float(sa[i_f]),
                                 "sigma_m_Pa": float(sm[i_f]), "se_Pa": cm["se"], "sut_Pa": cm["sut"]},
            "yield_sf": cm["sy"] / float(vm[i_p])}
    sf = {"fatigue": float(sf_n[i_f]), "bolt_fatigue": bolt["fatigue_sf"]}
    return LevelResult(level, mesh.n_nodes, mesh.n_elems, float(vm[i_p]), sf, data, mesh,
                       {"vm": vm, "goodman_sf": sf_n})


def plots(r: LevelResult, outdir, p):
    from fea.common import plot_field
    plot_field(r.mesh, r.fields["vm"] / 1.0e6, outdir / "main_cap_vm.png",
               "Main cap von Mises, preload + peak bearing load (quarter, mirrored)", "von Mises [MPa]",
               view=(-20, -60), mirror=(0, 1))
    return ["main_cap_vm.png"]


_ = roi_mask
