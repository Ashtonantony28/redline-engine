"""Conrod: compression at peak firing (worst WOT speed) minus inertia, tension at redline overrun TDC;
Euler/Johnson buckling of the shank in both planes; FEA stress; Goodman fatigue.

Model: quarter of the reconstructed rod (fea.geometry.conrod_solid), symmetry planes x = 0 (rod mid
plane, ux = 0) and y = 0 (rod plane of symmetry, uy = 0). The bearing loads are cosine-distributed
pressures over the loaded half of each bore (pin on the small end, crankpin on the big end), equal and
opposite, so the model is self-equilibrated; a single node at the big-end cap bottom fixes uz (its
reaction is zero up to round-off). Bearing shells, bolts/cap joint and the bore clearance are not
modelled (documented limits). Units SI, frame of cad/conrod.py: big end at z = 0, small end at z = L.
"""
from __future__ import annotations

import math

import numpy as np

from fea import ccx
from fea.common import (LevelResult, bore_faces, cosine_pressure, goodman_nodes, mesh_options, plane_nodes,
                        roi_mask, solve, workdir)
from fea.config import material
from fea.fatigue import johnson_euler
from fea.geometry import conrod_dims, conrod_solid, export
from fea.loads import conrod_loads
from fea.mesh import mesh_step

NAME = "conrod"
SECTION = "conrod"
FACE_TOL_FRAC = 0.02     # bore-face radius match tolerance as a share of the bore radius


def section_properties(p, z_cad: float) -> dict:
    """I-section of the shank at height z (CAD units) from spec/cad.yaml, SI. Bending 'in plane'
    = deflection in y (swing plane), about the x axis; 'out of plane' = deflection in x."""
    d = conrod_dims(p)
    k = p.k
    h = d["hb"] + (d["hs"] - d["hb"]) * z_cad / d["L"]          # half depth in y
    b, ft, tw = 2 * d["fw"], d["ft"], 2 * d["tw"]
    hi = h - ft
    A = b * 2 * h - (b - tw) * 2 * hi
    I_in = (b * (2 * h) ** 3 - (b - tw) * (2 * hi) ** 3) / 12.0
    I_out = (2 * ft * b**3 + 2 * hi * tw**3) / 12.0
    return {"A": A / k**2, "I_in": I_in / k**4, "I_out": I_out / k**4, "half_depth": h / k}


def buckling(p, fd, F_comp: float) -> dict:
    """Shank buckling with the smallest section (pocket end at the small eye) over the centre distance."""
    d = conrod_dims(p)
    z_min = d["L"] - d["se_R"] - p.L(fd["defeature"]["pocket_end_gap"])
    sp = section_properties(p, z_min)
    m = material(p.engine, p.cad["materials"]["conrod"])
    L = p.geom.conrod_length
    out = {"section": sp, "L": L}
    for plane, I, K in (("in_plane", sp["I_in"], fd["buckling"]["in_plane_K"]),
                        ("out_of_plane", sp["I_out"], fd["buckling"]["out_of_plane_K"])):
        P, regime, sr, cc = johnson_euler(sp["A"], I, K * L, m["E"], m["sy"])
        out[plane] = {"P_cr_N": P, "regime": regime, "slenderness": sr, "C_c": cc, "K": K, "sf": P / F_comp}
    out["sf"] = min(out["in_plane"]["sf"], out["out_of_plane"]["sf"])
    return out


def build(p, fspec):
    return conrod_solid(p, fspec[SECTION], quarter=True)


def analyse(p, fspec, level: int, loads: dict | None = None, tag: str = "") -> LevelResult:
    fd = fspec[SECTION]
    loads = loads or conrod_loads(p, fspec)
    step = export(build(p, fspec), workdir(fspec, NAME) / f"{NAME}{tag}.step")
    h, hmin, cp = mesh_options(fd["mesh"], level)
    mesh = mesh_step(step, 1.0 / p.k, h, hmin, cp, threads=int(fspec["solver"]["threads"]))
    d = conrod_dims(p)
    k = p.k
    L, be_r, se_r = d["L"] / k, d["be_r"] / k, d["se_r"] / k
    m = material(p.engine, p.cad["materials"]["conrod"])
    nu = fspec["materials"]["poisson"][p.cad["materials"]["conrod"]]
    deck = ccx.Deck(mesh, ccx.Material("steel_rod", m["E"], nu, m["rho"]))
    for n in plane_nodes(mesh, 0):
        deck.boundaries.append((int(n), 1, 1, 0.0))
    for n in plane_nodes(mesh, 1):
        deck.boundaries.append((int(n), 2, 2, 0.0))
    anchor_pt = np.array([0.0, 0.0, -d["be_R"] / k])
    anchor = int(np.argmin(np.linalg.norm(mesh.nodes - anchor_pt, axis=1)))
    deck.boundaries.append((anchor, 3, 3, 0.0))
    quarter = 0.25
    unit = 1.0
    steps = {}
    for case, se_dir, be_dir in (("compression", (0.0, -1.0), (0.0, 1.0)), ("tension", (0.0, 1.0), (0.0, -1.0))):
        s_sel, s_cos = bore_faces(mesh, (0.0, L), se_r, se_dir, FACE_TOL_FRAC * se_r)
        b_sel, b_cos = bore_faces(mesh, (0.0, 0.0), be_r, be_dir, FACE_TOL_FRAC * be_r)
        dl_s, _ = cosine_pressure(mesh, s_sel, s_cos, (0, 0, se_dir[1]), quarter * unit)
        dl_b, _ = cosine_pressure(mesh, b_sel, b_cos, (0, 0, be_dir[1]), quarter * unit)
        deck.steps.append(ccx.Step(dloads=dl_s + dl_b))
        steps[case] = len(deck.steps) - 1
    res = solve(deck, fspec, NAME, f"{NAME}{tag}_L{level}")
    S_c = res[steps["compression"]]["STRESS"][: mesh.n_nodes] * loads["compression_N"]
    S_t = res[steps["tension"]]["STRESS"][: mesh.n_nodes] * loads["tension_N"]
    roi = roi_mask(mesh, [anchor_pt], fd["roi_exclusion"])
    vm_c, vm_t = ccx.von_mises(S_c), ccx.von_mises(S_t)
    sf_n, sa, sm = goodman_nodes(0.5 * (S_t - S_c), 0.5 * (S_t + S_c), m["se"], m["sut"])
    sf_n = np.where(roi, sf_n, np.inf)
    i_f = int(np.argmin(sf_n))
    i_c = int(np.argmax(np.where(roi, vm_c, 0.0)))
    i_t = int(np.argmax(np.where(roi, vm_t, 0.0)))
    bk = buckling(p, fd, loads["compression_N"])
    peak = float(max(vm_c[i_c], vm_t[i_t]))      # monitored for convergence: governing case
    data = {
        "loads": {k2: v for k2, v in loads.items() if k2 != "per_rpm"},
        "buckling": bk,
        "peak_vm_compression_Pa": float(vm_c[i_c]), "peak_vm_compression_at_m": mesh.nodes[i_c].tolist(),
        "peak_vm_tension_Pa": float(vm_t[i_t]), "peak_vm_tension_at_m": mesh.nodes[i_t].tolist(),
        "fatigue_critical": {"at_m": mesh.nodes[i_f].tolist(), "sigma_a_Pa": float(sa[i_f]),
                             "sigma_m_Pa": float(sm[i_f]), "se_Pa": m["se"], "sut_Pa": m["sut"]},
        "yield_sf": m["sy"] / peak,
        "nominal_shank_stress_Pa": loads["compression_N"] / bk["section"]["A"],
    }
    sf = {"buckling": bk["sf"], "fatigue": float(sf_n[i_f])}
    return LevelResult(level, mesh.n_nodes, mesh.n_elems, peak, sf, data, mesh,
                       {"vm_compression": vm_c, "vm_tension": vm_t, "goodman_sf": sf_n})


def plots(r: LevelResult, outdir, p):
    from fea.common import plot_field
    plot_field(r.mesh, r.fields["vm_compression"] / 1.0e6, outdir / "conrod_vm_compression.png",
               "Conrod von Mises, peak firing compression (quarter model, mirrored)", "von Mises [MPa]",
               view=(15, 30), mirror=(0, 1))
    plot_field(r.mesh, r.fields["vm_tension"] / 1.0e6, outdir / "conrod_vm_tension.png",
               "Conrod von Mises, redline overrun tension", "von Mises [MPa]", view=(15, 30), mirror=(0, 1))
    return ["conrod_vm_compression.png", "conrod_vm_tension.png"]


def mass_check(p) -> dict:
    from fea.loads import cad_masses
    m = cad_masses(p)
    budget = p.engine["masses"]
    return {"conrod_kg": m["conrod"], "conrod_budget_kg": budget["conrod"],
            "conrod_rel_err": m["conrod"] / budget["conrod"] - 1.0}


_ = math
