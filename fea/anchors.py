"""Analytic anchors for the FEA chain (gmsh C3D10 mesh -> CalculiX -> .frd reader), CLAUDE.md rule 2:
each is checked in tests/test_m6_fea.py before the part models are trusted. Dimensions are test
fixtures (stated here), not engine parameters; steel-like E and nu.

1. Cantilever, tip load: Euler-Bernoulli tip deflection P L^3 / (3 E I) and bending stress
   P (L - x) c / I at x = L/4 (the clamped root itself is a stress singularity of the 3D model;
   St Venant: x = L/4 is 5 depths from the clamp). L/h = 20, so shear deflection is ~0.2%.
2. Thick-walled cylinder, internal pressure, plane strain (Lame): sigma_r = A - B / r^2,
   sigma_theta = A + B / r^2, A = p ri^2 / (ro^2 - ri^2), B = p ri^2 ro^2 / (ro^2 - ri^2).
3. Plate with a central hole in tension: Kt (gross) = (2 + (1 - d/W)^3) / (1 - d/W) (Heywood finite-
   width fit, Pilkey & Pilkey, Peterson's SCF 3rd ed., eq. 4.57), = 3.03 at d/W = 0.1; infinite plate 3.
4. 1-D steady conduction in a slab: T(0) fixed, flux q into x = L, all else adiabatic:
   T(x) = T0 + q x / k.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from fea import ccx
from fea.common import cosine_pressure, plane_nodes
from fea.geometry import cantilever, export, lame_quarter, plate_hole_quarter, slab
from fea.mesh import mesh_step

MM = 1.0e-3
E_STEEL, NU_STEEL = 210.0e+9, 0.3


def _mesh(shape, wd: Path, name: str, h: float, curvature_pts: int = 0):
    step = export(shape, wd / f"{name}.step")
    return mesh_step(step, MM, h, h / 10, curvature_pts)


def cantilever_case(wd: Path, h: float = 5.0) -> dict:
    L, b, hh, P = 400.0, 20.0, 20.0, 1000.0
    m = _mesh(cantilever(L, b, hh), wd, "cantilever", h)
    deck = ccx.Deck(m, ccx.Material("steel", E_STEEL, NU_STEEL))
    for n in plane_nodes(m, 0, 0.0):
        deck.boundaries.append((int(n), 1, 3, 0.0))
    tip = m.select_faces(lambda c, n: (np.abs(c[:, 0] - L * MM) < 1e-9) & (n[:, 0] > 0.99))
    ref = deck.add_ref_node([L * MM, 0.0, 0.0])
    deck.couplings.append((ref, tip, "distributing"))
    deck.steps.append(ccx.Step(cloads=[(ref, 3, -P)]))
    r = ccx.run(deck, wd, "cantilever")[0]
    I = b * hh**3 / 12 * MM**4
    tipn = plane_nodes(m, 0, L * MM)
    w_fea = float(r["DISP"][tipn, 2].mean())
    w_eb = -P * (L * MM) ** 3 / (3 * E_STEEL * I)
    x0 = L * MM / 4
    top = np.nonzero((np.abs(m.nodes[:, 0] - x0) < h * MM / 2) & (np.abs(m.nodes[:, 2] - hh * MM / 2) < 1e-9))[0]
    xs = m.nodes[top, 0]
    s_fea = float(np.polyval(np.polyfit(xs, r["STRESS"][top, 0], 1), x0))
    s_eb = P * (L * MM - x0) * (hh * MM / 2) / I
    return {"tip_fea_m": w_fea, "tip_eb_m": w_eb, "tip_rel_err": w_fea / w_eb - 1.0,
            "stress_fea_Pa": s_fea, "stress_eb_Pa": s_eb, "stress_rel_err": s_fea / s_eb - 1.0, "nodes": m.n_nodes}


def lame_case(wd: Path, h: float = 2.0) -> dict:
    ri, ro, t, p = 20.0, 40.0, 10.0, 100.0e+6
    m = _mesh(lame_quarter(ri, ro, t), wd, "lame", h, curvature_pts=24)
    deck = ccx.Deck(m, ccx.Material("steel", E_STEEL, NU_STEEL))
    for n in plane_nodes(m, 0):
        deck.boundaries.append((int(n), 1, 1, 0.0))
    for n in plane_nodes(m, 1):
        deck.boundaries.append((int(n), 2, 2, 0.0))
    for z in (0.0, t * MM):
        for n in plane_nodes(m, 2, z):
            deck.boundaries.append((int(n), 3, 3, 0.0))      # plane strain
    _e, _f, _cn, _six, c, _n, _a = m.boundary_faces()
    sel = np.nonzero(np.abs(np.hypot(c[:, 0], c[:, 1]) - ri * MM) < 0.02 * ri * MM)[0]
    e, f = m.boundary_faces()[:2]
    deck.steps.append(ccx.Step(dloads=[(int(e[s]), int(f[s]), p) for s in sel]))
    S = ccx.run(deck, wd, "lame")[0]["STRESS"][: m.n_nodes]
    x, y = m.nodes[:, 0], m.nodes[:, 1]
    r = np.hypot(x, y)
    cth, sth = x / r, y / r
    sr = S[:, 0] * cth**2 + S[:, 1] * sth**2 + 2 * S[:, 3] * sth * cth
    st = S[:, 0] * sth**2 + S[:, 1] * cth**2 - 2 * S[:, 3] * sth * cth
    A = p * ri**2 / (ro**2 - ri**2)
    B = p * ri**2 * ro**2 / (ro**2 - ri**2)
    sr_a = A - B / (r / MM) ** 2
    st_a = A + B / (r / MM) ** 2
    scale = A + B / ri**2
    return {"hoop_max_err_rel": float(np.abs(st - st_a).max() / scale),
            "radial_max_err_rel": float(np.abs(sr - sr_a).max() / scale),
            "hoop_inner_fea_Pa": float(st[np.argmin(r)]), "hoop_inner_lame_Pa": float(scale), "nodes": m.n_nodes}


def plate_hole_case(wd: Path, h: float = 4.0, h_hole: float = 0.3) -> dict:
    W, Hh, rad, t, sig = 100.0, 100.0, 10.0, 2.0, 100.0e+6
    step = export(plate_hole_quarter(W, Hh, rad, t), wd / "plate.step")
    refine = [{"expr": f"Sqrt(x*x+y*y)-{rad}", "size": h_hole, "dist_min": 0.0, "dist_max": 4 * rad}]
    m = mesh_step(step, MM, h, h_hole / 2, 0, refine=refine)
    deck = ccx.Deck(m, ccx.Material("steel", E_STEEL, NU_STEEL))
    for n in plane_nodes(m, 0):
        deck.boundaries.append((int(n), 1, 1, 0.0))
    for n in plane_nodes(m, 1):
        deck.boundaries.append((int(n), 2, 2, 0.0))
    z0 = plane_nodes(m, 2, 0.0)
    deck.boundaries.append((int(z0[0]), 3, 3, 0.0))
    sel = m.select_faces(lambda c, n: (np.abs(c[:, 0] - W * MM) < 1e-9) & (n[:, 0] > 0.99))
    dl, _ = cosine_pressure(m, sel, np.ones(len(sel)), (1, 0, 0), sig * Hh * MM * t * MM)
    deck.steps.append(ccx.Step(dloads=dl))
    S = ccx.run(deck, wd, "plate")[0]["STRESS"][: m.n_nodes]
    kt = float(S[:, 0].max() / sig)
    dw = 2 * rad / (2 * W)
    kt_ref = (2 + (1 - dw) ** 3) / (1 - dw)
    return {"kt_fea": kt, "kt_ref": kt_ref, "rel_err": kt / kt_ref - 1.0, "nodes": m.n_nodes}


def slab_case(wd: Path, h: float = 4.0) -> dict:
    L, a, k, q, T0 = 100.0, 10.0, 50.0, 1.0e+5, 300.0
    m = _mesh(slab(L, a), wd, "slab", h)
    deck = ccx.Deck(m, ccx.Material("steel", E_STEEL, NU_STEEL, k=k), t_ref=T0)
    e, f = m.boundary_faces()[:2]
    hot = m.select_faces(lambda c, n: (np.abs(c[:, 0] - L * MM) < 1e-9) & (n[:, 0] > 0.99))
    deck.steps.append(ccx.Step(kind="heat", fixed_temps=[(int(n), T0) for n in plane_nodes(m, 0, 0.0)],
                               dfluxes=[(int(e[s]), int(f[s]), q) for s in hot]))
    T = ccx.run(deck, wd, "slab")[0]["NDTEMP"][: m.n_nodes]
    T_a = T0 + q * m.nodes[:, 0] / k
    return {"max_err_K": float(np.abs(T - T_a).max()), "T_end_fea_K": float(T.max()),
            "T_end_exact_K": T0 + q * L * MM / k, "nodes": m.n_nodes}


def run_all(wd: Path) -> dict:
    return {"cantilever": cantilever_case(wd), "lame": lame_case(wd), "plate_hole": plate_hole_case(wd),
            "slab": slab_case(wd)}
