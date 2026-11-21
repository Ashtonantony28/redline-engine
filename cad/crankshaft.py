"""Flat-plane crankshaft: 5 mains, 4 crankpins, 8 counterweighted webs, snout and flange.

Frame: crank axis along x (front = -x, ADR 0001 body axes), z up the bore axis, built at
crank angle 0 (cylinder 1 at firing TDC). Crankpin of cylinder k sits at mechanical angle
alpha_k = (0 - firing phase_k) mod 360 from top, i.e. at (y, z) = r (sin alpha, cos alpha):
0/180/180/0 deg for firing order 1-3-4-2. The assembly rotates the shaft about x so a pin
follows (sin(theta - phase), cos(theta - phase)). See spec/cad.yaml `crankshaft`.
"""
from __future__ import annotations

import math

from cad.params import Params
from cad.primitives import cyl_x, prism_x

FAR = 3.0  # wedge extent as a multiple of the counterweight radius (any value > 1/cos works)



def pin_yz(p: Params, idx: int, crank_deg: float = 0.0):
    a = math.radians(p.pin_angle_deg(idx, crank_deg))
    return p.crank_r * math.sin(a), p.crank_r * math.cos(a)


def web(p: Params, idx: int, x0: float, x1: float, with_counterweight: bool = True):
    """One crank web between x0 and x1 for the throw of cylinder idx."""
    c = p.section("crankshaft")
    a = p.L(c["main_diameter"]) / 2 + p.L(c["web_main_margin"])
    b = p.L(c["pin_diameter"]) / 2 + p.L(c["web_pin_margin"])
    py, pz = pin_yz(p, idx)
    d = math.hypot(py, pz)
    u = (py / d, pz / d)
    n = (-u[1], u[0])
    cb = (a - b) / d
    sb = math.sqrt(1.0 - cb**2)
    hull = []
    for s in (1, -1):
        dy, dz = cb * u[0] + s * sb * n[0], cb * u[1] + s * sb * n[1]
        hull.append((a * dy, a * dz))
        hull.append((py + b * dy, pz + b * dz))
    hull = [hull[0], hull[1], hull[3], hull[2]]
    body = cyl_x(a, x0, x1) + cyl_x(b, x0, x1, py, pz) + prism_x(hull, x0, x1)
    if with_counterweight:
        body += counterweight(p, u, x0, x1)
    return body


def counterweight(p: Params, u, x0: float, x1: float):
    """Disc sector of radius counterweight_radius centred opposite the pin direction u."""
    c = p.section("crankshaft")
    R = p.L(c["counterweight_radius"])
    half = math.radians(c["counterweight_angle_deg"]) / 2
    far = FAR * R
    back = math.atan2(-u[0], -u[1])  # angle of -u measured from +z toward +y
    pts = [(0.0, 0.0)]
    steps = 4
    for i in range(steps + 1):
        ang = back - half + 2 * half * i / steps
        pts.append((far * math.sin(ang), far * math.cos(ang)))
    return cyl_x(R, x0, x1) & prism_x(pts, x0, x1)


def journal_spans(p: Params):
    c = p.section("crankshaft")
    mw = p.L(c["main_width"]) / 2
    return [(p.main_x(j) - mw, p.main_x(j) + mw) for j in range(p.n_mains)]


def crankshaft(p: Params, with_counterweights: bool = True):
    c = p.section("crankshaft")
    main_r = p.L(c["main_diameter"]) / 2
    pin_r = p.L(c["pin_diameter"]) / 2
    pw = p.L(c["pin_width"]) / 2
    mains = journal_spans(p)
    shaft = None
    for x0, x1 in mains:
        j = cyl_x(main_r, x0, x1)
        shaft = j if shaft is None else shaft + j
    for idx in range(p.cylinders):
        xc = p.cyl_x(idx)
        py, pz = pin_yz(p, idx)
        shaft += cyl_x(pin_r, xc - pw, xc + pw, py, pz)
        shaft += web(p, idx, mains[idx][1], xc - pw, with_counterweights)
        shaft += web(p, idx, xc + pw, mains[idx + 1][0], with_counterweights)
    front, rear = mains[0][0], mains[-1][1]
    shaft += cyl_x(p.L(c["snout_diameter"]) / 2, front - p.L(c["snout_length"]), front)
    # rear stub through the block end wall, then the flywheel flange outside it
    stub_end = rear + p.L(p.section("block")["end_wall"]) + p.L(c["flange_gap"])
    shaft += cyl_x(main_r, rear, stub_end)
    shaft += cyl_x(p.L(c["flange_diameter"]) / 2, stub_end, stub_end + p.L(c["flange_length"]))
    return shaft

