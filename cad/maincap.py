"""Main bearing cap (nodular iron) and its bolts: the ONE definition used by the CAD assembly and by the M6 FEA
(fea.geometry re-exports main_cap_dims / main_cap_solid), so the analysed cap is the cap in the engine.

Cap frame: crank axis along x through the origin, joint face at z = 0 (cap below), centred on x = 0; it is placed
at each main journal. The cap is narrower than the journal by main_cap.axial_clearance each side so the crank
webs run clear of the cap and bulkhead faces.
"""
from __future__ import annotations

from build123d import Pos, Rectangle, extrude

from cad.params import Params
from cad.piston import fillet_vertices
from cad.primitives import cyl_x, cyl_z

TOL = 1.0e-6


def main_bore_r(p: Params) -> float:
    return (p.L(p.section("crankshaft")["main_diameter"]) + p.L(p.section("block")["main_bore_over_journal"])) / 2


def cap_width(p: Params) -> float:
    return p.L(p.section("crankshaft")["main_width"]) - 2 * p.L(p.section("main_cap")["axial_clearance"])


def main_cap_dims(p: Params) -> dict:
    c = p.section("main_cap")
    bore_r = main_bore_r(p)
    return {"bore_r": bore_r, "depth": bore_r + p.L(c["depth_below_bore"]), "half_span": p.L(c["half_span"]),
            "bolt_y": p.L(c["bolt_spacing"]) / 2, "hole_r": p.L(c["bolt_hole_diameter"]) / 2,
            "head_r": p.L(c["head_bearing_diameter"]) / 2, "corner_r": p.L(c["corner_radius"]),
            "width": cap_width(p)}


def main_cap_solid(p: Params, half: bool = True):
    """Cap below the crank axis (joint face at z = 0), bore axis along x, bolts vertical at y = +-bolt_y."""
    from build123d import Plane
    d = main_cap_dims(p)
    sk = Pos(0, -d["depth"] / 2) * Rectangle(2 * d["half_span"], d["depth"])
    sk = fillet_vertices(sk, d["corner_r"], lambda v: v.Y < -TOL)
    cap = extrude(Plane.YZ * sk, amount=d["width"] / 2, both=True)
    cap -= cyl_x(d["bore_r"], -d["width"], d["width"])
    for s in (1, -1):
        cap -= cyl_z(d["hole_r"], -2 * d["depth"], d["depth"], 0.0, s * d["bolt_y"])
    if half:
        from cad.primitives import box
        cap = cap & box(-1e3, 1e3, 0.0, 1e3, -1e3, 1e3)
    return cap


def main_cap_bolt(p: Params):
    """Main-cap bolt at y = +bolt_y (instanced at -bolt_y by a 180 deg turn about z): flanged head bearing on the cap
    underside, shank through the cap clearance hole, threaded into the bulkhead (thread modelled at the nominal
    diameter, line-to-line with the tapped hole)."""
    d = main_cap_dims(p)
    c = p.section("main_cap")
    r = p.L(p.fea["main_cap"]["bolt_nominal_d"]) / 2
    top = p.L(c["thread_depth"])
    return (cyl_z(r, -d["depth"], top, 0.0, d["bolt_y"])
            + cyl_z(d["head_r"], -d["depth"] - p.L(c["bolt_head_height"]), -d["depth"], 0.0, d["bolt_y"]))
