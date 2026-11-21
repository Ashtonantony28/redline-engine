"""Simplified cylinder head: slab with pent-roof chambers and valve seats/port throats.

Frame: crank axis along x, z up; the head face sits on the gasket at Params.head_face_z.
The pent-roof ridge runs along x through each bore axis; intake on -y, exhaust on +y, each
valve axis normal to its roof face. The roof angle is solved from the spec compression
ratio when cad.yaml says `roof_angle_deg: auto` (Params.roof_angle_deg). Valve heads sit
flush with the roof when closed, so seats/throats are not part of the clearance volume.
"""
from __future__ import annotations

import math

from build123d import Cylinder, Pos, Rot

from cad.block import block_x_extent
from cad.params import Params
from cad.primitives import BOTTOM, box, cyl_z, prism_x


def chamber_cavity(p: Params, idx: int):
    """Pent-roof combustion chamber volume above the head face for cylinder idx."""
    R = p.bore_r
    x = p.cyl_x(idx)
    z0 = p.head_face_z
    roof = prism_x([(-R, z0), (R, z0), (0.0, z0 + p.roof_height)], x - R, x + R)
    return cyl_z(R, z0, z0 + p.roof_height, x) & roof


def roof_z(p: Params, y: float) -> float:
    return p.head_face_z + (p.bore_r - abs(y)) * math.tan(math.radians(p.roof_angle_deg))


def valve_cuts(p: Params, idx: int):
    """Seat counterbores and port throats along each valve axis (normal to its roof face)."""
    h = p.section("head")
    theta = p.roof_angle_deg
    out = None
    for _kind, sign, x, y, d in p.valve_layout():
        place = Pos(p.cyl_x(idx) + x, y, roof_z(p, y)) * Rot(-sign * theta, 0, 0)
        cut = place * Cylinder(d / 2, p.L(h["seat_depth"]), align=BOTTOM)
        cut += place * Cylinder(h["throat_fraction"] * d / 2, p.L(h["port_depth"]), align=BOTTOM)
        out = cut if out is None else out + cut
    return out


def head(p: Params):
    h = p.section("head")
    x0, x1 = block_x_extent(p)
    hw = p.L(p.section("block")["half_width"])
    z0 = p.head_face_z
    body = box(x0, x1, -hw, hw, z0, z0 + p.L(h["height"]))
    for i in range(p.cylinders):
        body -= chamber_cavity(p, i)
        body -= valve_cuts(p, i)
    return body
