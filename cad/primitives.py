"""Small build123d helpers shared by the parts (CAD units, mm)."""
from __future__ import annotations

from build123d import Align, Box, Cylinder, Pos, Rot

CENTER = (Align.CENTER, Align.CENTER, Align.CENTER)
BOTTOM = (Align.CENTER, Align.CENTER, Align.MIN)
QUARTER_TURN = 90.0


def cyl_z(r: float, z0: float, z1: float, x: float = 0.0, y: float = 0.0):
    """Cylinder with axis along z from z0 to z1."""
    return Pos(x, y, z0) * Cylinder(r, z1 - z0, align=BOTTOM)


def cyl_x(r: float, x0: float, x1: float, y: float = 0.0, z: float = 0.0):
    """Cylinder with axis along x from x0 to x1."""
    return Pos((x0 + x1) / 2, y, z) * Rot(0, QUARTER_TURN, 0) * Cylinder(r, x1 - x0, align=CENTER)


def tube_z(r_out: float, r_in: float, z0: float, z1: float, x: float = 0.0, y: float = 0.0):
    return cyl_z(r_out, z0, z1, x, y) - cyl_z(r_in, z0, z1, x, y)


def tube_x(r_out: float, r_in: float, x0: float, x1: float, y: float = 0.0, z: float = 0.0):
    return cyl_x(r_out, x0, x1, y, z) - cyl_x(r_in, x0, x1, y, z)


def box(x0: float, x1: float, y0: float, y1: float, z0: float, z1: float):
    """Axis-aligned box from its extents."""
    return Pos((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2) * Box(x1 - x0, y1 - y0, z1 - z0, align=CENTER)


def prism_x(pts, x0: float, x1: float):
    """Extrude a (y, z) polygon along +x from x0 to x1 (independent of polygon winding)."""
    from build123d import Plane, Polygon, extrude

    face = Plane.YZ.offset(x0) * Polygon(*pts, align=None)
    return extrude(face, amount=x1 - x0, dir=(1, 0, 0))
