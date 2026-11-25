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


def profile_face(start, segments, plane=None):
    """Closed planar face from a start point and segments: ("L", (u, v)) lines and ("A", (u, v) mid, (u, v) end)
    three-point arcs, in the plane's local (u, v) coordinates (default Plane.XZ: u = x, v = z)."""
    from build123d import Line, Plane, ThreePointArc, Wire, Face

    plane = plane or Plane.XZ
    g = lambda q: plane.from_local_coords((q[0], q[1]))
    edges, cur = [], start
    for seg in segments:
        if seg[0] == "L":
            if (abs(seg[1][0] - cur[0]) + abs(seg[1][1] - cur[1])) > 1e-9:
                edges.append(Line(g(cur), g(seg[1])))
            cur = seg[1]
        else:
            edges.append(ThreePointArc(g(cur), g(seg[1]), g(seg[2])))
            cur = seg[2]
    if (abs(start[0] - cur[0]) + abs(start[1] - cur[1])) > 1e-9:
        edges.append(Line(g(cur), g(start)))
    return Face(Wire(edges))


def revolve_rz(start, segments):
    """Solid of revolution about z from an (r, z) profile (see profile_face)."""
    from build123d import Axis, revolve

    return revolve(profile_face(start, segments), Axis.Z)


def tangent_fillet_arc(p_from, r_stem: float, z_stem: float):
    """Concave arc from p_from=(r, z) to the stem point (r_stem, z_stem), tangent to the stem there.
    Returns ("A", mid, end) for profile_face."""
    import math

    a = p_from[0] - r_stem
    b = z_stem - p_from[1]
    rho = (a * a + b * b) / (2 * a)
    cx, cz = r_stem + rho, z_stem
    a0 = math.atan2(p_from[1] - cz, p_from[0] - cx)
    a1 = math.pi                      # the stem point, straight in -r from the arc centre
    if a0 < 0:                        # take the short way round from a0 to pi
        a0 += 2 * math.pi
    am = 0.5 * (a0 + a1)
    return ("A", (cx + rho * math.cos(am), cz + rho * math.sin(am)), (r_stem, z_stem))
