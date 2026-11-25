"""DOHC 16V cylinder head: pent-roof chambers, siamesed curved ports, seat-insert and guide bores, spring pockets,
bucket towers, cam journal pedestals with tapped cap-bolt holes, spark-plug bores and head-bolt seats.

Frame: crank axis along x, z up; the head face sits on the gasket at Params.head_face_z. Intake on -y, exhaust on +y.
The roof faces are normal to their valve axes (valve angles from spec/valvetrain.yaml) and meet the head face at
the bore edge (Params.roof_z); the chamber footprint is the valve-head sweep plus a margin solved for the spec
compression ratio (Params.chamber_margin), joined by a core strip around the plug; outside it the head face is the
squish land. Valve-train stations along each valve axis come from cad/valvetrain.py `stack` (one geometry source
for the head bores and the valve-train parts). Dimensions: spec/cad.yaml `head` and `valvetrain`.

Simplifications (stated in docs/CAD-AUDIT.md): no head water jacket or oil drain passages, no cam cover, ports
are swept circles (round throat to round flange), plug threads are plain cylinders (zero-clearance fit).
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from build123d import Circle, Plane, Polygon, Spline, Vector, Wire, extrude, sweep

from cad.block import block_x_extent, head_bolt_xy, upper_half_width
from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z, prism_x
from cad import valvetrain as vt

BIG = 1.0e3


def hc(p: Params) -> dict:
    return p.section("head")


# --------------------------------------------------------------------------------------------- chamber
def footprint_face(p: Params, idx: int):
    """Plan footprint of the chamber of cylinder idx (valve sweeps + margin, core strip), clipped to the bore."""
    m = p.chamber_margin
    xc = p.cyl_x(idx)
    z0 = p.head_face_z
    face = None
    for hull in p.footprint_hulls(m):
        f = Plane.XY.offset(z0) * Polygon(*[(xc + float(u), float(v)) for u, v in hull], align=None)
        face = f if face is None else face + f
    cx, cy = p.core_half()
    face += Plane.XY.offset(z0) * Polygon((xc - cx, p.ridge_y - cy), (xc + cx, p.ridge_y - cy),
                                          (xc + cx, p.ridge_y + cy), (xc - cx, p.ridge_y + cy), align=None)
    return face


def chamber_cavity(p: Params, idx: int):
    """Pent-roof chamber volume above the head face (inside the footprint and the bore) for cylinder idx."""
    R = p.bore_r
    xc = p.cyl_x(idx)
    z0 = p.head_face_z
    H = p.roof_height
    prism = extrude(footprint_face(p, idx), amount=H + 1.0)
    roof = prism_x([(-R, z0), (R, z0), (p.ridge_y, z0 + H)], xc - R, xc + R)
    return prism & roof & cyl_z(R, z0, z0 + H + 1.0, xc)


# --------------------------------------------------------------------------------------------- layout
def side_y(p: Params, kind: str) -> float:
    """y of the port flange face on the kind's side of the head."""
    return p.side_sign(kind) * p.L(hc(p)[f"{kind}_face_y"])


def port_flange_z(p: Params, kind: str) -> float:
    return p.roof_z(p.valve_y(kind)) + p.L(hc(p)[f"{kind}_port_rise"])


def port_exit_r(p: Params, kind: str) -> float:
    g = p.gas
    d = g["intake"]["runner"]["diameter"] if kind == "intake" else g["exhaust"]["primary"]["diameter"]
    return p.L(d) / 2


def lower_top_z(p: Params) -> float:
    return p.head_face_z + p.L(hc(p)["lower_height"])


def cam_axis(p: Params, kind: str):
    return vt.cam_axis_yz(p, kind)


def journal_r(p: Params) -> float:
    c = vt.cam_c(p)
    return p.L(c["journal_diameter"]) / 2 + p.L(c["bearing_clearance"])


def head_bolt_seat_z(p: Params) -> float:
    """Head-bolt head seat: the head top sits head_bolt_journal_gap below the lowest cam journal bore (the bolts
    are under the cam journals and are reached with the camshafts out, as on many DOHC heads)."""
    hb = p.section("block")["head_bolt"]
    zj = min(cam_axis(p, k)[1] for k in ("intake", "exhaust")) - journal_r(p)
    return zj - p.L(hc(p)["head_bolt_journal_gap"]) - p.L(hb["head_height"])


def plug_xy(p: Params, idx: int):
    return p.cyl_x(idx), p.ridge_y


def plug_top_z(p: Params) -> float:
    return lower_top_z(p) + p.L(hc(p)["plug_tube_height"])


def port_bolt_xz(p: Params, kind: str, idx: int):
    """Manifold bolt positions (x, z) on the flange face for the port of cylinder idx."""
    h = hc(p)
    dx = port_exit_r(p, kind) + p.L(h["port_bolt_offset"])
    return [(p.cyl_x(idx) + s * dx, port_flange_z(p, kind)) for s in (-1, 1)]


# --------------------------------------------------------------------------------------------- ports
def port_path(p: Params, kind: str, idx: int, x_local: float, branch: int = 0):
    """Centre line of one port branch: from the seat-insert throat along the valve axis, curving out to the flange
    centre (both branches of a cylinder end on the same flange circle: siamesed port)."""
    h = hc(p)
    s0 = vt.stack(p, kind)["insert_top"] - p.L(h["port_insert_overlap"])
    P0 = Vector(*vt.axis_point(p, kind, s0, p.cyl_x(idx) + x_local))
    a = Vector(*p.valve_axis(kind))
    out = Vector(0.0, p.side_sign(kind), 0.0)
    # the two branches overrun the face by different amounts so their end faces are not coincident (robust fuse)
    Pe = Vector(p.cyl_x(idx), side_y(p, kind) + p.side_sign(kind) * p.L(h["port_flange_overrun"]) * (1 + branch),
                port_flange_z(p, kind))
    k0, k1 = h["port_tangent_scale"]
    L = (Pe - P0).length
    return P0, Pe, a, out, Wire([Spline(P0, Pe, tangents=[a, out], tangent_scalars=[k0, k1])]), L


def port_length(p: Params, kind: str) -> float:
    """Port centre-line length (mm), seat insert top to the flange face (outer branch of cylinder 1)."""
    lay = [lx for k, _s, lx, _y, _d in p.valve_layout() if k == kind]
    return float(port_path(p, kind, 0, lay[0])[4].length) - p.L(hc(p)["port_flange_overrun"])


def port_cut(p: Params, kind: str, idx: int):
    rt = hc(p)["throat_fraction"]
    branches = []
    for k2, _s, xl, _y, d in p.valve_layout():
        if k2 != kind:
            continue
        P0, Pe, a, o, path, _L = port_path(p, kind, idx, xl, len(branches))
        c0 = Plane(origin=P0, z_dir=a) * Circle(rt * d / 2)
        c1 = Plane(origin=Pe, z_dir=o) * Circle(port_exit_r(p, kind))
        branches.append(sweep([c0, c1], path=path, multisection=True))
    return branches[0].fuse(*branches[1:])


# --------------------------------------------------------------------------------------------- body
def _tower_profile(p: Params, kind: str):
    """(y, z) parallelogram around the bucket axis line from below the lower body top to the tower top."""
    st = vt.stack(p, kind)
    w = vt.bucket_r(p) + p.L(hc(p)["tower_wall"])
    a = p.valve_axis(kind)
    n = np.array([a[2], -a[1]])                       # in-plane normal to the axis (y, z)
    q0 = vt.axis_point(p, kind, st["floor"] - p.L(hc(p)["tower_root"]))[1:]
    q1 = vt.axis_point(p, kind, st["tower_top"])[1:]
    pts = [q0 + w * n, q1 + w * n, q1 - w * n, q0 - w * n]
    return [(float(u), float(v)) for u, v in pts]


def _lower_profile(p: Params):
    h = hc(p)
    z0 = p.head_face_z
    hd = upper_half_width(p)
    zl = z0 + p.L(h["deck_lip"])
    zf = zl + p.L(h["side_flare_height"])
    zt = lower_top_z(p)
    yi, ye = side_y(p, "intake"), side_y(p, "exhaust")
    return [(-hd, z0), (hd, z0), (hd, zl), (ye, zf), (ye, zt), (yi, zt), (yi, zf), (-hd, zl)]


def valve_cuts(p: Params, idx: int, kind: str, x_local: float):
    """Insert counterbore, guide bore, spring pocket and bucket bore along one valve axis (valve frame)."""
    st = vt.stack(p, kind)
    c = vt.vc(p)
    h = hc(p)
    cut = cyl_z(vt.insert_od_r(p, kind), -p.L(h["insert_bore_lead"]), st["insert_top"])
    cut += cyl_z(p.L(c["guide"]["outer_diameter"]) / 2, st["guide_bottom"] - 1.0, st["floor"] + 1.0)
    r_sp = p.L(c["spring_seat"]["outer_diameter"]) / 2 + p.L(c["spring_seat"]["radial_clearance"])
    bb = st["bucket_bottom"] - p.max_lift(kind) - p.L(h["bucket_bore_runout"])
    cut += cyl_z(r_sp, st["floor"], bb + 1.0)
    cut += cyl_z(vt.bucket_r(p) + p.L(h["bucket_clearance"]), bb, st["tower_top"] + BIG / 10)
    return vt.valve_frame(p, idx, kind, x_local) * cut


def lobe_relief(p: Params, kind: str, x0: float, x1: float):
    """Cylinder around the cam axis swept by the lobes (base circle + lift + clearance) over the head length; the
    journal pedestals are added after it is cut."""
    cy, cz = cam_axis(p, kind)
    r = p.L(p.vt["cam"]["base_circle_radius"]) + p.max_lift(kind) + p.L(hc(p)["lobe_clearance"])
    return cyl_x(r, x0 - 1.0, x1 + 1.0, cy, cz)


def pedestals(p: Params, kind: str):
    cc = vt.vc(p)["cam_cap"]
    c = vt.cam_c(p)
    cy, cz = cam_axis(p, kind)
    hw = p.L(cc["half_width"])
    jw = p.L(c["journal_width"])
    out = None
    for x in vt.journal_xs(p):
        b = box(x - jw / 2, x + jw / 2, cy - hw, cy + hw, lower_top_z(p) - 1.0, cz)
        out = b if out is None else out + b
    return out


@lru_cache(maxsize=None)
def head(p: Params):
    h = hc(p)
    x0, x1 = block_x_extent(p)
    body = prism_x(_lower_profile(p), x0, x1)
    body = body.fuse(*[prism_x(_tower_profile(p, k), x0, x1) for k in ("intake", "exhaust")])
    body = body.cut(*[lobe_relief(p, k, x0, x1) for k in ("intake", "exhaust")])
    adds = [pedestals(p, k) for k in ("intake", "exhaust")]
    for i in range(p.cylinders):
        x, y = plug_xy(p, i)
        adds.append(cyl_z(p.L(h["plug_boss_diameter"]) / 2, lower_top_z(p) - 1.0, plug_top_z(p), x, y))
    body = body.fuse(*adds)
    body &= extrude(Plane.XY.offset(p.head_face_z - 1.0) * _plan(p), amount=BIG)

    # cuts (one boolean): chambers, valve bores, ports, plug bores, cam journals, bolt holes
    z0 = p.head_face_z
    tools = []
    for i in range(p.cylinders):
        tools.append(chamber_cavity(p, i))
        for kind, _s, xl, _y, _d in p.valve_layout():
            tools.append(valve_cuts(p, i, kind, xl))
        for kind in ("intake", "exhaust"):
            tools.append(port_cut(p, kind, i))
            for (bx, bz) in port_bolt_xz(p, kind, i):
                tools.append(_port_bolt_hole(p, kind, bx, bz))
        x, y = plug_xy(p, i)
        zr = z0 + p.roof_height
        # the thread breaks through the sloping roof down past the plug tip (no head material under the thread end)
        z_tip = zr - p.L(h["spark_plug"]["tip_protrusion"])
        tools.append(cyl_z(p.L(h["plug_thread_diameter"]) / 2, z_tip - 1.0, zr + p.L(h["plug_thread_length"]), x, y))
        tools.append(cyl_z(p.L(h["plug_well_diameter"]) / 2, zr + p.L(h["plug_thread_length"]), plug_top_z(p) + 1.0,
                           x, y))
    cc = vt.vc(p)["cam_cap"]
    for kind in ("intake", "exhaust"):
        cy, cz = cam_axis(p, kind)
        tools.append(cyl_x(journal_r(p), x0 - 1.0, x1 + 1.0, cy, cz))
        for x in vt.journal_xs(p):
            for s in (-1, 1):
                tools.append(cyl_z(p.L(cc["bolt_diameter"]) / 2,
                                   cz - p.L(cc["thread_depth"]) - p.L(cc["thread_bottom_gap"]), cz + 1.0, x,
                                   cy + s * p.L(cc["bolt_offset"])))
    hb = p.section("block")["head_bolt"]
    seat = head_bolt_seat_z(p)
    for (x, y) in head_bolt_xy(p):
        tools.append(cyl_z(p.L(hb["diameter"]) / 2 + p.L(h["head_bolt_clearance"]), z0 - 1.0, seat + 1.0, x, y))
        # closed pocket: its roof stays head_bolt_pocket_roof under the journal bore (the journal saddle closes over
        # the bolt, as with a separate lower cam holder; modelled as one piece with the head)
        zroof = min(cam_axis(p, k)[1] for k in ("intake", "exhaust")) - journal_r(p) - p.L(h["head_bolt_pocket_roof"])
        tools.append(cyl_z(p.L(hb["head_diameter"]) / 2 + p.L(h["head_bolt_pocket_clearance"]), seat, zroof, x, y))
    return body.cut(*tools)


def _plan(p: Params):
    """Head plan: rounded rectangle from the intake to the exhaust port face (the deck-width lip and the side
    flare are in the section profile)."""
    from build123d import Pos, Rectangle
    from cad.piston import fillet_vertices
    x0, x1 = block_x_extent(p)
    yi, ye = side_y(p, "intake"), side_y(p, "exhaust")
    sk = Pos((x0 + x1) / 2, (yi + ye) / 2) * Rectangle(x1 - x0, ye - yi)
    return fillet_vertices(sk, p.L(p.section("block")["plan_corner_radius"]), lambda v: True)


def _port_bolt_hole(p: Params, kind: str, x: float, z: float):
    mf = p.section("manifolds")
    y = side_y(p, kind)
    s = p.side_sign(kind)
    depth = p.L(mf["bolt_thread_depth"]) + p.L(mf["bolt_thread_bottom_gap"])
    r = p.L(mf["bolt_diameter"]) / 2
    return _cyl_y(r, y - s * depth, y + s * 1.0, x, z)


def _cyl_y(r: float, y0: float, y1: float, x: float, z: float):
    from build123d import Pos, Rot, Cylinder
    from cad.primitives import CENTER
    lo, hi = min(y0, y1), max(y0, y1)
    return Pos(x, (lo + hi) / 2, z) * Rot(90, 0, 0) * Cylinder(r, hi - lo, align=CENTER)


# --------------------------------------------------------------------------------------------- spark plug
def spark_plug(p: Params):
    """Spark plug at the origin of its bore axis (vertical), seat on the bore step: thread, seat, hex, insulator."""
    h = hc(p)
    sp = h["spark_plug"]
    zr = p.head_face_z + p.roof_height
    rt = p.L(h["plug_thread_diameter"]) / 2
    zt = zr + p.L(h["plug_thread_length"])
    plug = cyl_z(rt, zr - p.L(sp["tip_protrusion"]), zt)
    plug += cyl_z(p.L(sp["electrode_diameter"]) / 2, zr - p.L(sp["tip_protrusion"]) - p.L(sp["electrode_length"]),
                  zr - p.L(sp["tip_protrusion"]))
    plug += cyl_z(p.L(sp["seat_diameter"]) / 2, zt, zt + p.L(sp["seat_height"]))
    zh = zt + p.L(sp["seat_height"])
    af = p.L(sp["hex_across_flats"]) / 2
    rv = af / math.cos(math.pi / 6)
    hexa = extrude(Plane.XY.offset(zh) * Polygon(*[(rv * math.cos(math.pi / 3 * k), rv * math.sin(math.pi / 3 * k))
                                                    for k in range(6)], align=None), amount=p.L(sp["hex_height"]))
    plug += hexa
    zi = zh + p.L(sp["hex_height"])
    plug += cyl_z(p.L(sp["insulator_diameter"]) / 2, zi, zi + p.L(sp["insulator_length"]))
    plug += cyl_z(p.L(sp["terminal_diameter"]) / 2, zi + p.L(sp["insulator_length"]),
                  zi + p.L(sp["insulator_length"]) + p.L(sp["terminal_length"]))
    return plug
