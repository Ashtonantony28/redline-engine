"""Aluminium cylinder block with cast-in iron liners, open-deck water jacket, five bulkheads with main-bearing
saddles and separate main caps, deep skirt to the pan rail, oil gallery with feeds to every saddle, head-bolt
bosses and tapped holes, breather windows, external ribs and rounded external edges; plus the parts that belong
to it: liners, head gasket, head bolts, main-cap bolts, oil pan and pan bolts.

Frame: crank axis along x (front = -x), z up the bore axis, deck at Params.deck_z. See spec/cad.yaml `block`,
`main_cap`, `head_gasket`, `oil_pan`. All geometry MODELLED.
"""
from __future__ import annotations

from functools import lru_cache
import math

from build123d import Plane, Polygon, Pos, Rectangle, Solid, extrude, fillet, offset

from cad.conrod import big_end_bore_r
from cad.crankshaft import journal_spans
from cad.maincap import cap_width, main_bore_r, main_cap_dims
from cad.params import Params
from cad.piston import fillet_vertices
from cad.primitives import box, cyl_x, cyl_z

TOL = 1.0e-6
BIG = 1.0e3


def bsec(p: Params) -> dict:
    return p.section("block")


def block_x_extent(p: Params):
    mains = journal_spans(p)
    e = p.L(bsec(p)["end_wall"])
    return mains[0][0] - e, mains[-1][1] + e


def swing_radius(p: Params) -> float:
    """Largest radius swept by the crank (counterweights) or a big-end bolt head / boss corner."""
    c = p.section("conrod")
    be_r = big_end_bore_r(p)
    boss_y = be_r + p.L(c["bolt_boss_width"])
    big_end = p.crank_r + math.hypot(boss_y, p.L(c["bolt_boss_height"]) / 2)
    from cad.conrod import bolt_y
    hb = p.section("conrod_bolt")
    head = p.crank_r + math.hypot(bolt_y(p) + p.L(hb["head_diameter"]) / 2,
                                  p.L(c["bolt_boss_height"]) / 2 + p.L(hb["head_height"]))
    cw = p.L(p.section("crankshaft")["counterweight_radius"])
    return max(big_end, head, cw)


def rail_z(p: Params) -> float:
    return -p.L(bsec(p)["skirt_below_axis"])


def liner_r(p: Params):
    """(inner, outer) radius of the iron liner."""
    return p.bore_r, p.bore_r + p.L(bsec(p)["liner_iron_thickness"])


def barrel_r(p: Params) -> float:
    return p.bore_r + p.L(bsec(p)["liner_thickness"])


def liner_bottom_z(p: Params) -> float:
    return p.deck_z - p.L(bsec(p)["bore_length"])


def upper_half_width(p: Params) -> float:
    """Half width of the upper block / deck: the jacket wall, or wider where the head-bolt bosses need it."""
    b = bsec(p)
    jacket = barrel_r(p) + p.L(b["water_jacket_width"]) + p.L(b["jacket_wall"])
    bolts = max(abs(y) for _x, y in head_bolt_xy(p)) + p.L(b["head_bolt"]["boss_diameter"]) / 2
    return max(jacket, bolts)


def outer_profile(p: Params, grow: float = 0.0):
    """(y, z) outer section of the block (pan rail lip, skirt, flare, upper barrel section), 2D-filleted."""
    b = bsec(p)
    hw, hu = p.L(b["half_width"]), upper_half_width(p)
    lip, lip_h = p.L(b["pan_rail_lip"]), p.L(b["pan_rail_height"])
    zr, zb, zt, deck = rail_z(p), p.L(b["flare_bottom_z"]), p.L(b["flare_top_z"]), p.deck_z
    right = [(hw + lip, zr), (hw + lip, zr + lip_h), (hw, zr + lip_h), (hw, zb), (hu, zt), (hu, deck)]
    pts = right + [(-y, z) for y, z in reversed(right)]
    sk = Polygon(*pts, align=None)
    r = p.L(b["profile_fillet"])
    sk = fillet_vertices(sk, r, lambda v: abs(v.Y - zb) < TOL or abs(v.Y - zt) < TOL)
    # rail-lip re-entrant corner: radius limited by the lip width
    sk = fillet_vertices(sk, min(r, lip / 2), lambda v: abs(abs(v.X) - hw) < TOL and abs(v.Y - (zr + lip_h)) < TOL)
    return offset(sk, grow) if grow else sk


def deck_outline(p: Params):
    """Plan outline of the deck (rounded rectangle), shared with the head gasket and the head face."""
    x0, x1 = block_x_extent(p)
    hu = upper_half_width(p)
    sk = Pos((x0 + x1) / 2, 0) * Rectangle(x1 - x0, 2 * hu)
    return fillet_vertices(sk, p.L(bsec(p)["plan_corner_radius"]), lambda v: True)


def head_bolt_xy(p: Params) -> list:
    """Head bolts: one under each cam journal (x = main bulkheads, y = the cam axis of that side), so the cam-cap
    bolts beside them stay in solid pedestal material (cad/head.py)."""
    from cad.valvetrain import cam_axis_yz
    return [(p.main_x(j), cam_axis_yz(p, k)[0]) for j in range(p.n_mains) for k in ("intake", "exhaust")]


def bulkhead_span(p: Params, j: int):
    w = cap_width(p)
    return p.main_x(j) - w / 2, p.main_x(j) + w / 2


def cyl_between(p0, p1, r: float):
    import numpy as np
    a, b = np.asarray(p0, float), np.asarray(p1, float)
    d = b - a
    L = float(np.linalg.norm(d))
    return Solid.make_cylinder(r, L, Plane(origin=tuple(a), z_dir=tuple(d / L)))


def gallery_yz(p: Params):
    g = bsec(p)["gallery"]
    return p.L(g["y"]), p.L(g["z"])


@lru_cache(maxsize=None)
def block(p: Params):
    b = bsec(p)
    x0, x1 = block_x_extent(p)
    deck = p.deck_z
    zr = rail_z(p)
    body = extrude(Plane.YZ.offset(x0) * outer_profile(p), amount=x1 - x0)
    try:  # round every external end-face edge (front/rear walls)
        ends = [e for e in body.edges() if abs(e.center().X - x0) < TOL or abs(e.center().X - x1) < TOL]
        body = fillet(ends, p.L(b["end_edge_fillet"]))
    except Exception:  # cosmetic only
        pass
    # external ribs at the bulkheads (skirt and flare), head-bolt bosses, gallery boss
    rib_t = p.L(b["rib_thickness"])
    rib_band = (outer_profile(p, p.L(b["rib_height"])) - outer_profile(p)) & \
        (Pos(0, (zr + p.L(b["flare_top_z"])) / 2) * Rectangle(BIG, p.L(b["flare_top_z"]) - zr - 2 * p.L(b["pan_rail_height"])))
    for j in range(1, p.n_mains - 1):
        mx = p.main_x(j)
        body += extrude(Plane.YZ.offset(mx - rib_t / 2) * rib_band, amount=rib_t)
    hb = b["head_bolt"]
    for (x, y) in head_bolt_xy(p):
        body += cyl_z(p.L(hb["boss_diameter"]) / 2, p.L(hb["boss_bottom_z"]), deck, x, y)
    gy, gz = gallery_yz(p)
    g = b["gallery"]
    body += cyl_x(p.L(g["boss_diameter"]) / 2, x0, x1, gy, gz)
    # deck plan: rounded corners at the deck level as well
    body &= extrude(Plane.XY.offset(zr - 1.0) * _plan(p), amount=deck - zr + 2.0)

    # crankcase bays between the bulkheads, open at the bottom
    sw = swing_radius(p) + p.L(b["crank_swing_clearance"])
    inner = outer_profile(p, -p.L(b["crankcase_wall"]))
    roof = deck - p.L(b["water_jacket_depth"]) - p.L(b["bay_roof_gap"])
    inner = inner & (Pos(0, (zr - 1.0 + roof) / 2) * Rectangle(BIG, roof - zr + 1.0))
    lb = liner_bottom_z(p)
    for i in range(p.cylinders):
        bx0, bx1 = bulkhead_span(p, i)[1], bulkhead_span(p, i + 1)[0]
        bay = extrude(Plane.YZ.offset(bx0) * inner, amount=bx1 - bx0)
        bay += cyl_x(sw, bx0, bx1) & box(bx0, bx1, -BIG, BIG, zr - 1.0, roof)
        bay -= cyl_z(barrel_r(p), lb, roof + 1.0, p.cyl_x(i))
        bay -= cyl_x(p.L(g["boss_diameter"]) / 2, bx0 - 1, bx1 + 1, gy, gz)
        body -= bay
    # main saddles and the cap register: below the joint face the bulkhead is open between the cap faces
    body -= cyl_x(main_bore_r(p), x0 - 1.0, x1 + 1.0)
    d = main_cap_dims(p)
    for j in range(p.n_mains):
        cx0, cx1 = bulkhead_span(p, j)
        body -= box(cx0, cx1, -d["half_span"], d["half_span"], zr - 1.0, 0.0)
        # main-cap bolt threads (nominal diameter)
        for s in (-1, 1):
            body -= cyl_z(p.L(p.fea["main_cap"]["bolt_nominal_d"]) / 2, -1.0, p.L(p.section("main_cap")["thread_depth"]),
                          p.main_x(j), s * d["bolt_y"])
        # oil feed from the gallery to the top of the saddle
        body -= cyl_between((p.main_x(j), gy, gz), (p.main_x(j), 0.0, main_bore_r(p) - 1.0), p.L(g["feed_diameter"]) / 2)
    # breather windows in the inner bulkheads (exhaust side; the feeds run on the intake side)
    w = b["bulkhead_window"]
    win = fillet_vertices(Pos(p.L(w["y0"] + w["y1"]) / 2, p.L(w["z0"] + w["z1"]) / 2)
                          * Rectangle(p.L(w["y1"] - w["y0"]), p.L(w["z1"] - w["z0"])), p.L(w["corner_radius"]),
                          lambda v: True)
    for j in range(1, p.n_mains - 1):
        cx0, cx1 = bulkhead_span(p, j)
        body -= extrude(Plane.YZ.offset(cx0 - 1.0) * win, amount=cx1 - cx0 + 2.0)
    # bores: liner bore to the deck, skirt relief below the liner
    skirt_bdc = p.rod_len - p.crank_r - p.L(p.section("piston")["skirt_below_pin"])
    relief_bot = skirt_bdc - p.L(b["piston_bdc_clearance"])
    for i in range(p.cylinders):
        x = p.cyl_x(i)
        body -= cyl_z(liner_r(p)[1], lb, deck + 1.0, x)
        body -= cyl_z(p.bore_r + p.L(b["bore_relief_clearance"]), relief_bot, lb + 1.0, x)
    if b["water_jacket"]:
        body -= water_jacket(p)
    # oil gallery, head-bolt threads
    body -= cyl_x(p.L(g["diameter"]) / 2, x0 - 1.0, x1 + 1.0, gy, gz)
    for (x, y) in head_bolt_xy(p):
        body -= cyl_z(p.L(hb["diameter"]) / 2, deck - p.L(hb["thread_depth"]), deck + 1.0, x, y)
    # intake-manifold stay bolts in the rear wall (cad/manifolds.py bracket_bolts)
    from cad.manifolds import bracket_bolts
    ms = p.section("manifolds")
    for (x, y, z) in bracket_bolts(p)["block"]:
        body -= cyl_x(p.L(ms["bolt_diameter"]) / 2, x - p.L(ms["intake"]["bracket"]["thread_depth"]), x + 1.0, y, z)
    # pan-bolt threads in the rail lip
    for (x, y) in pan_bolt_xy(p):
        body -= cyl_z(p.L(p.section("oil_pan")["bolt_diameter"]) / 2, zr - 1.0, zr + p.L(p.section("oil_pan")["thread_depth"]),
                      x, y)
    return body


def _plan(p: Params):
    """Plan outline of the whole block: rounded rectangle at the skirt width (lip included)."""
    x0, x1 = block_x_extent(p)
    b = bsec(p)
    hw = p.L(b["half_width"]) + p.L(b["pan_rail_lip"])
    sk = Pos((x0 + x1) / 2, 0) * Rectangle(x1 - x0, 2 * hw)
    return fillet_vertices(sk, p.L(b["plan_corner_radius"]), lambda v: True)


def water_jacket(p: Params):
    """Open-deck siamesed jacket: band around the outside of the joined barrels."""
    b = bsec(p)
    deck = p.deck_z
    z0 = deck - p.L(b["water_jacket_depth"])
    r_in = barrel_r(p)
    r_out = r_in + p.L(b["water_jacket_width"])
    outer = inner = None
    for i in range(p.cylinders):
        o = cyl_z(r_out, z0, deck + 1.0, p.cyl_x(i))
        n = cyl_z(r_in, z0 - 1.0, deck + 2.0, p.cyl_x(i))
        outer = o if outer is None else outer + o
        inner = n if inner is None else inner + n
    return outer - inner


def liner(p: Params):
    """Cast-in iron liner (local: bore axis = z through the origin), deck-flush."""
    ri, ro = liner_r(p)
    return cyl_z(ro, liner_bottom_z(p), p.deck_z) - cyl_z(ri, liner_bottom_z(p) - 1.0, p.deck_z + 1.0)


# --------------------------------------------------------------------------------------------- gasket, head bolts
def head_gasket(p: Params):
    """Compressed multi-layer steel gasket between deck and head face: fire openings, bolt holes, coolant holes
    over the jacket band (spec/cad.yaml head_gasket)."""
    gk = p.section("head_gasket")
    deck = p.deck_z
    t = p.head_face_z - deck
    g = extrude(Plane.XY.offset(deck) * deck_outline(p), amount=t)
    for i in range(p.cylinders):
        g -= cyl_z(p.bore_r + p.L(gk["fire_ring_clearance"]), deck - 1, deck + t + 1, p.cyl_x(i))
    hb = bsec(p)["head_bolt"]
    for (x, y) in head_bolt_xy(p):
        g -= cyl_z(p.L(hb["diameter"]) / 2 + p.L(gk["bolt_hole_clearance"]), deck - 1, deck + t + 1, x, y)
    rc = barrel_r(p) + p.L(bsec(p)["water_jacket_width"]) / 2
    n = int(gk["coolant_holes_per_cylinder"])
    for i in range(p.cylinders):
        for k in range(n):
            a = 2 * math.pi * (k + 0.5) / n
            x, y = p.cyl_x(i) + rc * math.cos(a), rc * math.sin(a)
            # only where the jacket is open (not across the siamese webs)
            if all(math.hypot(x - p.cyl_x(m), y) >= barrel_r(p) for m in range(p.cylinders)):
                g -= cyl_z(p.L(gk["coolant_hole_diameter"]) / 2, deck - 1, deck + t + 1, x, y)
    return g


def head_bolt(p: Params):
    """Head bolt at (0, 0) in plan (instanced at each bolt position): flanged head on the head's bolt seat, shank
    through the head and gasket, threaded into the block (nominal diameter, line-to-line with the tapped hole)."""
    hb = bsec(p)["head_bolt"]
    from cad.head import head_bolt_seat_z
    seat = head_bolt_seat_z(p)
    bottom = p.deck_z - p.L(hb["thread_depth"]) + p.L(hb["thread_bottom_gap"])
    return cyl_z(p.L(hb["diameter"]) / 2, bottom, seat) + cyl_z(p.L(hb["head_diameter"]) / 2, seat,
                                                                seat + p.L(hb["head_height"]))


# --------------------------------------------------------------------------------------------- oil pan
def pan_bolt_xy(p: Params) -> list:
    op = p.section("oil_pan")
    x0, x1 = block_x_extent(p)
    hw = p.L(bsec(p)["half_width"]) + p.L(bsec(p)["pan_rail_lip"]) / 2
    n = int(op["bolts_per_side"])
    m = p.L(op["bolt_end_margin"])
    xs = [x0 + m + (x1 - x0 - 2 * m) * k / (n - 1) for k in range(n)]
    return [(x, s * hw) for x in xs for s in (-1, 1)]


def _pan_body(p: Params, inset: float, top_extra: float):
    """Pan volume (outer at inset 0): (y, z) section with rounded bottom corners intersected with the (x, z) side
    profile (shallow front, deep rear sump)."""
    op = p.section("oil_pan")
    x0, x1 = block_x_extent(p)
    zr = rail_z(p) - p.L(op["flange_thickness"])
    e = p.L(op["wall_inset"]) + inset
    hw = p.L(bsec(p)["half_width"]) - p.L(op["wall_inset"]) - inset
    zs = zr - p.L(op["shallow_depth"]) + inset
    zd = zr - p.L(op["sump_depth"]) + inset
    sx0 = p.L(op["sump_start_x"])
    slope = p.L(op["sump_slope_run"])
    top = zr + top_extra
    xz = Polygon((x0 + e, top), (x0 + e, zs), (sx0 - slope, zs), (sx0, zd), (x1 - e, zd), (x1 - e, top), align=None)
    xz = fillet_vertices(xz, p.L(op["corner_radius"]), lambda v: v.Y < top - TOL)
    yz = Polygon((-hw, top), (-hw, zd - 1.0), (hw, zd - 1.0), (hw, top), align=None)
    a = extrude(Plane.XZ * xz, amount=BIG, both=True)
    bsec_ = extrude(Plane.YZ.offset(x0 - 1.0) * yz, amount=x1 - x0 + 2.0)
    body = a & bsec_
    try:
        bot = [e_ for e_ in body.edges() if e_.center().Z < zs + 1.0 and abs(abs(e_.center().Y) - hw) < 1.0]
        body = fillet(bot, p.L(op["corner_radius"]) - inset)
    except Exception:
        pass
    return body


@lru_cache(maxsize=None)
def oil_pan(p: Params):
    op = p.section("oil_pan")
    zr = rail_z(p)
    t = p.L(op["wall"])
    ft = p.L(op["flange_thickness"])
    outer = _pan_body(p, 0.0, 0.0)
    inner = _pan_body(p, t, ft + 1.0)
    flange = extrude(Plane.XY.offset(zr - ft) * _plan(p), amount=ft)
    pan = (outer + flange) - inner
    for (x, y) in pan_bolt_xy(p):
        pan -= cyl_z(p.L(op["bolt_diameter"]) / 2 + p.L(op["bolt_hole_clearance"]), zr - ft - 1, zr + 1, x, y)
    # drain boss + plug hole at the sump floor
    return pan


def pan_bolt(p: Params):
    """Pan bolt at (0, 0) in plan, head under the pan flange, thread into the block rail lip."""
    op = p.section("oil_pan")
    zr = rail_z(p)
    ft = p.L(op["flange_thickness"])
    r = p.L(op["bolt_diameter"]) / 2
    return cyl_z(r, zr - ft, zr + p.L(op["thread_depth"]) - p.L(op["thread_bottom_gap"])) + \
        cyl_z(p.L(op["bolt_head_diameter"]) / 2, zr - ft - p.L(op["bolt_head_height"]), zr - ft)
