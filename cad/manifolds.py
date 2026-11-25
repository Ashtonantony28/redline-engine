"""Intake manifold (plenum + 4 runners + common flange) and 4-2-1 exhaust header (primaries, 2-into-1 merges,
secondaries, collector, outlet flange), bolted to the head port faces.

Sizes from spec/gasdyn.yaml: runner diameter 36 mm and plenum-to-valve-seat length 550 mm (the CAD runner length is
solved so head port + runner = 550 mm), plenum volume 4.0 L; exhaust primaries 38 mm x 450 mm (pairs 1+4 and 2+3,
lengths solved equal, measured from the head flange face), secondaries 45 mm x 450 mm, outlet 55 mm (tailpipe_front
diameter). Routing and wall/flange sizes: spec/cad.yaml `manifolds`. Pipes are swept circles along splines; each
duct is built as an outer solid minus a bore solid so walls are uniform.
"""
from __future__ import annotations

from functools import lru_cache
import math

from build123d import Circle, Plane, Spline, Vector, Wire, sweep

from cad.head import port_bolt_xz, port_exit_r, port_flange_z, port_length, side_y, _cyl_y
from cad.params import Params
from cad.primitives import box, cyl_x

SOLVE_STEPS = 40


def mc(p: Params) -> dict:
    return p.section("manifolds")


def _pipe(points, tangents, r, scal=None):
    pts = [Vector(*q) for q in points]
    sp = Spline(*pts, tangents=[Vector(*t) for t in tangents], tangent_scalars=scal) if scal else \
        Spline(*pts, tangents=[Vector(*t) for t in tangents])
    path = Wire([sp])
    t0 = Vector(*tangents[0]).normalized()
    sec = Plane(origin=pts[0], z_dir=t0) * Circle(r)
    return sweep(sec, path=path), path.length


def _tube(points, tangents, r_in, wall, scal=None):
    outer, L = _pipe(points, tangents, r_in + wall, scal)
    inner, _ = _pipe(points, tangents, r_in, scal)
    return outer, inner, L


# --------------------------------------------------------------------------------------------- intake
# Swan-neck runners: out of the flange (-y), a 90-deg bend up, a vertical leg, a 90-deg bend back (+y) over the
# head and a horizontal leg through the plenum's runner-side wall to a bellmouth lip inside. Centre lines are
# straight legs + tangent bends (FilletPolyline: cylinders + tori, robust booleans). The vertical leg is solved so
# head port + flange + runner (to the bellmouth lip) = the gasdyn runner length exactly.
def _im(p: Params) -> dict:
    return mc(p)["intake"]


def runner_cad_length(p: Params) -> float:
    """Runner centre-line length from the flange outer face to the bellmouth lip."""
    return p.L(p.gas["intake"]["runner"]["length"]) - port_length(p, "intake") - p.L(mc(p)["flange_thickness"])


def _runner_y(p: Params) -> dict:
    m = _im(p)
    y0 = side_y(p, "intake") - p.L(mc(p)["flange_thickness"])     # flange outer face
    R = p.L(m["curl_radius"])
    yc = y0 - p.L(m["runner_out"]) - R                               # corner of the two bends
    yw = p.L(m["plenum_wall_y"])
    return {"y0": y0, "R": R, "yc": yc, "yw": yw, "yend": yw + p.L(m["bellmouth_insert"])}


@lru_cache(maxsize=None)
def runner_vertical(p: Params) -> float:
    """Straight vertical leg between the two bends (closed form from the target length)."""
    g, m = _runner_y(p), _im(p)
    fixed = p.L(m["runner_out"]) + math.pi * g["R"] + (g["yend"] - (g["yc"] + g["R"]))
    s = runner_cad_length(p) - fixed
    if s <= 0.0:
        raise ValueError(f"intake runner layout already {fixed:.1f} mm > {runner_cad_length(p):.1f} mm: shorten the curl")
    return s


def runner_top_z(p: Params) -> float:
    return port_flange_z(p, "intake") + 2.0 * _runner_y(p)["R"] + runner_vertical(p)


def _runner_points(p: Params, idx: int, y_end: float):
    g = _runner_y(p)
    x, zf, zt = p.cyl_x(idx), port_flange_z(p, "intake"), runner_top_z(p)
    return [(x, g["y0"], zf), (x, g["yc"], zf), (x, g["yc"], zt), (x, y_end, zt)]


def runner_length(p: Params, idx: int = 0) -> float:
    """Measured CAD centre-line length, flange outer face to bellmouth lip (checks the closed form)."""
    g = _runner_y(p)
    return _poly_pipe(_runner_points(p, idx, g["yend"]), 1.0, g["R"])[1]


def plenum_geom(p: Params) -> dict:
    """Inner plenum box: x0..x1, y0..y1, z0..z1 (CAD mm) and inner edge radius; height solved from the volume."""
    m = _im(p)
    r = port_exit_r(p, "intake")
    x0 = p.cyl_x(0) - r - p.L(m["plenum_end_margin"])
    x1 = p.cyl_x(p.cylinders - 1) + r + p.L(m["plenum_end_margin"])
    w, rf = p.L(m["plenum_width"]), p.L(m["plenum_fillet"])
    V = p.gas["intake"]["plenum"]["volume"] * p.k**3
    # box with every edge rounded = core box (a, b, c) grown by rf (Steiner):
    # V = abc + 2 rf (ab + bc + ca) + pi rf^2 (a + b + c) + 4/3 pi rf^3, linear in the core height c
    a, b = (x1 - x0) - 2 * rf, w - 2 * rf
    c = (V - 2 * rf * a * b - math.pi * rf**2 * (a + b) - 4.0 / 3.0 * math.pi * rf**3) / \
        (a * b + 2 * rf * (a + b) + math.pi * rf**2)
    h = c + 2 * rf
    zc = runner_top_z(p)
    y0 = p.L(m["plenum_wall_y"])
    return {"x0": x0, "x1": x1, "y0": y0, "y1": y0 + w, "z0": zc - h / 2, "z1": zc + h / 2, "rf": rf,
            "yc": y0 + w / 2, "zc": zc}


def _rbox(x0, x1, y0, y1, z0, z1, rf):
    from build123d import fillet
    b = box(x0, x1, y0, y1, z0, z1)
    return fillet(b.edges(), rf)


def _bell(p: Params, x: float, outer: bool):
    """Conical bellmouth inside the plenum, bore radius at the runner-side wall flaring by bellmouth_flare at the
    lip (outer solid, or the bore to subtract)."""
    from build123d import Cone, Pos, Rot
    from cad.primitives import BOTTOM
    g, m = _runner_y(p), _im(p)
    r = port_exit_r(p, "intake")
    fl, bw = p.L(m["bellmouth_flare"]), p.L(m["bellmouth_wall"])
    Lb = g["yend"] - g["yw"]
    def rad(y):
        return r + fl * (y - g["yw"]) / Lb
    if outer:
        ya, yb, add = g["yw"] - p.L(mc(p)["wall"]) / 2, g["yend"], bw
    else:
        ya, yb, add = g["yw"] - 1.0, g["yend"] + 1.0, 0.0
    return Pos(x, ya, runner_top_z(p)) * Rot(-90, 0, 0) * Cone(rad(ya) + add, rad(yb) + add, yb - ya, align=BOTTOM)


def throttle_x(p: Params) -> dict:
    t = _im(p)["throttle"]
    pg = plenum_geom(p)
    wall = p.L(mc(p)["wall"])
    xn = pg["x1"] + wall + p.L(t["neck_length"])           # neck flange face = throttle body flange face
    return {"face": xn, "inlet": xn + p.L(t["body_length"]), "rb": p.L(p.gas["intake"]["inlet_duct"]["diameter"]) / 2}


def bracket_geom(p: Params) -> dict:
    from cad.block import block_x_extent
    b = _im(p)["bracket"]
    pg = plenum_geom(p)
    wall = p.L(mc(p)["wall"])
    xr = block_x_extent(p)[1]                               # block (and head) rear face
    t = p.L(b["thickness"])
    zpad1 = pg["z0"] - wall                                  # plenum outer floor
    zpad0 = zpad1 - p.L(b["pad_height"])
    xp1 = pg["x1"] - pg["rf"]                                # pad on the flat floor (inside the edge fillet)
    xp0 = xp1 - p.L(b["pad_length"])
    hw = p.L(b["width"]) / 2
    return {"xr": xr, "t": t, "y0": pg["yc"] - hw, "y1": pg["yc"] + hw, "yc": pg["yc"], "zb": p.L(b["bottom_z"]),
            "zpad0": zpad0, "zpad1": zpad1, "xp0": xp0, "xp1": xp1, "ztab0": zpad0 - t}


def bracket_bolts(p: Params) -> dict:
    """Placements (x, y, z) of the stay bolts: 'block' (along -x into the rear wall), 'pad' (along +z into the pad)."""
    b, g = _im(p)["bracket"], bracket_geom(p)
    dy = p.L(b["pad_bolt_dy"])
    return {"block": [(g["xr"], g["yc"], p.L(z)) for z in b["block_bolt_z"]],
            "pad": [(0.5 * (g["xp0"] + g["xp1"]), g["yc"] + s * dy, g["zpad0"]) for s in (-1, 1)]}


@lru_cache(maxsize=None)
def intake_manifold(p: Params):
    m, mi = mc(p), _im(p)
    wall = p.L(m["wall"])
    r = port_exit_r(p, "intake")
    g = _runner_y(p)
    pg = plenum_geom(p)
    outer = [_rbox(pg["x0"] - wall, pg["x1"] + wall, pg["y0"] - wall, pg["y1"] + wall, pg["z0"] - wall,
                   pg["z1"] + wall, pg["rf"] + wall)]
    inner = [_rbox(pg["x0"], pg["x1"], pg["y0"], pg["y1"], pg["z0"], pg["z1"], pg["rf"])]
    for i in range(p.cylinders):
        # outer tube ends half-way into the plenum wall, the bore runs through it
        outer.append(_poly_pipe(_runner_points(p, i, g["yw"] - wall / 2), r + wall, g["R"])[0])
        bore_pts = _runner_points(p, i, g["yw"] + 1.0)
        bore_pts[0] = (bore_pts[0][0], bore_pts[0][1] + 1.0, bore_pts[0][2])
        inner.append(_poly_pipe(bore_pts, r, g["R"])[0])
    # throttle neck with its flange on the rear end wall
    tx, t = throttle_x(p), mi["throttle"]
    rb, ft = tx["rb"], p.L(t["flange_thickness"])
    outer.append(cyl_x(rb + wall, pg["x1"] - 1.0, tx["face"], pg["yc"], pg["zc"]))
    outer.append(cyl_x(rb + p.L(t["flange_over"]), tx["face"] - ft, tx["face"], pg["yc"], pg["zc"]))
    inner.append(cyl_x(rb, pg["x1"] - 1.0, tx["face"] + 1.0, pg["yc"], pg["zc"]))
    # stay pad on the plenum floor
    bg = bracket_geom(p)
    outer.append(box(bg["xp0"], bg["xp1"], bg["y0"], bg["y1"], bg["zpad0"], bg["zpad1"] + wall / 2))
    # common flange plate on the head intake face
    yf = side_y(p, "intake")
    tf = p.L(m["flange_thickness"])
    hz = p.L(m["flange_half_height"])
    zf = port_flange_z(p, "intake")
    xa = p.cyl_x(0) - p.L(m["flange_end_margin"])
    xb = p.cyl_x(p.cylinders - 1) + p.L(m["flange_end_margin"])
    outer.append(box(xa, xb, yf - tf, yf, zf - hz, zf + hz))
    body = outer[0].fuse(*outer[1:]).clean()
    body = body.cut(*inner)
    # bellmouths inside the plenum
    bells_o = [_bell(p, p.cyl_x(i), True) for i in range(p.cylinders)]
    bells_i = [_bell(p, p.cyl_x(i), False) for i in range(p.cylinders)]
    body = body.fuse(*bells_o).cut(*bells_i)
    tools = []
    for i in range(p.cylinders):
        tools.append(_cyl_y(r, yf - tf - 1.0, yf + 1.0, p.cyl_x(i), zf))
        for (bx, bz) in port_bolt_xz(p, "intake", i):
            tools.append(_cyl_y(p.L(m["bolt_diameter"]) / 2 + p.L(m["bolt_hole_clearance"]), yf - tf - 1.0, yf + 1.0,
                                bx, bz))
    b = mi["bracket"]
    for (x, y, z) in bracket_bolts(p)["pad"]:
        tools.append(_zcyl(p.L(m["bolt_diameter"]) / 2, z - 1.0, z + p.L(b["pad_thread_depth"]), x, y))
    return body.cut(*tools)


@lru_cache(maxsize=None)
def throttle_body(p: Params):
    """Throttle body: flanged tube on the plenum neck, shaft across the bore, blade at blade_open_deg."""
    from build123d import Cylinder, Pos, Rot
    t = _im(p)["throttle"]
    pg, tx = plenum_geom(p), throttle_x(p)
    rb, w, ft = tx["rb"], p.L(t["wall"]), p.L(t["flange_thickness"])
    rfl = rb + p.L(t["flange_over"])
    y, z, x0, x1 = pg["yc"], pg["zc"], tx["face"], tx["inlet"]
    xs = 0.5 * (x0 + x1)
    body = cyl_x(rb + w, x0, x1, y, z).fuse(cyl_x(rfl, x0, x0 + ft, y, z), cyl_x(rfl, x1 - ft, x1, y, z))
    body = body.cut(cyl_x(rb, x0 - 1.0, x1 + 1.0, y, z))
    rs = p.L(t["shaft_diameter"]) / 2
    hl = rb + w + p.L(t["shaft_over"])
    shaft = _cyl_y(rs, y - hl, y + hl, xs, z)
    blade = Pos(xs, y, z) * Rot(0, 90.0 - float(t["blade_open_deg"]), 0) * Cylinder(rb, p.L(t["blade_thickness"]))
    blade = blade & cyl_x(rb, x0, x1, y, z)
    return body.fuse(shaft, blade).clean()


@lru_cache(maxsize=None)
def intake_bracket(p: Params):
    """Steel stay: vertical plate on the block rear face (and the head/gasket above it), tab under the plenum pad."""
    m, b, g = mc(p), _im(p)["bracket"], bracket_geom(p)
    plate = box(g["xr"], g["xr"] + g["t"], g["y0"], g["y1"], g["zb"], g["zpad0"])
    tab = box(g["xp0"], g["xr"] + g["t"], g["y0"], g["y1"], g["ztab0"], g["zpad0"])
    body = plate.fuse(tab).clean()
    hole = p.L(m["bolt_diameter"]) / 2 + p.L(m["bolt_hole_clearance"])
    tools = [cyl_x(hole, g["xr"] - 1.0, g["xr"] + g["t"] + 1.0, y, z) for (_x, y, z) in bracket_bolts(p)["block"]]
    tools += [_zcyl(hole, g["ztab0"] - 1.0, g["zpad0"] + 1.0, x, y) for (x, y, _z) in bracket_bolts(p)["pad"]]
    return body.cut(*tools)


def bracket_bolt(p: Params, kind: str):
    """Stay bolt in its local frame (origin on the clamped face of the threaded part): 'block' along -x, 'pad' +z."""
    m, b, g = mc(p), _im(p)["bracket"], bracket_geom(p)
    r = p.L(m["bolt_diameter"]) / 2
    rh, hh = p.L(m["bolt_head_diameter"]) / 2, p.L(m["bolt_head_height"])
    if kind == "block":
        return cyl_x(r, -p.L(b["thread_depth"]), g["t"]) + cyl_x(rh, g["t"], g["t"] + hh)
    return _zcyl(r, -g["t"], p.L(b["pad_thread_depth"]), 0.0, 0.0) + _zcyl(rh, -g["t"] - hh, -g["t"], 0.0, 0.0)


def _zcyl(r, z0, z1, x, y):
    from cad.primitives import cyl_z
    return cyl_z(r, z0, z1, x, y)


# --------------------------------------------------------------------------------------------- exhaust
# Exhaust ducts are swept along straight legs joined by tangent bends of radius bend_radius (cylinders + tori:
# analytic surfaces, so the 2-into-1 unions are robust). Where two pipes enter a merge cone or the collector they
# run parallel, side by side (offset +-merge_offset / +-collector_mouth_over), into a mouth wider than both.
def _poly_pipe(pts, r, bend):
    from build123d import FilletPolyline
    path = FilletPolyline(*[Vector(*q) for q in pts], radius=bend)
    t0 = (Vector(*pts[1]) - Vector(*pts[0])).normalized()
    return sweep(Plane(origin=Vector(*pts[0]), z_dir=t0) * Circle(r), path=path), path.length


def _merge_centre(p: Params, pair: int):
    """Merge-cone mouth centre of a pair (pairs 1+4 and 2+3: 360 deg apart in firing order 1-3-4-2). Pair 1+4 (outer
    cylinders) merges further outboard, pair 2+3 deeper, so the two sets of primaries never cross."""
    me = mc(p)["exhaust"]
    yf = side_y(p, "exhaust") + p.L(mc(p)["flange_thickness"])
    zf = port_flange_z(p, "exhaust")
    return 0.0, yf + p.L(me["merge_y"]) + (1 - pair) * p.L(me["merge_pair_dy"]), \
        zf - p.L(me["merge_drop"]) - pair * p.L(me["inner_pair_extra_drop"])


def _pair_side(p: Params, idx: int):
    return (0 if idx in (0, p.cylinders - 1) else 1), (-1.0 if p.cyl_x(idx) < 0 else 1.0)


def _primary_points(p: Params, idx: int, h: float):
    """Flange -> out -> down by h (solved) -> diagonal -> vertical entry into the merge mouth."""
    me = mc(p)["exhaust"]
    yf = side_y(p, "exhaust") + p.L(mc(p)["flange_thickness"])
    zf = port_flange_z(p, "exhaust")
    pair, side = _pair_side(p, idx)
    xm, ym, zm = _merge_centre(p, pair)
    x = p.cyl_x(idx)
    a = p.L(me["primary_out"])
    xe = xm + side * p.L(me["merge_offset"])
    z3 = zm + p.L(me["entry_length"])
    return [(x, yf, zf), (x, yf + a, zf), (x, yf + a, zf - h), (xe, ym, z3), (xe, ym, zm)]


def _bend(p: Params) -> float:
    return p.L(mc(p)["exhaust"]["bend_radius"])


@lru_cache(maxsize=None)
def primary_drop(p: Params, idx: int) -> float:
    """Vertical drop h of primary idx solved for the spec primary length (length grows monotonically with h)."""
    target = p.L(p.gas["exhaust"]["primary"]["length"])
    r = port_exit_r(p, "exhaust")
    me = mc(p)["exhaust"]
    pair, _ = _pair_side(p, idx)
    zf = port_flange_z(p, "exhaust")
    z3 = _merge_centre(p, pair)[2] + p.L(me["entry_length"])
    lo, hi = 2.0 * _bend(p), zf - z3 - 2.0 * _bend(p)
    L = lambda h: _poly_pipe(_primary_points(p, idx, h), r, _bend(p))[1]
    if not L(lo) <= target <= L(hi):
        raise ValueError(f"primary {idx + 1}: spec length {target:.0f} mm outside the routable range "
                         f"{L(lo):.0f}..{L(hi):.0f} mm (adjust manifolds.exhaust)")
    for _ in range(SOLVE_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if L(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def primary_length(p: Params, idx: int) -> float:
    return _poly_pipe(_primary_points(p, idx, primary_drop(p, idx)), port_exit_r(p, "exhaust"), _bend(p))[1]


def _collector_centre(p: Params):
    """Collector mouth: below the pair merges, at the y midway between them; depth so each secondary (vertical,
    diagonal, vertical entry) has the spec secondary length."""
    me = mc(p)["exhaust"]
    _, y0, z0 = _merge_centre(p, 0)
    _, y1, z1 = _merge_centre(p, 1)
    yc = 0.5 * (y0 + y1)
    return 0.0, yc, min(z0, z1) - p.L(me["merge_cone_length"]) - p.L(me["secondary_drop"])


def _secondary_points(p: Params, pair: int, w: float):
    """Merge-cone outlet -> short vertical -> swing out sideways by w in x (solved: the two secondaries start at
    different heights but must have the same spec length) -> vertical entry into the collector mouth."""
    me = mc(p)["exhaust"]
    xm, ym, zm = _merge_centre(p, pair)
    xc, yc, zc = _collector_centre(p)
    side = 1.0 if ym > yc else -1.0
    ye = yc + side * p.L(me["collector_mouth_over"])
    z0 = zm - p.L(me["merge_cone_length"])
    z1 = z0 - p.L(me["secondary_vertical"])
    z3 = zc + p.L(me["entry_length"])
    sx = -1.0 if pair == 0 else 1.0
    return [(xm, ym, z0), (xm, ym, z1), (xm + sx * w, 0.5 * (ym + ye), 0.5 * (z1 + z3)), (xc, ye, z3), (xc, ye, zc)]


@lru_cache(maxsize=None)
def secondary_h(p: Params, pair: int) -> float:
    """Sideways swing w of a secondary solved for the spec secondary length (length grows with w)."""
    target = p.L(p.gas["exhaust"]["secondary"]["length"])
    r = p.L(p.gas["exhaust"]["secondary"]["diameter"]) / 2
    lo, hi = p.L(mc(p)["exhaust"]["min_secondary_swing"]), p.L(mc(p)["exhaust"]["max_secondary_swing"])
    L = lambda w: _poly_pipe(_secondary_points(p, pair, w), r, _bend(p))[1]
    if not L(lo) <= target <= L(hi):
        raise ValueError(f"secondary {pair}: spec length {target:.0f} mm outside {L(lo):.0f}..{L(hi):.0f} mm")
    for _ in range(SOLVE_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if L(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


@lru_cache(maxsize=None)
def exhaust_header(p: Params):
    m = mc(p)
    me = m["exhaust"]
    wall = p.L(m["wall"])
    bend = _bend(p)
    r1 = port_exit_r(p, "exhaust")
    r2 = p.L(p.gas["exhaust"]["secondary"]["diameter"]) / 2
    r3 = p.L(p.gas["exhaust"]["tailpipe_front"]["diameter"]) / 2
    outer, inner = [], []
    for i in range(p.cylinders):
        pts = _primary_points(p, i, primary_drop(p, i))
        outer.append(_poly_pipe(pts, r1 + wall, bend)[0])
        inner.append(_poly_pipe(pts, r1, bend)[0])
    Lc = p.L(me["merge_cone_length"])
    for pair in (0, 1):
        xm, ym, zm = _merge_centre(p, pair)
        rm = r1 + p.L(me["merge_offset"]) + p.L(me["mouth_margin"])   # wider than both entering pipes
        outer.append(_cone_z(rm + wall, r2 + wall, zm, zm - Lc, xm, ym))
        inner.append(_cone_z(rm, r2, zm + 0.5, zm - Lc - 0.5, xm, ym))
        pts = _secondary_points(p, pair, secondary_h(p, pair))
        outer.append(_poly_pipe(pts, r2 + wall, bend)[0])
        inner.append(_poly_pipe(pts, r2, bend)[0])
    xc, yc, zc = _collector_centre(p)
    Lk = p.L(me["collector_length"])
    rcm = r2 + p.L(me["collector_mouth_over"]) + p.L(me["mouth_margin"])
    outer.append(_cone_z(rcm + wall, r3 + wall, zc, zc - Lk, xc, yc))
    inner.append(_cone_z(rcm, r3, zc + 0.5, zc - Lk - 0.5, xc, yc))
    Lo = p.L(me["outlet_length"])
    outer.append(_zc(r3 + wall, zc - Lk - Lo, zc - Lk, xc, yc))
    inner.append(_zc(r3, zc - Lk - Lo - 1.0, zc - Lk + 1.0, xc, yc))
    tf = p.L(m["flange_thickness"])
    outer.append(_zc(r3 + p.L(me["outlet_flange_over"]), zc - Lk - Lo, zc - Lk - Lo + tf, xc, yc))
    # per-port flanges on the head exhaust face
    yf = side_y(p, "exhaust")
    hz = p.L(m["flange_half_height"])
    zf = port_flange_z(p, "exhaust")
    holes = []
    for i in range(p.cylinders):
        bx = port_bolt_xz(p, "exhaust", i)
        x0 = bx[0][0] - p.L(m["exhaust_flange_bolt_margin"])
        x1 = bx[1][0] + p.L(m["exhaust_flange_bolt_margin"])
        outer.append(box(x0, x1, yf, yf + tf, zf - hz, zf + hz))
        holes.append(_cyl_y(r1, yf - 1.0, yf + tf + 1.0, p.cyl_x(i), zf))
        for (x, z) in bx:
            holes.append(_cyl_y(p.L(m["bolt_diameter"]) / 2 + p.L(m["bolt_hole_clearance"]), yf - 1.0, yf + tf + 1.0, x, z))
    body = outer[0].fuse(*outer[1:])
    return body.cut(*inner, *holes)


def _zc(r, z0, z1, x, y):
    from cad.primitives import cyl_z
    return cyl_z(r, min(z0, z1), max(z0, z1), x, y)


def _cone_z(r_top, r_bot, z_top, z_bot, x, y):
    from build123d import Cone, Pos
    from cad.primitives import BOTTOM
    return Pos(x, y, z_bot) * Cone(r_bot, r_top, z_top - z_bot, align=BOTTOM)


def manifold_bolt(p: Params, kind: str, x: float, z: float):
    """Manifold bolt along y: tapped into the head (zero-clearance thread fit), head on the flange outer face."""
    m = mc(p)
    s = p.side_sign(kind)
    yf = side_y(p, kind)
    t = p.L(m["flange_thickness"])
    r = p.L(m["bolt_diameter"]) / 2
    y_in = yf - s * p.L(m["bolt_thread_depth"])
    y_out = yf + s * t
    return _cyl_y(r, y_in, y_out, x, z) + _cyl_y(p.L(m["bolt_head_diameter"]) / 2, y_out,
                                                  y_out + s * p.L(m["bolt_head_height"]), x, z)
