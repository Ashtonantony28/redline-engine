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
# 4-2-1 header routed like a production long-tube header: each primary leaves its port along +y, bends down and drops
# in a vertical plane parallel to the block side, then swings diagonally in x into a vertical 2-into-1 merge cone
# beside the block skirt at x = 0 (pairs 1+4 and 2+3 are 360 deg apart in the firing order 1-3-4-2 and symmetric
# about the block middle). Pair 1+4 drops in the inboard plane, pair 2+3 in an outboard plane (longer out leg), so
# no primaries cross. Each secondary drops from its cone, bends rearward and slopes down to the collector, whose axis
# runs along +x (towards the flywheel end) at about the sump floor; the outlet flange faces rearward. Centre lines are
# straight legs joined by tangent torus bends (FilletPolyline: analytic surfaces, robust booleans); every bend radius
# is at least min_bend_to_od x the tube outer diameter. The free leg of every duct is solved so the centre-line
# length equals the spec/gasdyn.yaml primary / secondary length.
def _poly_pipe(pts, r, bend):
    from build123d import FilletPolyline
    path = FilletPolyline(*[Vector(*q) for q in pts], radius=bend)
    t0 = (Vector(*pts[1]) - Vector(*pts[0])).normalized()
    return sweep(Plane(origin=Vector(*pts[0]), z_dir=t0) * Circle(r), path=path), path.length


def _corner_tangents(pts, bend):
    """Tangent length R tan(turn/2) at every interior corner of a polyline (0 at the ends)."""
    t = [0.0]
    for a, b, c in zip(pts, pts[1:], pts[2:]):
        u = (Vector(*b) - Vector(*a)).normalized()
        v = (Vector(*c) - Vector(*b)).normalized()
        turn = math.acos(max(-1.0, min(1.0, u.dot(v))))
        t.append(bend * math.tan(turn / 2.0))
    return t + [0.0]


def poly_length(pts, bend):
    """Exact centre-line length of a polyline whose corners are replaced by tangent arcs of radius bend, or None when
    a leg is shorter than the tangents of its two corners (the bends would overlap: not routable)."""
    if any((Vector(*b) - Vector(*a)).length < 1e-6 for a, b in zip(pts, pts[1:])):
        return None
    t = _corner_tangents(pts, bend)
    L = 0.0
    for k, (a, b) in enumerate(zip(pts, pts[1:])):
        s = (Vector(*b) - Vector(*a)).length
        if s < t[k] + t[k + 1] - 1e-9:
            return None
        L += s
    for k in range(1, len(pts) - 1):
        turn = 2.0 * math.atan(t[k] / bend) if bend > 0 else 0.0
        L -= 2.0 * t[k] - bend * turn
    return L


def _solve_leg(f, lo, hi, target, what):
    """Bisect the scalar leg parameter q in [lo, hi] for poly_length(f(q)) = target. The routable q form one interval;
    it is found on a fine scan first, then the length (monotonic in q) is bisected inside it."""
    n = 400
    qs = [lo + (hi - lo) * k / n for k in range(n + 1)]
    ok = [(q, f(q)) for q in qs]
    ok = [(q, L) for q, L in ok if L is not None]
    if not ok:
        raise ValueError(f"{what}: no routable layout (adjust spec/cad.yaml manifolds.exhaust)")
    (qa, La), (qb, Lb) = ok[0], ok[-1]
    if not min(La, Lb) <= target <= max(La, Lb):
        raise ValueError(f"{what}: spec length {target:.1f} mm outside the routable range {min(La, Lb):.1f}.."
                         f"{max(La, Lb):.1f} mm (adjust spec/cad.yaml manifolds.exhaust)")
    rising = Lb > La
    for _ in range(SOLVE_STEPS):
        mid = 0.5 * (qa + qb)
        Lm = f(mid)
        if Lm is None:
            raise ValueError(f"{what}: routable range is not an interval")
        if (Lm < target) == rising:
            qa = mid
        else:
            qb = mid
    return 0.5 * (qa + qb)


def _ex(p: Params) -> dict:
    return mc(p)["exhaust"]


def exhaust_radii(p: Params) -> dict:
    """Bore radii of the primaries, secondaries and outlet (spec/gasdyn.yaml diameters)."""
    g = p.gas["exhaust"]
    return {"primary": port_exit_r(p, "exhaust"), "secondary": p.L(g["secondary"]["diameter"]) / 2,
            "outlet": p.L(g["tailpipe_front"]["diameter"]) / 2}


@lru_cache(maxsize=None)
def bend_radii(p: Params) -> dict:
    """Centre-line bend radius of the primaries and secondaries; each must be >= min_bend_to_od x the tube OD."""
    me, w, r = _ex(p), p.L(mc(p)["wall"]), exhaust_radii(p)
    out = {k: p.L(me[f"{k}_bend_radius"]) for k in ("primary", "secondary")}
    for k, R in out.items():
        od = 2.0 * (r[k] + w)
        if R < float(me["min_bend_to_od"]) * od - 1e-9:
            raise ValueError(f"{k} bend radius {R:.1f} mm < {me['min_bend_to_od']} x OD {od:.1f} mm")
    return out


def _flange_out_y(p: Params) -> float:
    return side_y(p, "exhaust") + p.L(mc(p)["flange_thickness"])


def _pair_side(p: Params, idx: int):
    """(pair, side): pair 0 = cylinders 1+4 (inboard plane), pair 1 = 2+3 (outboard); side = sign of the cylinder x."""
    return (0 if idx in (0, p.cylinders - 1) else 1), (-1.0 if p.cyl_x(idx) < 0 else 1.0)


def primary_out(p: Params, pair: int) -> float:
    """Straight leg of a pair's primaries out of the head flange (the pair's drop plane is at flange face + this)."""
    me = _ex(p)
    return p.L(me["primary_out"]) + pair * p.L(me["outer_pair_extra_out"])


def merge_centre(p: Params, pair: int):
    """Merge-cone mouth centre (x, y, z): x = 0, y = the pair's drop plane, merge_drop (+ the outer pair's extra)
    below the port centre."""
    me = _ex(p)
    z = port_flange_z(p, "exhaust") - p.L(me["merge_drop"]) - pair * p.L(me["outer_pair_drop_offset"])
    return 0.0, _flange_out_y(p) + primary_out(p, pair), z


def _primary_points(p: Params, idx: int, h: float):
    """Flange face -> out along +y -> down by h (solved) -> diagonal in x -> vertical entry into the merge mouth."""
    me = _ex(p)
    yf, zf = _flange_out_y(p), port_flange_z(p, "exhaust")
    pair, side = _pair_side(p, idx)
    xm, ym, zm = merge_centre(p, pair)
    x = p.cyl_x(idx)
    xe = xm + side * p.L(me["merge_offset"])
    z3 = zm + p.L(me["entry_length"])
    return [(x, yf, zf), (x, ym, zf), (x, ym, zf - h), (xe, ym, z3), (xe, ym, zm)]


@lru_cache(maxsize=None)
def primary_drop(p: Params, idx: int) -> float:
    """Vertical drop h of primary idx solved for the spec primary length (length grows with h)."""
    R = bend_radii(p)["primary"]
    zf = port_flange_z(p, "exhaust")
    z3 = merge_centre(p, _pair_side(p, idx)[0])[2] + p.L(_ex(p)["entry_length"])
    return _solve_leg(lambda h: poly_length(_primary_points(p, idx, h), R), 0.0, zf - z3,
                      p.L(p.gas["exhaust"]["primary"]["length"]), f"primary {idx + 1}")


def primary_points(p: Params, idx: int):
    return _primary_points(p, idx, primary_drop(p, idx))


def primary_length(p: Params, idx: int) -> float:
    """Measured CAD centre-line length (FilletPolyline wire), head flange face to merge-cone mouth."""
    return _poly_pipe(primary_points(p, idx), exhaust_radii(p)["primary"], bend_radii(p)["primary"])[1]


def collector_centre(p: Params):
    """Collector mouth centre: collector_from_rear_face behind the block rear face, collector_y outboard of the head
    flange face, collector_below_pan_floor below the sump floor; the collector axis runs along +x."""
    from cad.block import block_x_extent, pan_floor_z
    me = _ex(p)
    return (block_x_extent(p)[1] + p.L(me["collector_from_rear_face"]), _flange_out_y(p) + p.L(me["collector_y"]),
            pan_floor_z(p) - p.L(me["collector_below_pan_floor"]))


def _secondary_points(p: Params, pair: int, zb: float):
    """Merge-cone outlet -> vertical drop to zb (solved) -> rearward slope -> straight entry along +x into the
    collector mouth beside the other secondary (offset +- collector_mouth_over in y)."""
    me = _ex(p)
    xm, ym, zm = merge_centre(p, pair)
    xc, yc, zc = collector_centre(p)
    ye = yc + (1.0 if pair else -1.0) * p.L(me["collector_mouth_over"])
    z0 = zm - p.L(me["merge_cone_length"])
    return [(xm, ym, z0), (xm, ym, zb), (xc - p.L(me["collector_entry_length"]), ye, zc), (xc, ye, zc)]


@lru_cache(maxsize=None)
def secondary_bend_z(p: Params, pair: int) -> float:
    """Height zb of a secondary's rearward bend solved for the spec secondary length (length falls as zb rises)."""
    R = bend_radii(p)["secondary"]
    z0 = merge_centre(p, pair)[2] - p.L(_ex(p)["merge_cone_length"])
    zc = collector_centre(p)[2]
    return _solve_leg(lambda zb: poly_length(_secondary_points(p, pair, zb), R), zc, z0,
                      p.L(p.gas["exhaust"]["secondary"]["length"]), f"secondary {pair}")


def secondary_points(p: Params, pair: int):
    return _secondary_points(p, pair, secondary_bend_z(p, pair))


def secondary_length(p: Params, pair: int) -> float:
    """Measured CAD centre-line length, merge-cone outlet to collector mouth."""
    return _poly_pipe(secondary_points(p, pair), exhaust_radii(p)["secondary"], bend_radii(p)["secondary"])[1]


def _overlap(pts, d0: float, d1: float):
    """The same centre line with its end legs extended straight by d0 before the start and d1 past the end, so a duct
    sinks into the flange / cone it butts against and the union is one solid (a butt joint on a shared face is not
    fused reliably). Lengths are always measured on the unextended line."""
    a, b = Vector(*pts[0]), Vector(*pts[1])
    c, e = Vector(*pts[-2]), Vector(*pts[-1])
    s = a - (b - a).normalized() * d0
    t = e + (e - c).normalized() * d1
    return [tuple(s)] + list(pts[1:-1]) + [tuple(t)]


def outlet_flange_x(p: Params) -> float:
    """Rear face of the outlet flange (downpipe joint)."""
    me = _ex(p)
    return collector_centre(p)[0] + p.L(me["collector_length"]) + p.L(me["outlet_length"])


@lru_cache(maxsize=None)
def exhaust_header(p: Params):
    m = mc(p)
    me = _ex(p)
    wall = p.L(m["wall"])
    R = bend_radii(p)
    r = exhaust_radii(p)
    r1, r2, r3 = r["primary"], r["secondary"], r["outlet"]
    ov = p.L(me["joint_overlap"])
    if not 0.0 < p.L(me["mouth_margin"]) < wall:
        raise ValueError("mouth_margin must be in (0, wall): the entering pipe walls must sink into the cone wall")
    outer, inner = [], []
    for i in range(p.cylinders):
        pts = _overlap(primary_points(p, i), ov, ov)
        outer.append(_poly_pipe(pts, r1 + wall, R["primary"])[0])
        inner.append(_poly_pipe(pts, r1, R["primary"])[0])
    Lc = p.L(me["merge_cone_length"])
    for pair in (0, 1):
        xm, ym, zm = merge_centre(p, pair)
        rm = r1 + p.L(me["merge_offset"]) + p.L(me["mouth_margin"])   # clears both entering bores
        outer.append(_cone_z(rm + wall, r2 + wall, zm, zm - Lc, xm, ym))
        inner.append(_cone_z(rm, r2, zm + 0.5, zm - Lc - 0.5, xm, ym))
        pts = _overlap(secondary_points(p, pair), ov, ov)
        outer.append(_poly_pipe(pts, r2 + wall, R["secondary"])[0])
        inner.append(_poly_pipe(pts, r2, R["secondary"])[0])
    # collector cone and outlet along +x, rearward-facing outlet flange with its bolt holes
    xc, yc, zc = collector_centre(p)
    Lk = p.L(me["collector_length"])
    rcm = r2 + p.L(me["collector_mouth_over"]) + p.L(me["mouth_margin"])
    outer.append(_cone_x(rcm + wall, r3 + wall, xc, xc + Lk, yc, zc))
    inner.append(_cone_x(rcm, r3, xc - 0.5, xc + Lk + 0.5, yc, zc))
    xo = outlet_flange_x(p)
    outer.append(cyl_x(r3 + wall, xc + Lk - ov, xo, yc, zc))
    inner.append(cyl_x(r3, xc + Lk - 1.0, xo + 1.0, yc, zc))
    tf = p.L(m["flange_thickness"])
    over = p.L(me["outlet_flange_over"])
    outer.append(cyl_x(r3 + over, xo - tf, xo, yc, zc))
    holes = []
    rbc = r3 + 0.5 * (wall + over)                       # bolt circle midway between tube wall and flange rim
    nb = int(me["outlet_flange_bolts"])
    for k in range(nb):
        a = 2.0 * math.pi * (k + 0.5) / nb
        holes.append(cyl_x(p.L(m["bolt_diameter"]) / 2 + p.L(m["bolt_hole_clearance"]), xo - tf - 1.0, xo + 1.0,
                           yc + rbc * math.cos(a), zc + rbc * math.sin(a)))
    # welded stay pad on the 1+4 secondary (its bore is cut with the tube's), tapped for the stay foot bolt
    sg = stay_geom(p)
    outer.append(box(sg["xs"] - sg["pl"] / 2, sg["xs"] + sg["pl"] / 2, sg["yp"] - sg["pw"] / 2, sg["yp"] + sg["pw"] / 2,
                     sg["zaxis"], sg["zp"]))
    for (x, y, z) in stay_bolts(p)["pad"]:
        holes.append(_zcyl(p.L(m["bolt_diameter"]) / 2, z - p.L(me["stay"]["pad_thread_depth"]), z + 1.0, x, y))
    # per-port flanges on the head exhaust face
    yf = side_y(p, "exhaust")
    hz = p.L(m["flange_half_height"])
    zf = port_flange_z(p, "exhaust")
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


# Header stay (spec manifolds.exhaust.stay): the header hangs ~0.5 m from the head flanges, so a pressed-steel stay
# ties the 1+4 secondary to the block: a tab bolted to the external skirt rib of one bulkhead (solid bulkhead behind
# it for the thread), an arm in a y-z plane (the tab's front edge, clear of the bolt heads) dropping outboard past the
# pan-rail lip, and a foot bolted to a pad welded on top of the secondary.
def _secondary_axis_at_x(p: Params, pair: int, x: float):
    """Point of a secondary's centre line at a given x (on its rearward run)."""
    from build123d import FilletPolyline
    w = FilletPolyline(*[Vector(*q) for q in secondary_points(p, pair)], radius=bend_radii(p)["secondary"])
    lo, hi = 0.0, 1.0
    for _ in range(SOLVE_STEPS + 20):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if w.position_at(mid).X < x else (lo, mid)
    q = w.position_at(0.5 * (lo + hi))
    return q.X, q.Y, q.Z


@lru_cache(maxsize=None)
def stay_geom(p: Params) -> dict:
    from cad.block import bsec, rail_z
    st, w = _ex(p)["stay"], p.L(mc(p)["wall"])
    b = bsec(p)
    xs = p.main_x(int(st["bulkhead"]))
    yr = p.L(b["half_width"]) + p.L(b["rib_height"])           # rib face
    _x, yp, zaxis = _secondary_axis_at_x(p, 0, xs)
    zp = zaxis + exhaust_radii(p)["secondary"] + w + p.L(st["pad_height"])
    t, tw = p.L(st["thickness"]), p.L(st["tab_width"])
    z0, z1 = (p.L(v) for v in st["tab_z"])
    lip_top = rail_z(p) + p.L(b["pan_rail_height"])
    if z0 <= lip_top:
        raise ValueError("stay tab must sit above the pan-rail lip")
    return {"xs": xs, "yr": yr, "yp": yp, "zaxis": zaxis, "zp": zp, "t": t, "tw": tw, "z0": z0, "z1": z1,
            "xa0": xs - tw / 2, "pl": p.L(st["pad_length"]), "pw": p.L(st["pad_width"])}


def stay_bolts(p: Params) -> dict:
    """Placements (x, y, z) of the stay bolts: 'block' (origin on the rib face, thread along -y), 'pad' (origin on
    the pad top, thread along -z)."""
    g, st = stay_geom(p), _ex(p)["stay"]
    return {"block": [(g["xs"], g["yr"], p.L(z)) for z in st["bolt_z"]], "pad": [(g["xs"], g["yp"], g["zp"])]}


@lru_cache(maxsize=None)
def exhaust_stay(p: Params):
    m, st, g = mc(p), _ex(p)["stay"], stay_geom(p)
    t, yr, yp, zp = g["t"], g["yr"], g["yp"], g["zp"]
    xa0, xa1 = g["xa0"], g["xa0"] + g["tw"]
    tab = box(xa0, xa1, yr, yr + t, g["z0"], g["z1"])
    fl = p.L(st["foot_length"])
    foot = box(xa0, xa1, yp - fl / 2, yp + fl / 2, zp, zp + t)
    ky, kz = (p.L(v) for v in st["knee"])
    arm_pts = [(yr, g["z1"]), (yr, g["z0"]), (yr + ky, g["z0"] + kz), (yp - fl / 2, zp), (yp + fl / 2, zp),
               (yp + fl / 2, zp + p.L(st["arm_end_height"])), (yr + p.L(st["arm_top_width"]), g["z1"])]
    from build123d import Plane as _Pl, Polygon, extrude
    arm = extrude(_Pl.YZ.offset(xa0) * Polygon(*arm_pts, align=None), amount=t)
    body = tab.fuse(foot, arm).clean()
    hole = p.L(m["bolt_diameter"]) / 2 + p.L(m["bolt_hole_clearance"])
    tools = [_cyl_y(hole, yr - 1.0, yr + t + 1.0, x, z) for (x, _y, z) in stay_bolts(p)["block"]]
    tools += [_zcyl(hole, zp - 1.0, zp + t + 1.0, x, y) for (x, y, _z) in stay_bolts(p)["pad"]]
    return body.cut(*tools)


def stay_bolt(p: Params, kind: str):
    """Stay bolt in its local frame (origin on the clamped face of the threaded part): 'block' along -y, 'pad' -z."""
    m, st, g = mc(p), _ex(p)["stay"], stay_geom(p)
    r, t = p.L(m["bolt_diameter"]) / 2, g["t"]
    rh, hh = p.L(m["bolt_head_diameter"]) / 2, p.L(m["bolt_head_height"])
    if kind == "block":
        return _cyl_y(r, -p.L(st["thread_depth"]), t, 0.0, 0.0) + _cyl_y(rh, t, t + hh, 0.0, 0.0)
    return _zcyl(r, -p.L(st["pad_thread_depth"]), t, 0.0, 0.0) + _zcyl(rh, t, t + hh, 0.0, 0.0)


def _cone_z(r_top, r_bot, z_top, z_bot, x, y):
    from build123d import Cone, Pos
    from cad.primitives import BOTTOM
    return Pos(x, y, z_bot) * Cone(r_bot, r_top, z_top - z_bot, align=BOTTOM)


def _cone_x(r0, r1, x0, x1, y, z):
    """Cone along +x from radius r0 at x0 to r1 at x1."""
    from build123d import Cone, Pos, Rot
    from cad.primitives import BOTTOM
    return Pos(x0, y, z) * Rot(0, 90, 0) * Cone(r0, r1, x1 - x0, align=BOTTOM)


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
