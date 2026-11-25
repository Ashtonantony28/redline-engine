"""Front end and flywheel: crank timing sprocket, torsional damper pulley (hub, rubber, inertia ring with poly-V
grooves), cam pulleys, belt tensioner and idler, the timing belt itself (one closed band on the real path: arcs on
every pulley joined by common tangents), and the flywheel with its ring gear and bolts.

The belt runs in the plane x = belt_plane_x; pulleys sit inside the loop (the path is the convex hull of the pulley
rims), the crank sprocket at the crank axis and the cam pulleys on the cam axes with the 2:1 tooth ratio of a
four-stroke (pitch radii from the tooth counts and belt pitch, spec/cad.yaml front). Teeth are not modelled: the
belt is a smooth band and the pulley rims sit at the belt inner face. Flywheel and damper ring are sized to the
inertias in spec/engine.yaml masses.flywheel_inertia and spec/dynamics.yaml torsional.damper.ring_inertia.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from build123d import Plane, Polygon, Pos, extrude

from cad.params import Params
from cad.primitives import cyl_x, revolve_rz, tube_x

SOLVE_STEPS = 40


def fc(p: Params) -> dict:
    return p.section("front")


def belt_plane_x(p: Params) -> float:
    """x of the belt centre plane: block front face minus the sprocket set-back and half the pulley width."""
    from cad.block import block_x_extent
    f = fc(p)
    return block_x_extent(p)[0] - p.L(f["setback"]) - cam_pulley_width(p) / 2


def cam_pulley_width(p: Params) -> float:
    f = fc(p)
    return p.L(f["belt_width"]) + 2 * p.L(f["flange_thickness"])


def pitch_r(p: Params, teeth: int) -> float:
    return teeth * p.L(fc(p)["belt_pitch"]) / (2 * math.pi)


def rim_r(p: Params, teeth: int) -> float:
    """Pulley rim radius = pitch radius minus the belt pitch-line offset (where the belt inner face sits)."""
    return pitch_r(p, teeth) - p.L(fc(p)["pitch_line_offset"])


def crank_teeth(p: Params) -> int:
    return int(fc(p)["crank_teeth"])


def cam_teeth(p: Params) -> int:
    return 2 * crank_teeth(p)          # four-stroke: cam turns once per two crank turns


# --------------------------------------------------------------------------------------------- pulleys
def _flanged_pulley_profile(p, r_rim, r_bore, w_rim, hub):
    """(r, x) profile of a flanged pulley about the x axis, x centred on 0: rim width w_rim between flanges."""
    f = fc(p)
    tf, hf = p.L(f["flange_thickness"]), p.L(f["flange_height"])
    rw = r_rim - p.L(f["rim_depth"])
    hw = w_rim / 2 + tf
    pts = [(r_bore, -hw - hub), (rw * 0.55, -hw - hub), (rw * 0.55, -hw), (r_rim + hf, -hw), (r_rim + hf, -w_rim / 2),
           (r_rim, -w_rim / 2), (r_rim, w_rim / 2), (r_rim + hf, w_rim / 2), (r_rim + hf, hw), (rw * 0.55, hw),
           (r_bore, hw)]
    return pts


def _revolve_x(pts_rx):
    """Solid of revolution about the x axis from (r, x) points."""
    from build123d import Rot
    solid = revolve_rz((pts_rx[0][0], pts_rx[0][1]), [("L", q) for q in pts_rx[1:]])
    return Rot(0, 90, 0) * solid      # local z (profile axis) -> global x


def crank_sprocket(p: Params):
    """Crank timing pulley on the snout (bore = snout), centred on the belt plane (crank frame)."""
    f = fc(p)
    rs = p.L(p.section("crankshaft")["snout_diameter"]) / 2
    pts = _flanged_pulley_profile(p, rim_r(p, crank_teeth(p)), rs, p.L(f["belt_width"]), 0.0)
    return Pos(belt_plane_x(p), 0, 0) * _revolve_x([(r, -x) for r, x in pts])


def cam_pulley(p: Params):
    """Cam pulley in the cam frame (cam axis = x through the origin), bore on the cam nose."""
    f = fc(p)
    rn = p.L(p.section("valvetrain")["cam"]["nose_diameter"]) / 2
    r_rim = rim_r(p, cam_teeth(p))
    pts = _flanged_pulley_profile(p, r_rim, rn, p.L(f["belt_width"]), 0.0)
    pul = Pos(belt_plane_x(p), 0, 0) * _revolve_x([(r, -x) for r, x in pts])
    # lightening holes in the web
    n = int(f["cam_pulley_holes"])
    rh = p.L(f["cam_pulley_hole_diameter"]) / 2
    rc = 0.5 * (rn + r_rim - p.L(f["rim_depth"]))
    for k in range(n):
        a = 2 * math.pi * k / n
        pul -= cyl_x(rh, belt_plane_x(p) - 50, belt_plane_x(p) + 50, rc * math.sin(a), rc * math.cos(a))
    return pul


def idler_pulley(p: Params, r_rim: float):
    """Plain idler/tensioner pulley (belt runs on its rim) with its mounting boss to the block front face; in a
    local frame centred on the pulley axis (axis along x), boss toward +x."""
    from cad.block import block_x_extent
    f = fc(p)
    w = cam_pulley_width(p)
    bx = belt_plane_x(p)
    rb = p.L(f["idler_boss_diameter"]) / 2
    body = tube_x(r_rim, r_rim - p.L(f["idler_rim_thickness"]), bx - w / 2, bx + w / 2)
    body += tube_x(r_rim - p.L(f["idler_rim_thickness"]) + 0.5, rb, bx - p.L(f["idler_web"]) / 2, bx + p.L(f["idler_web"]) / 2)
    body += cyl_x(rb, bx - w / 2, block_x_extent(p)[0])
    return body


# --------------------------------------------------------------------------------------------- belt path
def _tangent(c1, r1, c2, r2):
    """External tangent (both circles on the left of travel c1 -> c2, loop counter-clockwise in (y, z))."""
    d = np.asarray(c2) - np.asarray(c1)
    L = np.linalg.norm(d)
    u = d / L
    a = math.acos((r1 - r2) / L)
    ang = math.atan2(u[1], u[0])
    # outward normal on the right side of travel (circles inside a CCW loop lie on the left)
    th = ang - a
    n = np.array([math.cos(th), math.sin(th)])
    return np.asarray(c1) + r1 * n, np.asarray(c2) + r2 * n


@lru_cache(maxsize=None)
def belt_circles(p: Params):
    """[(name, centre (y, z), rim radius)] in counter-clockwise loop order around the (y, z) plane."""
    from cad.valvetrain import cam_axis_yz
    f = fc(p)
    yi, zi = cam_axis_yz(p, "intake")
    ye, ze = cam_axis_yz(p, "exhaust")
    rc, rk = rim_r(p, crank_teeth(p)), rim_r(p, cam_teeth(p))
    rt = p.L(f["tensioner_diameter"]) / 2
    circles = [("crank", (0.0, 0.0), rc), ("idler", tuple(p.L(v) for v in f["idler_yz"]), rt),
               ("exhaust", (ye, ze), rk), ("intake", (yi, zi), rk),
               ("tensioner", tuple(p.L(v) for v in f["tensioner_yz"]), rt)]
    return circles


def belt_path(p: Params):
    """Tangent points and arcs of the belt inner face: list of (circle, entry point, exit point)."""
    cs = belt_circles(p)
    n = len(cs)
    tang = []
    for i in range(n):
        _a, c1, r1 = cs[i]
        _b, c2, r2 = cs[(i + 1) % n]
        tang.append(_tangent(c1, r1, c2, r2))
    out = []
    for i in range(n):
        out.append((cs[i], tang[i - 1][1], tang[i][0]))
    return out, tang


def belt_length(p: Params) -> float:
    out, tang = belt_path(p)
    L = sum(np.linalg.norm(b - a) for a, b in tang)
    for (_n, c, r), pin, pout in out:
        a0 = math.atan2(pin[1] - c[1], pin[0] - c[0])
        a1 = math.atan2(pout[1] - c[1], pout[0] - c[0])
        L += r * ((a1 - a0) % (2 * math.pi))
    return L


def belt(p: Params):
    """Closed belt band of thickness belt_thickness outside the pulley rims, width belt_width (static)."""
    from build123d import Face, Line, ThreePointArc, Wire, offset
    f = fc(p)
    out, tang = belt_path(p)
    bx = belt_plane_x(p)
    pl = Plane.YZ.offset(bx)
    g = lambda q: pl.from_local_coords((float(q[0]), float(q[1])))
    edges = []
    n = len(out)
    for i in range(n):
        (_nm, c, r), pin, pout = out[i]
        a0 = math.atan2(pin[1] - c[1], pin[0] - c[0])
        a1 = a0 + ((math.atan2(pout[1] - c[1], pout[0] - c[0]) - a0) % (2 * math.pi))
        am = 0.5 * (a0 + a1)
        mid = (c[0] + r * math.cos(am), c[1] + r * math.sin(am))
        edges.append(ThreePointArc(g(pin), g(mid), g(pout)))
        edges.append(Line(g(tang[i][0]), g(tang[i][1])))
    inner = Face(Wire(edges))
    outer = offset(inner, p.L(f["belt_thickness"]))
    band = outer - inner
    w = p.L(f["belt_width"])
    return extrude(Pos(-w / 2, 0, 0) * band, amount=w, dir=(1, 0, 0))


# --------------------------------------------------------------------------------------------- damper
def damper_x0(p: Params) -> float:
    """x of the damper hub rear face = front face of the crank sprocket hub."""
    return belt_plane_x(p) - cam_pulley_width(p) / 2


def damper_hub(p: Params):
    d = fc(p)["damper"]
    rs = p.L(p.section("crankshaft")["snout_diameter"]) / 2
    x0 = damper_x0(p)
    L = p.L(d["hub_length"])
    r_in = p.L(d["ring_inner_diameter"]) / 2 - p.L(d["rubber_thickness"])
    web = p.L(d["web_thickness"])
    hub = tube_x(p.L(d["hub_diameter"]) / 2, rs, x0 - L, x0)
    xw0 = x0 - L / 2 - web / 2
    hub += tube_x(r_in, p.L(d["hub_diameter"]) / 2 - 0.5, xw0, xw0 + web)
    hub += tube_x(r_in, r_in - p.L(d["hub_rim_thickness"]), x0 - L / 2 - p.L(d["ring_width"]) / 2,
                  x0 - L / 2 + p.L(d["ring_width"]) / 2)
    return hub


def damper_rubber(p: Params):
    d = fc(p)["damper"]
    x0 = damper_x0(p)
    ri = p.L(d["ring_inner_diameter"]) / 2
    xm = x0 - p.L(d["hub_length"]) / 2
    w = p.L(d["ring_width"])
    return tube_x(ri, ri - p.L(d["rubber_thickness"]), xm - w / 2, xm + w / 2)


def _damper_ring(p: Params, r_out: float):
    d = fc(p)["damper"]
    x0 = damper_x0(p)
    ri = p.L(d["ring_inner_diameter"]) / 2
    xm = x0 - p.L(d["hub_length"]) / 2
    w = p.L(d["ring_width"])
    ng = int(d["poly_v_ribs"])
    depth = p.L(d["poly_v_depth"])
    pitch = p.L(d["poly_v_pitch"])
    # (r, x) profile: poly-V grooves across the outer face, centred on the ring
    pts = [(ri, xm - w / 2), (r_out, xm - w / 2)]
    xg0 = xm - ng * pitch / 2
    for k in range(ng):
        xa = xg0 + k * pitch
        pts += [(r_out, xa), (r_out - depth, xa + pitch / 2)]
    pts += [(r_out, xg0 + ng * pitch), (r_out, xm + w / 2), (ri, xm + w / 2)]
    return _revolve_x(pts)


@lru_cache(maxsize=None)
def damper_ring(p: Params):
    """Inertia ring, outer radius solved for spec/dynamics.yaml torsional.damper.ring_inertia."""
    from cad.massprops import mass_props
    d = fc(p)["damper"]
    target = float(p.dyn["torsional"]["damper"]["ring_inertia"])
    ri = p.L(d["ring_inner_diameter"]) / 2
    lo, hi = ri + p.L(d["poly_v_depth"]) + 1.0, 3 * ri

    def inertia(ro):
        mp = mass_props(_damper_ring(p, ro), p.density("damper_ring"), p.k)
        return float(mp.inertia_about((0, 0, 0))[0, 0])
    for _ in range(SOLVE_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if inertia(mid) < target else (lo, mid)
    return _damper_ring(p, 0.5 * (lo + hi))


# --------------------------------------------------------------------------------------------- flywheel
def flange_x(p: Params):
    """(x0, x1) of the crank rear flange (cad/crankshaft.py)."""
    from cad.crankshaft import journal_spans
    c = p.section("crankshaft")
    stub_end = journal_spans(p)[-1][1] + p.L(p.section("block")["end_wall"]) + p.L(c["flange_gap"])
    return stub_end, stub_end + p.L(c["flange_length"])


def _flywheel(p: Params, rim_width: float):
    fw = fc(p)["flywheel"]
    x0 = flange_x(p)[1]
    r_pilot = p.L(fw["pilot_diameter"]) / 2
    r_hub = p.L(fw["hub_diameter"]) / 2
    r_rim = p.L(fw["rim_inner_diameter"]) / 2
    r_out = p.L(fw["outer_diameter"]) / 2
    t_hub, t_web = p.L(fw["hub_thickness"]), p.L(fw["web_thickness"])
    step = p.L(fw["ring_gear_seat"])
    ring_w = p.L(fw["ring_gear_width"])
    pts = [(r_pilot, x0), (r_hub, x0), (r_hub, x0 + t_hub - t_web), (r_rim, x0 + t_hub - t_web), (r_rim, x0),
           (r_out - step, x0), (r_out - step, x0 + ring_w), (r_out, x0 + ring_w), (r_out, x0 + rim_width),
           (r_pilot, x0 + rim_width)]
    # clutch face side (+x) is flat; the -x side carries the hub, web pocket and ring-gear seat
    fwh = _revolve_x(pts)
    n = int(fw["bolts"])
    rb = p.L(fw["bolt_circle_diameter"]) / 2
    rh = p.L(fw["bolt_diameter"]) / 2 + p.L(fw["bolt_hole_clearance"])
    for k in range(n):
        a = 2 * math.pi * (k + 0.5) / n
        fwh -= cyl_x(rh, x0 - 1, x0 + rim_width + 1, rb * math.sin(a), rb * math.cos(a))
    return fwh


def ring_gear(p: Params):
    """Starter ring gear (trapezoid teeth), shrunk on the flywheel seat."""
    fw = fc(p)["flywheel"]
    x0 = flange_x(p)[1]
    r_in = p.L(fw["outer_diameter"]) / 2 - p.L(fw["ring_gear_seat"])
    n = int(fw["ring_gear_teeth"])
    r_root = p.L(fw["outer_diameter"]) / 2 + p.L(fw["ring_gear_root_over_rim"])
    r_tip = r_root + p.L(fw["ring_gear_tooth_height"])
    pts = []
    for k in range(n):
        a0 = 2 * math.pi * k / n
        da = 2 * math.pi / n
        for frac, r in ((0.0, r_root), (0.3, r_tip), (0.55, r_tip), (0.85, r_root)):
            a = a0 + frac * da
            pts.append((r * math.cos(a), r * math.sin(a)))
    outer = Plane.YZ.offset(x0) * Polygon(*pts, align=None)
    gear = extrude(outer, amount=p.L(fw["ring_gear_width"]))
    gear -= cyl_x(r_in, x0 - 1, x0 + p.L(fw["ring_gear_width"]) + 1)
    return gear


def flywheel_bolt(p: Params, k: int = 0, x_face: float | None = None):
    """Flywheel bolt k of the bolt circle: tapped through the crank flange (zero-clearance fit = thread), passing
    the flywheel hub in a clearance hole, head on the flat clutch-side face (x_face; the solved hub face by
    default)."""
    fw = fc(p)["flywheel"]
    fx0, fx1 = flange_x(p)
    rb = p.L(fw["bolt_circle_diameter"]) / 2
    r = p.L(fw["bolt_diameter"]) / 2
    a = 2 * math.pi * (k + 0.5) / int(fw["bolts"])
    y, z = rb * math.sin(a), rb * math.cos(a)
    xf = fx1 + rim_width(p) if x_face is None else x_face
    return cyl_x(r, fx0, xf, y, z) + cyl_x(p.L(fw["bolt_head_diameter"]) / 2, xf, xf + p.L(fw["bolt_head_height"]), y, z)


def bolt_yz(p: Params):
    fw = fc(p)["flywheel"]
    rb = p.L(fw["bolt_circle_diameter"]) / 2
    n = int(fw["bolts"])
    return [(rb * math.sin(2 * math.pi * (k + 0.5) / n), rb * math.cos(2 * math.pi * (k + 0.5) / n)) for k in range(n)]


@lru_cache(maxsize=None)
def rim_width(p: Params) -> float:
    """Flywheel width solved so flywheel + ring gear + bolts reach spec/engine.yaml masses.flywheel_inertia
    (inertia about the crank axis does not depend on where along x the bolt heads sit)."""
    from cad.massprops import mass_props
    fw = fc(p)["flywheel"]
    target = float(p.engine["masses"]["flywheel_inertia"])
    rest = mass_props(ring_gear(p), p.density("ring_gear"), p.k).inertia_about((0, 0, 0))[0, 0]
    for k in range(int(fw["bolts"])):
        bolt = mass_props(flywheel_bolt(p, k, flange_x(p)[1]), p.density("flywheel_bolt"), p.k)
        rest += bolt.inertia_about((0, 0, 0))[0, 0]
    lo, hi = p.L(fw["ring_gear_width"]) + 1.0, p.L(fw["max_rim_width"])

    def inertia(w):
        return float(mass_props(_flywheel(p, w), p.density("flywheel"), p.k).inertia_about((0, 0, 0))[0, 0]) + rest
    for _ in range(SOLVE_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if inertia(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


@lru_cache(maxsize=None)
def flywheel(p: Params):
    return _flywheel(p, rim_width(p))
