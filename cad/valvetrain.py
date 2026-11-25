"""Valvetrain CAD (direct-acting buckets, shim under bucket, DOHC 16 V): valves, seat inserts, guides, spring
seats, helical springs, retainers, split collets, shims, buckets, camshafts with lobes generated from the kin.cam
flat-follower profiles (the M3 design), cam caps and cam-cap bolts.

Valve-local frame: origin at the centre of the closed valve face (in its roof face), +z up the valve axis
(Params.valve_axis); the valve opens toward -z. Every station along the axis comes from `stack()`:

    valve face (0) -> seat insert top -> spring-seat floor -> spring (installed height, spec/valvetrain.yaml)
    -> retainer -> stem tip (valve_length) -> shim -> bucket crown -> lash -> cam base circle -> cam axis

Cam-local frame: cam axis along x through the origin, built at crank angle 0 with the VVT parked; the lobe of
cylinder k points at its bucket when the crank reaches that lobe's nose angle (kin.cam centre + firing phase).
The cam turns by -theta/2 about +x (same sense as the crank, ADR 0002; belt drive 2:1). All values MODELLED.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from build123d import Circle, Location, Plane, Pos, Rot, Spline, extrude, sweep

from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z, profile_face, revolve_rz, tangent_fillet_arc, tube_z

MASS_SOLVE_STEPS = 40     # bisection steps of the valve tulip-height solve
SPRING_POINTS_PER_TURN = 24
CAM_PROFILE_POINTS = 360


def vc(p: Params) -> dict:
    return p.section("valvetrain")


def _r(p, key, sec=None):
    return p.L((sec or vc(p))[key]) / 2


# --------------------------------------------------------------------------------------------- stack
@lru_cache(maxsize=None)
def stack(p: Params, kind: str) -> dict:
    """Stations along the valve axis (mm from the closed valve face) for `kind`."""
    c = vc(p)
    sp = p.vt["spring"]
    L = p.L(c["valve_length"])
    rt = c["retainer"]
    s = {"tip": L}
    s["retainer_top"] = L - p.L(c["tip_protrusion"])
    s["retainer_bottom"] = s["retainer_top"] - p.L(rt["thickness"])
    s["spigot_bottom"] = s["retainer_bottom"] - p.L(rt["spigot_height"])
    s["spring_top"] = s["retainer_bottom"]
    s["spring_bottom"] = s["spring_top"] - p.L(sp["installed_height"])
    s["floor"] = s["spring_bottom"] - p.L(c["spring_seat"]["thickness"])
    g = c["guide"]
    s["guide_bottom"] = s["floor"] - p.L(g["below_floor"])
    s["guide_top"] = s["floor"] + p.L(g["above_floor"])
    s["insert_top"] = p.L(c["insert"]["height"])
    s["shim_top"] = L + p.L(c["shim"]["thickness"])
    s["bucket_top"] = s["shim_top"] + p.L(c["bucket"]["crown"])
    s["bucket_bottom"] = s["bucket_top"] - p.L(c["bucket"]["height"])
    s["tower_top"] = s["bucket_top"] - p.L(c["bucket"]["exposed"])
    s["collet_top"] = s["retainer_top"]
    s["collet_bottom"] = s["retainer_top"] - p.L(c["collet"]["height"])
    s["cam_axis"] = s["bucket_top"] + p.L(p.vt["cam"][kind]["lash"]) + p.L(p.vt["cam"]["base_circle_radius"])
    return s


def valve_frame(p: Params, idx: int, kind: str, x_local: float) -> Location:
    o = p.valve_origin(idx, kind, x_local)
    return Pos(*o) * Rot(-p.side_sign(kind) * p.valve_angle_deg(kind), 0, 0)


def axis_point(p: Params, kind: str, s: float, x: float = 0.0) -> np.ndarray:
    """World point at station s on a valve axis of cylinder-local x (y, z do not depend on x)."""
    y = p.valve_y(kind)
    return np.array([x, y, p.roof_z(y)]) + s * p.valve_axis(kind)


def cam_axis_yz(p: Params, kind: str):
    q = axis_point(p, kind, stack(p, kind)["cam_axis"])
    return float(q[1]), float(q[2])


def stem_r(p: Params) -> float:
    return p.L(p.engine["valves"]["stem_diameter"]) / 2


# --------------------------------------------------------------------------------------------- valve
def _valve_solid(p: Params, kind: str, tulip: float):
    c = vc(p)
    r = p.L(p.engine["valves"][f"{kind}_head_diameter"]) / 2
    m, w = p.L(c["head_margin"]), p.L(c["face_width"])
    rs, L = stem_r(p), stack(p, kind)["tip"]
    g = c["collet_groove"]
    gz1 = L - p.L(g["from_tip"])
    gz0 = gz1 - p.L(g["width"])
    gd = p.L(g["depth"])
    seg = [("L", (r, 0.0)), ("L", (r, m)), ("L", (r - w, m + w)), tangent_fillet_arc((r - w, m + w), rs, tulip),
           ("L", (rs, gz0)), ("L", (rs - gd, gz0)), ("L", (rs - gd, gz1)), ("L", (rs, gz1)), ("L", (rs, L)),
           ("L", (0.0, L))]
    return revolve_rz((0.0, 0.0), seg)


@lru_cache(maxsize=None)
def valve(p: Params, kind: str):
    """Valve with the tulip height solved so its mass equals spec/valvetrain.yaml valve_assembly mass."""
    c = vc(p)
    target = p.vt["valve_assembly"][f"{kind}_valve_mass"] / p.density("valve") * p.k**3
    lo = p.L(c["head_margin"]) + p.L(c["face_width"]) + 0.5
    hi = stack(p, kind)["guide_bottom"]
    if _valve_solid(p, kind, hi).volume < target:
        return _valve_solid(p, kind, hi)
    for _ in range(MASS_SOLVE_STEPS):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _valve_solid(p, kind, mid).volume < target else (lo, mid)
    return _valve_solid(p, kind, 0.5 * (lo + hi))


def seat_insert(p: Params, kind: str):
    """Seat ring: 45 deg seat on the valve face line, throat bore above (port start), OD in the head counterbore."""
    c = vc(p)
    ins = c["insert"]
    r = p.L(p.engine["valves"][f"{kind}_head_diameter"]) / 2
    m = p.L(c["head_margin"])
    gap = p.L(ins["radial_gap"])
    rt = p.section("head")["throat_fraction"] * r
    R = r + p.L(ins["wall"])
    h = stack(p, kind)["insert_top"]
    # seat line r = r - (s - m): from (r + gap, m - gap) to the throat radius
    pts = [(r + gap, 0.0), (R, 0.0), (R, h), (rt, h), (rt, m + (r - rt)), (r + gap, m - gap)]
    seg = [("L", q) for q in pts[1:]]
    return revolve_rz(pts[0], seg)


def insert_od_r(p: Params, kind: str) -> float:
    return p.L(p.engine["valves"][f"{kind}_head_diameter"]) / 2 + p.L(vc(p)["insert"]["wall"])


def guide(p: Params, kind: str):
    g = vc(p)["guide"]
    s = stack(p, kind)
    return tube_z(_r(p, "outer_diameter", g), stem_r(p) + p.L(g["clearance"]), s["guide_bottom"], s["guide_top"])


def spring_seat(p: Params, kind: str = "intake"):
    c = vc(p)
    s = stack(p, kind)
    r_in = _r(p, "outer_diameter", c["guide"]) + p.L(c["spring_seat"]["radial_clearance"])
    return tube_z(_r(p, "outer_diameter", c["spring_seat"]), r_in, s["floor"], s["spring_bottom"])


# --------------------------------------------------------------------------------------------- spring
@lru_cache(maxsize=None)
def spring_design(p: Params) -> dict:
    """Wire diameter and active coils from rate, solid height and mean diameter (k = G d^4 / (8 D^3 n_a),
    solid height = (n_a + 2 n_dead) d), all from spec/valvetrain.yaml + spec/cad.yaml spring."""
    sp, c = p.vt["spring"], vc(p)["spring"]
    k = float(sp["rate"])
    D = float(c["mean_diameter"])
    G = float(p.material("valve_spring")["G"])
    nd = float(c["dead_coils_each_end"])
    Hs = float(sp["solid_height"])

    def d_of(n):
        return (8.0 * k * D**3 * n / G) ** 0.25
    lo, hi = 0.5, 30.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if (mid + 2 * nd) * d_of(mid) < Hs else (lo, mid)
    n_a = 0.5 * (lo + hi)
    d = d_of(n_a)
    mass_active = float(p.material("valve_spring")["rho"]) * math.pi * d**2 / 4 * math.pi * D * n_a
    return {"wire_d_m": d, "active_coils": n_a, "total_coils": n_a + 2 * nd, "mean_d_m": D, "dead_each_end": nd,
            "rate_check_N_per_m": G * d**4 / (8 * D**3 * n_a), "solid_height_m": (n_a + 2 * nd) * d,
            "active_mass_kg": mass_active, "active_mass_spec_kg": sp["mass"]}


@lru_cache(maxsize=None)
def spring(p: Params, length: float):
    """Helical spring with closed, ground ends, `length` mm between its seats, bottom seat at z = 0."""
    sd = spring_design(p)
    d, D = p.L(sd["wire_d_m"]), p.L(sd["mean_d_m"])
    na, nd = sd["active_coils"], sd["dead_each_end"]
    nt = na + 2 * nd
    pa = (length - 2 * nd * d) / na
    t = np.linspace(0.0, nt, int(math.ceil(nt * SPRING_POINTS_PER_TURN)) + 1)
    z = np.where(t < nd, t * d, np.where(t < nd + na, nd * d + (t - nd) * pa, nd * d + na * pa + (t - nd - na) * d))
    ang = 2 * math.pi * t
    pts = [(D / 2 * math.cos(a), D / 2 * math.sin(a), zz) for a, zz in zip(ang, z)]
    path = Spline(*pts)
    tan = path.tangent_at(0.0)
    prof = Plane(origin=pts[0], z_dir=tan) * Circle(d / 2)
    coil = sweep(prof, path=path, is_frenet=True)
    return coil & cyl_z(D / 2 + d, 0.0, length)


def spring_located_length(p: Params, kind: str, lift: float) -> float:
    s = stack(p, kind)
    return s["spring_top"] - s["spring_bottom"] - lift


# --------------------------------------------------------------------------------------------- retainer, collets
def _cone_r(p: Params, kind: str, s: float) -> float:
    c = vc(p)
    st = stack(p, kind)
    return _r(p, "bore_top_diameter", c["retainer"]) - (st["retainer_top"] - s) * math.tan(
        math.radians(c["retainer"]["cone_half_angle_deg"]))


def retainer(p: Params, kind: str = "intake"):
    c = vc(p)
    rt = c["retainer"]
    s = stack(p, kind)
    R, Rs = _r(p, "outer_diameter", rt), _r(p, "spigot_diameter", rt)
    rb = _cone_r(p, kind, s["collet_bottom"])
    pts = [(rb, s["spigot_bottom"]), (Rs, s["spigot_bottom"]), (Rs, s["retainer_bottom"]), (R, s["retainer_bottom"]),
           (R, s["retainer_top"]), (_cone_r(p, kind, s["retainer_top"]), s["retainer_top"]), (rb, s["collet_bottom"])]
    return revolve_rz(pts[0], [("L", q) for q in pts[1:]])


def collet(p: Params, kind: str = "intake"):
    """One of the two collet halves (y > 0 side, split gap to its twin), bead in the stem groove."""
    c = vc(p)
    s = stack(p, kind)
    g = c["collet_groove"]
    rs = stem_r(p)
    gz1 = s["tip"] - p.L(g["from_tip"])
    gz0 = gz1 - p.L(g["width"])
    gd = p.L(g["depth"])
    z0, z1 = s["collet_bottom"], s["collet_top"]
    pts = [(rs, z0), (_cone_r(p, kind, z0), z0), (_cone_r(p, kind, z1), z1), (rs, z1), (rs, gz1), (rs - gd, gz1),
           (rs - gd, gz0), (rs, gz0)]
    ring = revolve_rz(pts[0], [("L", q) for q in pts[1:]])
    half_gap = p.L(c["collet"]["split_gap"]) / 2
    R = _r(p, "bore_top_diameter", c["retainer"]) + 1.0
    return ring & box(-R, R, half_gap, R, z0 - 1.0, z1 + 1.0)


def shim(p: Params, kind: str = "intake"):
    s = stack(p, kind)
    return cyl_z(_r(p, "diameter", vc(p)["shim"]), s["tip"], s["shim_top"])


def bucket_r(p: Params) -> float:
    return p.L(p.vt["tappet"]["diameter"]) / 2


def bucket(p: Params, kind: str = "intake"):
    """Inverted cup: crown on the shim, skirt in the bore; shim pocket under the crown."""
    b = vc(p)["bucket"]
    s = stack(p, kind)
    R = bucket_r(p)
    ri = R - p.L(b["wall"])
    cup = tube_z(R, ri, s["bucket_bottom"], s["shim_top"]) + cyl_z(R, s["shim_top"], s["bucket_top"])
    # shim locating boss under the crown, recessed for the shim
    rb = _r(p, "diameter", vc(p)["shim"]) + p.L(b["boss_wall"])
    cup += tube_z(rb, _r(p, "diameter", vc(p)["shim"]), s["tip"] - p.L(b["boss_depth"]), s["shim_top"])
    return cup


# --------------------------------------------------------------------------------------------- camshafts
@lru_cache(maxsize=None)
def cam_profiles(p: Params):
    from kin.cam import design_cams
    return design_cams(p.engine, p.vt)


def lobe_nose_deg(p: Params, kind: str, idx: int) -> float:
    """Crank angle [deg] of maximum lift of the `kind` lobe of cylinder idx (VVT parked)."""
    return math.degrees(cam_profiles(p)[kind].centre_crank) + p.firing_phase_deg(idx)


def lobe_points(p: Params, kind: str, nose_dir) -> list:
    """Closed lobe contour (y, z) in the cam frame with the nose pointing along nose_dir (unit (y, z))."""
    from kin.cam import flat_follower_contour
    prof = cam_profiles(p)[kind]
    phi = np.linspace(-math.pi, math.pi, CAM_PROFILE_POINTS, endpoint=False)
    d = prof.derivatives(phi)
    X, Z = flat_follower_contour(prof.base_radius, phi, d[0], d[1])
    X, Z = X * p.k, Z * p.k
    n = np.asarray(nose_dir, float)
    t = np.array([-n[1], n[0]])
    return [(float(xx * n[0] + zz * t[0]), float(xx * n[1] + zz * t[1])) for xx, zz in zip(X, Z)]


def cam_c(p: Params) -> dict:
    return vc(p)["cam"]


def journal_xs(p: Params) -> list:
    return [p.main_x(j) for j in range(p.n_mains)]


def camshaft_x_extent(p: Params):
    c = cam_c(p)
    from cad.front import belt_plane_x, cam_pulley_width
    x_front = belt_plane_x(p) - cam_pulley_width(p) / 2
    x_rear = journal_xs(p)[-1] + p.L(c["journal_width"]) / 2 + p.L(c["rear_overhang"])
    return x_front, x_rear


def lobe_nose_dir(p: Params, kind: str, idx: int):
    """Unit (y, z) of the lobe nose in the cam body frame at crank 0 (VVT parked)."""
    d = -p.valve_axis(kind)[1:]                    # from the cam axis toward the bucket
    b = math.radians(lobe_nose_deg(p, kind, idx) / 2.0)
    return (d[0] * math.cos(b) - d[1] * math.sin(b), d[0] * math.sin(b) + d[1] * math.cos(b))


@lru_cache(maxsize=None)
def camshaft(p: Params, kind: str):
    c = cam_c(p)
    rj, rc, rn, rb = (_r(p, k, c) for k in ("journal_diameter", "core_diameter", "nose_diameter", "bore_diameter"))
    x_front, x_rear = camshaft_x_extent(p)
    xs = journal_xs(p)
    jw = p.L(c["journal_width"])
    shaft = cyl_x(rc, xs[0] - jw / 2 - p.L(c["front_step"]), x_rear)
    shaft += cyl_x(rn, x_front, xs[0] - jw / 2 - p.L(c["front_step"]))
    for x in xs:
        shaft += cyl_x(rj, x - jw / 2, x + jw / 2)
    w = p.L(p.vt["cam"]["lobe_width"])
    for idx in range(p.cylinders):
        nd = lobe_nose_dir(p, kind, idx)
        for k2, _sign, xl, _y, _d in p.valve_layout():
            if k2 != kind:
                continue
            x = p.cyl_x(idx) + xl
            shaft += extrude(Plane.YZ.offset(x - w / 2) * _spline_face(lobe_points(p, kind, nd)), amount=w)
    shaft -= cyl_x(rb, x_front - 1.0, x_rear + 1.0)
    return shaft


def _spline_face(pts):
    from build123d import Face, Wire
    return Face(Wire([Spline(*[(u, v) for u, v in pts], periodic=True)]))


def cam_cap(p: Params):
    """Cam bearing cap in the cam frame (axis along x at the origin, split plane z = 0, cap above)."""
    cc = vc(p)["cam_cap"]
    c = cam_c(p)
    rj = _r(p, "journal_diameter", c) + p.L(c["bearing_clearance"])
    hw, h = p.L(cc["half_width"]), _r(p, "journal_diameter", c) + p.L(cc["height_over_journal"])
    w = p.L(c["journal_width"])
    cap = box(-w / 2, w / 2, -hw, hw, 0.0, h)
    from build123d import fillet as _fil
    try:
        top = [e for e in cap.edges() if abs(e.center().Z - h) < 1e-6 and abs(abs(e.center().Y) - hw) < 1e-6]
        cap = _fil(top, p.L(cc["edge_radius"]))
    except Exception:  # cosmetic: keep the square cap if the fillet fails
        pass
    cap -= cyl_x(rj, -w, w)
    bo = p.L(cc["bolt_offset"])
    for s in (-1, 1):
        cap -= cyl_z(_r(p, "bolt_diameter", cc) + p.L(cc["bolt_hole_clearance"]), -1.0, h + 1.0, 0.0, s * bo)
    return cap


def cam_cap_bolt(p: Params):
    """Flanged cap screw, head on the cap top (z = cap height), shank down into the head (cam frame)."""
    cc = vc(p)["cam_cap"]
    c = cam_c(p)
    h = _r(p, "journal_diameter", c) + p.L(cc["height_over_journal"])
    r = _r(p, "bolt_diameter", cc)
    return cyl_z(r, -p.L(cc["thread_depth"]), h) + cyl_z(_r(p, "bolt_head_diameter", cc), h,
                                                          h + p.L(cc["bolt_head_height"]))
