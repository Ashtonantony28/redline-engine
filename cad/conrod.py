"""Connecting rod: big end with bolt bosses, tapered I-beam shank, small end.

Local frame: big-end bore axis along x through the origin, small-end bore axis along x at
z = conrod_length (exact centre distance from spec/engine.yaml). The rod is split at z = 0
(the big-end axis) into the rod body and a separate cap held by two bolts on y = +-bolt_y
(`conrod_bolt` in spec/cad.yaml): the bolt shank is tapped into the rod body (modelled as a
zero-clearance cylindrical fit = thread engagement, allow-listed in spec/cad_fits.yaml), passes
the cap in a clearance hole, and its head seats on a spot face cut into the cap. The big-end
bearing shells are separate parts (cad/bearings.py). `conrod(p)` returns the uncut one-piece
envelope (rod body + cap, used by fea/); `conrod_parts(p)` returns the split parts.
"""
from __future__ import annotations

from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z, prism_x

BIG = 1.0e3  # half-space extent (mm), larger than any part


def big_end_bore_r(p: Params) -> float:
    c = p.section("conrod")
    return (p.L(p.section("crankshaft")["pin_diameter"]) + p.L(c["big_end_bore_over_pin"])) / 2


def small_end_bore_r(p: Params) -> float:
    return p.L(p.section("wrist_pin")["outer_diameter"]) / 2



def conrod(p: Params):
    c = p.section("conrod")
    L = p.rod_len
    be_r = big_end_bore_r(p)
    be_R = be_r + p.L(c["big_end_wall"])
    be_w = p.L(c["big_end_width"])
    se_r = small_end_bore_r(p)
    se_R = se_r + p.L(c["small_end_wall"])
    se_w = p.L(c["small_end_width"])

    big = cyl_x(be_R, -be_w / 2, be_w / 2)
    boss_y = be_r + p.L(c["bolt_boss_width"])
    boss_h = p.L(c["bolt_boss_height"])
    big += box(-be_w / 2, be_w / 2, -boss_y, boss_y, -boss_h / 2, boss_h / 2)
    small = cyl_x(se_R, -se_w / 2, se_w / 2, z=L)

    # shank: tapered outline in (y, z), I-beam pockets on both x faces
    hb = p.L(c["shank_width_big_end"]) / 2
    hs = p.L(c["shank_width_small_end"]) / 2
    fw = p.L(c["flange_width"]) / 2
    ft = p.L(c["flange_thickness"])
    tw = p.L(c["web_thickness"]) / 2
    z_b, z_s = 0.0, L  # run the outline centre-to-centre; the eyes are unioned over it
    outline = [(-hb, z_b), (hb, z_b), (hs, z_s), (-hs, z_s)]
    shank = prism_x(outline, -fw, fw)
    zp0, zp1 = be_R, L - se_R  # pockets only between the eyes
    taper = (hs - hb) / (z_s - z_b)
    hp0, hp1 = hb + taper * zp0 - ft, hb + taper * zp1 - ft
    pocket = [(-hp0, zp0), (hp0, zp0), (hp1, zp1), (-hp1, zp1)]
    for x0, x1 in ((tw, fw), (-fw, -tw)):
        shank -= prism_x(pocket, x0, x1)

    rod = big + small + shank
    rod -= cyl_x(be_r, -be_w, be_w)
    rod -= cyl_x(se_r, -se_w, se_w, z=L)
    return rod


def bolt_sec(p: Params) -> dict:
    return p.section("conrod_bolt")


def bolt_y(p: Params) -> float:
    """Bolt axis: midway between the big-end bore and the outer face of the bolt boss."""
    c = p.section("conrod")
    return big_end_bore_r(p) + p.L(c["bolt_boss_width"]) / 2


def cap_face_z(p: Params) -> float:
    """Underside of the cap bolt bosses (the bolt-head spot face is cut up from here)."""
    return -p.L(p.section("conrod")["bolt_boss_height"]) / 2


def _bolt_cyl(p: Params, r: float, z0: float, z1: float, s: int):
    return cyl_z(r, z0, z1, 0.0, s * bolt_y(p))


def conrod_bolt(p: Params, side: int = 1):
    """Rod bolt in rod-local coordinates: head in the cap spot face, shank tapped into the rod body."""
    b = bolt_sec(p)
    r = p.L(b["diameter"]) / 2
    seat = cap_face_z(p) + p.L(b["spot_face_depth"])
    head = _bolt_cyl(p, p.L(b["head_diameter"]) / 2, seat - p.L(b["head_height"]), seat, side)
    return head + _bolt_cyl(p, r, seat, p.L(b["thread_depth"]), side)


def conrod_parts(p: Params) -> dict:
    """Split connecting rod: {'conrod': body, 'conrod_cap': cap, 'conrod_bolt_a/b': bolts} (rod-local)."""
    b = bolt_sec(p)
    rod = conrod(p)
    body = rod & box(-BIG, BIG, -BIG, BIG, 0.0, BIG)
    cap = rod & box(-BIG, BIG, -BIG, BIG, -BIG, 0.0)
    r = p.L(b["diameter"]) / 2
    seat = cap_face_z(p) + p.L(b["spot_face_depth"])
    for s in (-1, 1):
        body -= _bolt_cyl(p, r, -1.0, p.L(b["thread_depth"]) + p.L(b["thread_bottom_gap"]), s)
        cap -= _bolt_cyl(p, r + p.L(b["hole_clearance"]), -BIG, 1.0, s)
        cap -= _bolt_cyl(p, p.L(b["head_diameter"]) / 2 + p.L(b["spot_face_clearance"]), -BIG, seat, s)
    return {"conrod": body, "conrod_cap": cap, "conrod_bolt_a": conrod_bolt(p, -1), "conrod_bolt_b": conrod_bolt(p, 1)}
