"""Connecting rod: big end with bolt bosses, tapered I-beam shank, small end.

Local frame: big-end bore axis along x through the origin, small-end bore axis along x at
z = conrod_length (exact centre distance from spec/engine.yaml). The rod, cap and bolts are
one solid (steel); bearing shells are not modelled. See spec/cad.yaml `conrod`.
"""
from __future__ import annotations

from cad.params import Params
from cad.primitives import box, cyl_x, prism_x


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
