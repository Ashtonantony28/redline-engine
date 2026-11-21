"""Simplified cylinder block: 4 bores at bore spacing, deck, main bulkheads with bearing
bores, open crankcase bays sized to the crank/big-end swing, optional siamesed water jacket.
Frame: crank axis along x, z up, deck at Params.deck_z. See spec/cad.yaml `block`.
"""
from __future__ import annotations

import math

from cad.conrod import big_end_bore_r
from cad.crankshaft import journal_spans
from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z


def block_x_extent(p: Params):
    mains = journal_spans(p)
    e = p.L(p.section("block")["end_wall"])
    return mains[0][0] - e, mains[-1][1] + e


def swing_radius(p: Params) -> float:
    """Largest radius swept by the crank (counterweights) or a big-end bolt boss corner."""
    c = p.section("conrod")
    boss_y = big_end_bore_r(p) + p.L(c["bolt_boss_width"])
    big_end = p.crank_r + math.hypot(boss_y, p.L(c["bolt_boss_height"]) / 2)
    cw = p.L(p.section("crankshaft")["counterweight_radius"])
    return max(big_end, cw)


def block(p: Params):
    b = p.section("block")
    x0, x1 = block_x_extent(p)
    hw = p.L(b["half_width"])
    z_bot = -p.L(b["skirt_below_axis"])
    deck = p.deck_z
    bore_bot = deck - p.L(b["bore_length"])
    # the skirt leaves the liner near BDC: relieve the bulkheads below the liner along the bore
    skirt_bdc = p.rod_len - p.crank_r - p.L(p.section("piston")["skirt_below_pin"])
    bore_cut_bot = min(bore_bot, skirt_bdc - p.L(b["piston_bdc_clearance"]))
    body = box(x0, x1, -hw, hw, z_bot, deck)

    # crankcase bays between the bulkheads, open at the bottom (oil pan not modelled)
    sw = swing_radius(p) + p.L(b["crank_swing_clearance"])
    mains = journal_spans(p)
    for i in range(p.cylinders):
        bx0, bx1 = mains[i][1], mains[i + 1][0]
        body -= cyl_x(sw, bx0, bx1)
        body -= box(bx0, bx1, -sw, sw, z_bot, 0.0)
        body -= box(bx0, bx1, -p.bore_r, p.bore_r, 0.0, bore_bot)
    # main bearing bores through the bulkheads
    main_r = (p.L(p.section("crankshaft")["main_diameter"]) + p.L(b["main_bore_over_journal"])) / 2
    body -= cyl_x(main_r, x0, x1)
    # bores
    for i in range(p.cylinders):
        body -= cyl_z(p.bore_r, bore_cut_bot, deck, p.cyl_x(i))
    if b["water_jacket"]:
        body -= water_jacket(p)
    return body


def water_jacket(p: Params):
    """Open-deck siamesed jacket: band around the outside of the joined liners."""
    b = p.section("block")
    deck = p.deck_z
    z0 = deck - p.L(b["water_jacket_depth"])
    r_in = p.bore_r + p.L(b["liner_thickness"])
    r_out = r_in + p.L(b["water_jacket_width"])
    outer = inner = None
    for i in range(p.cylinders):
        o = cyl_z(r_out, z0, deck, p.cyl_x(i))
        n = cyl_z(r_in, z0, deck, p.cyl_x(i))
        outer = o if outer is None else outer + o
        inner = n if inner is None else inner + n
    return outer - inner
