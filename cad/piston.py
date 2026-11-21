"""Piston, wrist pin, rings and clips. Local frame: pin axis along x through the origin,
crown up (+z). See spec/cad.yaml `piston`, `wrist_pin`, `rings`, `clips`."""
from __future__ import annotations


from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z, tube_x, tube_z


def piston_radius(p: Params) -> float:
    return (p.bore - p.L(p.section("piston")["skirt_clearance"])) / 2


def groove_spans(p: Params):
    """[(z_top, z_bottom)] of each ring groove, crown at z = compression height."""
    c = p.section("piston")
    z = p.compression_height - p.L(c["top_land"])
    spans = []
    for i, h in enumerate(c["groove_heights"]):
        if i:
            z -= p.L(c["ring_land"])
        spans.append((z, z - p.L(h)))
        z -= p.L(h)
    return spans


def piston(p: Params):
    c = p.section("piston")
    R = piston_radius(p)
    top = p.compression_height
    bottom = -p.L(c["skirt_below_pin"])
    crown_under = top - p.L(c["crown_thickness"])
    belt_bottom = groove_spans(p)[-1][1] - p.L(c["belt_extra_below_rings"])
    r_belt = R - p.L(c["ring_belt_thickness"])
    r_skirt = R - p.L(c["skirt_wall"])
    flat_x = p.L(c["skirt_flat_half_width"])
    span = 2 * R  # generous extent for cutting boxes

    body = cyl_z(R, bottom, top)
    body -= cyl_z(r_skirt, bottom, belt_bottom)
    body -= cyl_z(r_belt, belt_bottom, crown_under)
    # relieve the skirt on the pin sides, then close it with side panels
    for s in (-1, 1):
        x0, x1 = sorted((s * flat_x, s * span))
        body -= box(x0, x1, -span, span, bottom, belt_bottom)
    panel_t = p.L(c["side_panel_thickness"])
    boss_r = p.L(c["boss_radius"])
    panels = None
    for s in (-1, 1):
        x0, x1 = sorted((s * flat_x, s * (flat_x - panel_t)))
        pan = box(x0, x1, -span, span, bottom, belt_bottom) & cyl_z(R, bottom, top)
        panels = pan if panels is None else panels + pan
    body += panels
    # pin bosses, tied up into the crown
    rib = p.L(c["boss_rib_thickness"])
    gap = p.L(p.section("conrod")["small_end_width"]) / 2 + p.L(c["boss_side_gap"])
    for s in (-1, 1):
        x0, x1 = sorted((s * gap, s * flat_x))
        boss = cyl_x(boss_r, x0, x1) + box(x0, x1, -rib / 2, rib / 2, 0.0, crown_under)
        body += boss & cyl_z(R, bottom, top)
    body -= cyl_x(p.L(p.section("wrist_pin")["outer_diameter"]) / 2, -span, span)
    # ring grooves
    r_groove = R - p.L(c["groove_depth"])
    for zt, zb in groove_spans(p):
        body -= tube_z(span, r_groove, zb, zt)
    for cp in clips(p):  # circlip grooves in the pin bore
        body -= cp
    if c["valve_reliefs"]:
        body -= valve_reliefs(p, top)
    return body


def valve_reliefs(p: Params, crown_z: float):
    """Shallow cylindrical pockets under each valve head (optional, cad.yaml)."""
    c = p.section("piston")
    depth = p.L(c["valve_relief_depth"])
    out = None
    for _kind, _sign, x, y, d in p.valve_layout():
        pocket = cyl_z((d + p.L(c["valve_relief_margin"])) / 2, crown_z - depth, crown_z, x, y)
        out = pocket if out is None else out + pocket
    return out & cyl_z(piston_radius(p), crown_z - depth, crown_z)


def pin_half_length(p: Params) -> float:
    """The pin ends just inside the relieved boss faces, retained by the clips."""
    flat_x = p.L(p.section("piston")["skirt_flat_half_width"])
    return flat_x - p.L(p.section("wrist_pin")["clip_inset"])


def wrist_pin(p: Params):
    w = p.section("wrist_pin")
    half = pin_half_length(p)
    return tube_x(p.L(w["outer_diameter"]) / 2, p.L(w["inner_diameter"]) / 2, -half, half)


def rings(p: Params):
    """Compression and oil rings seated in their grooves (outer face at the bore)."""
    c, r = p.section("piston"), p.section("rings")
    r_out = p.bore_r
    r_in = piston_radius(p) - p.L(c["groove_depth"]) + p.L(r["back_clearance"])
    side = p.L(r["side_clearance"])
    out = []
    for zt, zb in groove_spans(p):
        out.append(tube_z(r_out, r_in, zb + side / 2, zt - side / 2))
    return out


def clips(p: Params):
    """Wire circlips retaining the pin, one each side, modelled as square-section annuli."""
    w, cl = p.section("wrist_pin"), p.section("clips")
    half = pin_half_length(p)
    t = p.L(cl["wire_width"])
    r_out = (p.L(w["outer_diameter"]) + p.L(cl["outer_diameter_over_pin"])) / 2
    r_in = r_out - t
    return [tube_x(r_out, r_in, s * half + (0 if s > 0 else -t), s * half + (t if s > 0 else 0))
            for s in (-1, 1)]


def piston_parts(p: Params) -> dict:
    """{name: (shape, material key)} for one piston assembly, local frame."""
    parts = {"piston": (piston(p), "piston"), "wrist_pin": (wrist_pin(p), "wrist_pin")}
    for i, rg in enumerate(rings(p)):
        parts[f"ring_{i + 1}"] = (rg, "rings")
    for i, cp in enumerate(clips(p)):
        parts[f"clip_{i + 1}"] = (cp, "clips")
    return parts


