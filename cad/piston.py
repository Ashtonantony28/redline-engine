"""Piston, wrist pin, rings and clips. Local frame: pin axis along x through the origin,
crown up (+z). See spec/cad.yaml `piston`, `wrist_pin`, `rings`, `clips`."""
from __future__ import annotations


from build123d import Axis, Circle, Plane, Polygon, Pos, Rectangle, extrude, fillet, revolve

from cad.params import Params
from cad.primitives import box, cyl_x, cyl_z, tube_x, tube_z

TOL = 1.0e-6          # CAD units, vertex matching in 2D profiles
BIG = 1.0e3           # CAD units, extent of cutting boxes and planform strips
ON_CIRCLE_TOL = 1.0e-3


def fillet_vertices(sketch, radius: float, pred):
    """2D fillet of the sketch vertices matching pred (unchanged if none match)."""
    vs = [v for v in sketch.vertices() if pred(v)]
    return fillet(vs, radius) if vs else sketch


def piston_radius(p: Params) -> float:
    return (p.bore - p.L(p.section("piston")["skirt_clearance"])) / 2


def boss_gap(p: Params) -> float:
    """|x| of the pin-boss inner faces: half the conrod small-end width plus the side gap."""
    return p.L(p.section("conrod")["small_end_width"]) / 2 + p.L(p.section("piston")["boss_side_gap"])


def bridge_bottom(p: Params) -> float:
    """z of the crown-bridge underside: conrod small-end outer radius plus the bridge clearance."""
    from cad.conrod import small_end_bore_r
    return (small_end_bore_r(p) + p.L(p.section("conrod")["small_end_wall"])
            + p.L(p.section("piston")["bridge_clearance"]))


def boss_pillars(p: Params, R: float, top: float, cu: float, bb: float, rb: float):
    """Pin-boss pillars, belt ties and crown bridge (M6 redesign, ADR 0006, fea/design_history.json; shared by the
    CAD piston and the FEA reconstruction so both carry the same geometry).

    Each boss is a pillar: the boss cylinder carried straight up into the crown underside with vertical
    flanks tangent to it at the pin-axis level (pillar width = boss diameter), flaring into the crown with
    a tangent haunch of pillar_crown_radius, so the gas load goes from the crown directly into the pin-bore
    top. boss_belt_tie runs each pillar on outward past the skirt flat into the ring belt above the belt
    bottom (crown edge in the pin plane tied to the boss). crown_bridge joins the two pillars under the
    crown as a vault whose underside is the conrod small-end envelope (cylinder about the pin axis of
    radius `bridge_bottom`), flattened at bridge_side_height towards the flanks (no knife edge where the
    vault would meet the flanks). Re-entrant corners carry pillar_blend_radius blends. Built as the
    intersection of a (y, z) prism and an (x, z) prism with 2D-filleted profiles (3D fillets on this
    boolean topology are not reliable), plus planform webs (tie flanks into the belt wall) and a revolved
    vault cutter (bridge)."""
    c = p.section("piston")
    br, f, fb = p.L(c["boss_radius"]), p.L(c["pillar_crown_radius"]), p.L(c["pillar_blend_radius"])
    x0, x1 = boss_gap(p), p.L(c["skirt_flat_half_width"])
    tie, bridge, zb = bool(c["boss_belt_tie"]), bool(c["crown_bridge"]), bridge_bottom(p)
    xe = R + 1.0 if tie else x1
    lo, cap = -br, f
    yz = Circle(br) + Pos(0, (cu + cap) / 2) * Rectangle(2 * br, cu + cap) \
        + Pos(0, cu + cap / 2) * Rectangle(2 * br + 4 * f, cap)
    yz = fillet_vertices(yz, f, lambda v: abs(abs(v.X) - br) < TOL and abs(v.Y - cu) < TOL)
    right = [(x1, lo)]
    right += [(x1, bb), (xe, bb), (xe, cu + cap)] if tie else [(x1, cu), (x1 + 2 * fb, cu), (x1 + 2 * fb, cu + cap)]
    pts = [(x0, lo)] + right + [(-x, z) for x, z in reversed(right)] + [(-x0, lo)]
    zf = p.L(c["bridge_side_height"])
    if bridge:
        pts += [(-x0, zf), (x0, zf)]
    elif x0 - 2 * fb > 0.0:      # inner pillar faces flare into the crown with the blend radius
        pts += [(-x0, cu), (-x0 + 2 * fb, cu), (-x0 + 2 * fb, cu + cap), (x0 - 2 * fb, cu + cap), (x0 - 2 * fb, cu),
                (x0, cu)]
    else:
        raise ValueError("pillar_blend_radius too large for the boss gap without a crown bridge")
    xz = Polygon(*pts, align=None)
    corner_z = bb if tie else cu
    xz = fillet_vertices(xz, fb, lambda v: abs(abs(v.X) - x1) < TOL and abs(v.Y - corner_z) < TOL)
    if bridge:
        xz = fillet_vertices(xz, fb, lambda v: abs(abs(v.X) - x0) < TOL and abs(v.Y - zf) < TOL)
    else:
        xz = fillet_vertices(xz, fb, lambda v: abs(abs(v.X) - x0) < TOL and abs(v.Y - cu) < TOL)
    xm = xe + 3 * f
    big = 4 * R
    p_yz = extrude(Plane.YZ.offset(-xm) * yz, amount=2 * xm)
    p_xz = extrude(Plane.XZ * xz, amount=2 * br + 4 * f, both=True)
    out = p_yz & p_xz & box(-xm, xm, -big, big, lo, cu + cap) & cyl_z(R, lo - 1.0, top)
    if bridge:
        vault = fillet_vertices(Pos(0, zb / 2) * Rectangle(2 * x0, zb), fb, lambda v: abs(v.Y - zb) < TOL)
        out -= revolve(Plane.XY * vault, Axis.X)
    if tie:
        rect = Rectangle(2 * (R + 2.0), 2 * br)
        sharp = rect + (Circle(R + 2.0) - Circle(rb))
        rnd = fillet_vertices(sharp, fb, lambda v: abs(abs(v.Y) - br) < TOL
                              and abs((v.X ** 2 + v.Y ** 2) ** 0.5 - rb) < ON_CIRCLE_TOL)
        webs = (rnd - sharp) & Rectangle(2 * (R + 2.0), 2 * (br + fb + 1.0))
        out += extrude(Plane.XY.offset(bb) * webs, amount=cu + cap - bb) & cyl_z(R, bb, top)
    return out


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


def piston_dims(p: Params) -> dict:
    c = p.section("piston")
    R = piston_radius(p)
    top = p.compression_height
    spans = groove_spans(p)
    return {"R": R, "top": top, "crown_under": top - p.L(c["crown_thickness"]),
            "belt_bottom": spans[-1][1] - p.L(c["belt_extra_below_rings"]),
            "r_belt": R - p.L(c["ring_belt_thickness"]), "r_groove": R - p.L(c["groove_depth"]),
            "groove1_top": spans[0][0], "boss_r": p.L(c["boss_radius"]), "boss_gap": boss_gap(p),
            "panel_to_crown": bool(c["side_panel_to_crown"]),
            "flat_x": p.L(c["skirt_flat_half_width"]), "pin_r": p.L(p.section("wrist_pin")["outer_diameter"]) / 2}


def piston_blank(p: Params, df: dict | None = None):
    """Piston body shared by the CAD and the FEA (ADR 0006): crown, top land, ring belt (grooves filled), skirt,
    side panels, pillars, pin bore, with the drawing blends of spec/fea.yaml piston.defeature (df: radii in CAD
    units; default from spec/fea.yaml). The CAD piston cuts the ring grooves, circlip grooves and valve pockets
    from it (`piston`); the FEA analyses it as is (quarter)."""
    d = piston_dims(p)
    if df is None:
        df = {k: p.L(v) for k, v in p.fea["piston"]["defeature"].items()}
    R, top, cu, bb, rb, rg, g1 = (d[k] for k in ("R", "top", "crown_under", "belt_bottom", "r_belt", "r_groove",
                                                  "groove1_top"))
    c = p.section("piston")
    zbot = -p.L(c["skirt_below_pin"])
    rs, fx, pt = R - p.L(c["skirt_wall"]), d["flat_x"], p.L(c["side_panel_thickness"])
    # axisymmetric crown + top land + grooved belt + skirt, profile in the (r, z) half plane. Ring grooves
    # filled (the lands between grooves stiffen the belt like hoops; cutting the belt back to the groove root
    # would leave a 1.5 mm wall and grossly understate the crown edge support).
    prof = Polygon((0, top), (R, top), (R, zbot), (rs, zbot), (rs, bb), (rb, bb), (rb, cu), (0, cu), align=None)
    prof = fillet_vertices(prof, df["crown_under_fillet"], lambda v: abs(v.X - rb) < TOL and abs(v.Y - cu) < TOL)
    prof = fillet_vertices(prof, df["skirt_belt_fillet"], lambda v: abs(v.X - rs) < TOL and abs(v.Y - bb) < TOL)
    body = revolve(Plane.XZ * prof, Axis.Z)
    # skirt relief on the pin sides below the belt (cutter edge rounded: blend at the relief corner)
    fr = df["panel_fillet"]
    relief = fillet_vertices(Polygon((fx, zbot - 1.0), (BIG, zbot - 1.0), (BIG, bb), (fx, bb), align=None), fr,
                           lambda v: abs(v.X - fx) < TOL and abs(v.Y - bb) < TOL)
    relief = extrude(Plane.XZ * relief, amount=BIG, both=True)
    body -= relief
    body -= relief.mirror(Plane.YZ)
    # side panels joining the skirt halves at the reliefs; with side_panel_to_crown they run on up to the crown
    # underside (walls tying crown, ring belt, skirt and pillars). Blends: tangent flare into the crown
    # underside (xz profile), planform webs into the skirt bore (below the belt bottom), the ring-belt
    # inner wall and the pillar flanks (above), and a flare under the belt bottom between belt and skirt bore.
    br = d["boss_r"]
    ztop = cu if d["panel_to_crown"] else bb
    pp = Polygon((fx - pt, zbot), (fx, zbot), (fx, ztop + fr), (fx - pt - 2 * fr, ztop + fr), (fx - pt - 2 * fr, ztop),
                 (fx - pt, ztop), align=None)
    if d["panel_to_crown"]:
        pp = fillet_vertices(pp, fr, lambda v: abs(v.X - (fx - pt)) < TOL and abs(v.Y - ztop) < TOL)
    # fillet material only (no caps), so the piece is fully fused to panel face and belt bottom
    under = Pos(fx - pt - fr / 2, bb - fr / 2) * Rectangle(fr, fr) - Pos(fx - pt - fr, bb - fr) * Circle(fr)
    for mir in (False, True):
        pan = extrude(Plane.XZ * pp, amount=BIG, both=True)
        und = extrude(Plane.XZ * under, amount=BIG, both=True) & (cyl_z(rs, bb - 2 * fr, bb + fr)
                                                                   - cyl_z(rb, bb - 2 * fr, bb + fr))
        piece = (pan & cyl_z(R, zbot, ztop + fr)) + und
        body += piece.mirror(Plane.YZ) if mir else piece
    strip = Rectangle(pt, 2 * (R + 2.0))
    strips = Pos(fx - pt / 2, 0) * strip + Pos(-(fx - pt / 2), 0) * strip
    inner_x = lambda v: abs(abs(v.X) - (fx - pt)) < TOL

    def webs(other, pred, clip, radius=fr):
        sharp = strips + other
        return (fillet_vertices(sharp, radius, lambda v: pred(v) and inner_x(v)) - sharp) & clip

    on_r = lambda r: (lambda v: abs((v.X ** 2 + v.Y ** 2) ** 0.5 - r) < 1.0e-3)
    # panel to skirt bore below the belt: its own (larger) blend, df panel_skirt_fillet (ADR 0006)
    w = webs(Circle(R + 2.0) - Circle(rs), on_r(rs), Circle(rs), df.get("panel_skirt_fillet", fr))
    body += extrude(Plane.XY.offset(zbot) * w, amount=bb - zbot)
    if d["panel_to_crown"]:
        w = webs(Circle(R + 2.0) - Circle(rb), on_r(rb), Circle(rb))
        body += extrude(Plane.XY.offset(bb) * w, amount=cu + fr - bb)
    w = webs(Rectangle(2 * fx, 2 * br), lambda v: abs(abs(v.Y) - br) < TOL, Rectangle(2 * fx, 2 * (br + fr + 1.0)))
    body += extrude(Plane.XY * w, amount=ztop + fr)
    # pin-boss pillars, belt ties and crown bridge: the CAD geometry itself (cad.piston.boss_pillars), which
    # already carries its blends
    body += boss_pillars(p, R, top, cu, bb, rb)
    body -= cyl_x(d["pin_r"], -2 * R, 2 * R)
    return body




def piston(p: Params):
    """CAD piston: the shared blank with ring grooves, circlip grooves and valve pockets cut."""
    c = p.section("piston")
    R = piston_radius(p)
    span = 2 * R
    body = piston_blank(p)
    r_groove = R - p.L(c["groove_depth"])
    for zt, zb in groove_spans(p):
        body -= tube_z(span, r_groove, zb, zt)
    for cp in clips(p):  # circlip grooves in the pin bore
        body -= cp
    pk = valve_pockets(p)
    if pk is not None:
        body -= pk
    return body


def valve_pockets(p: Params):
    """Crown valve pockets (spec/valvetrain.yaml combustion_chamber pocket depths), piston-local frame.

    Each pocket is a cylinder on its valve axis, diameter = valve head + valve_relief_margin, whose floor is normal
    to the valve axis at valve_to_crown_at_tdc + pocket_depth below the closed valve face (piston at TDC): exactly
    the gap the M3 clearance model (kin/valvetrain.py valve_pocket_gap) uses. With a flat crown the floor plane
    leaves a crescent, deepest under the valve's outer edge."""
    from cad.valvetrain import valve_frame
    c = p.section("piston")
    if not c["valve_pockets"]:
        return None
    cc = p.vt["combustion_chamber"]
    g0 = p.L(cc["valve_to_crown_at_tdc"])
    out = None
    to_local = Pos(-p.cyl_x(0), 0.0, -p.tdc_pin_z)
    for kind, _sign, x, _y, d in p.valve_layout():
        depth = p.L(cc[f"{kind}_pocket_depth"])
        if depth <= 0.0:
            continue
        r = (d + p.L(c["valve_relief_margin"])) / 2
        cut = to_local * valve_frame(p, 0, kind, x) * cyl_z(r, -(g0 + depth), p.L(c["valve_pocket_cutter_height"]))
        out = cut if out is None else out + cut
    if out is None:
        return None
    R = piston_radius(p)
    top = p.compression_height
    return out & cyl_z(R, top - 2 * R, top)


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


