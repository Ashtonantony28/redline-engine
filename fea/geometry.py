"""Meshing-ready solids (build123d, CAD units = mm) for the FEA, exported to STEP and meshed by fea.mesh.

The M1 CAD solids (cad/) are built for mass properties and clearances: they have sharp re-entrant
corners (shank-to-eye, pocket corners, pin-to-web) where elastic peak stress is singular, and their
boolean topology does not fillet reliably. Each part is therefore reconstructed here from the SAME
spec/cad.yaml parameters with the blend radii of spec/fea.yaml (or spec/cad.yaml for the crank fillets)
and the documented defeaturing. `cad_volume_check` compares each reconstruction with the CAD solid.
Symmetric parts are cut to the symmetric portion that is analysed (conrod and piston: quarter).
"""
from __future__ import annotations

from pathlib import Path

from build123d import (Box, Circle, Plane, Polygon, Pos, Rectangle, Unit, export_step, extrude, fillet,
                       revolve, Axis)

from cad.conrod import big_end_bore_r, small_end_bore_r
from cad.crankshaft import web
from cad.params import Params
from cad.piston import piston_blank, piston_dims  # noqa: F401  (piston_dims re-exported for fea.parts.piston)
from cad.primitives import box, cyl_x, cyl_z

TOL = 1.0e-6          # CAD units, vertex / face position matching
BIG = 1.0e3           # CAD units, extent of half-space cutting boxes
MIN_STEP = 0.05       # CAD units, eye overhang below which the flange is flush with the eyes
FLANK_FILLET_MAX_SHARE = 0.95   # flank fillet must stay within the eye overhang


def _halfspace(x=None, y=None, z=None):
    """Box covering coordinate >= 0 for each axis given as True (CAD units)."""
    x0, x1 = (0.0, BIG) if x else (-BIG, BIG)
    y0, y1 = (0.0, BIG) if y else (-BIG, BIG)
    z0, z1 = (0.0, BIG) if z else (-BIG, BIG)
    return box(x0, x1, y0, y1, z0, z1)


def _yz(sketch):
    """Place an (X = y, Y = z) sketch on the global YZ plane."""
    return Plane.YZ * sketch


def _fillet_where(sk, radius, pred):
    vs = [v for v in sk.vertices() if pred(v)]
    return fillet(vs, radius) if vs else sk


def _tool_with_floor_fillet(sk2d, x_floor, x_far, radius):
    """Prism of a (y, z) sketch from x_floor to x_far with its floor edges (at x_floor) filleted: removed
    from a part it leaves concave fillets of `radius` at the floor."""
    tool = extrude(Plane.YZ.offset(x_floor) * sk2d, amount=x_far - x_floor)
    floor = [e for e in tool.edges() if abs(e.center().X - x_floor) < TOL]
    return fillet(floor, radius)


# ---------------------------------------------------------------------------------------------
# conrod
# ---------------------------------------------------------------------------------------------
def conrod_dims(p: Params) -> dict:
    c = p.section("conrod")
    be_r = big_end_bore_r(p)
    se_r = small_end_bore_r(p)
    return {
        "L": p.rod_len, "be_r": be_r, "be_R": be_r + p.L(c["big_end_wall"]), "be_w": p.L(c["big_end_width"]),
        "se_r": se_r, "se_R": se_r + p.L(c["small_end_wall"]), "se_w": p.L(c["small_end_width"]),
        "boss_y": be_r + p.L(c["bolt_boss_width"]), "boss_h": p.L(c["bolt_boss_height"]),
        "hb": p.L(c["shank_width_big_end"]) / 2, "hs": p.L(c["shank_width_small_end"]) / 2,
        "fw": p.L(c["flange_width"]) / 2, "ft": p.L(c["flange_thickness"]), "tw": p.L(c["web_thickness"]) / 2,
    }


def conrod_solid(p: Params, fd: dict, quarter: bool = True):
    """Conrod in the CAD local frame (big end at the origin, small end at z = L, bores along x)."""
    d = conrod_dims(p)
    df = {k: p.L(v) for k, v in fd["defeature"].items()}
    if abs(d["be_w"] - d["se_w"]) > TOL:
        raise ValueError("FEA conrod reconstruction assumes equal big- and small-end widths")
    L, hb, hs = d["L"], d["hb"], d["hs"]
    boss_half = d["boss_h"] / 2
    eyeboss = _fillet_where(Circle(d["be_R"]) + Rectangle(2 * d["boss_y"], d["boss_h"]), df["boss_eye_blend"],
                            lambda v: abs(abs(v.Y) - boss_half) < TOL)
    prof = eyeboss + Pos(0, L) * Circle(d["se_R"]) + Polygon((-hb, 0), (hb, 0), (hs, L), (-hs, L), align=None)
    prof = _fillet_where(prof, df["eye_shank_blend"], lambda v: abs(v.Y) > boss_half + TOL)
    W = d["be_w"] / 2
    rod = extrude(_yz(prof), amount=W, both=True)
    # flanks: shank narrower (flange width) than the eyes
    Y = d["boss_y"] + df["eye_shank_blend"]
    band = Pos(0, L / 2) * Rectangle(2 * Y, L) - eyeboss - Pos(0, L) * Circle(d["se_R"])
    far = W + 1.0
    step_w = W - d["fw"]                       # eye overhang beyond the flange face
    if step_w > MIN_STEP:
        r_flank = min(df["flank_fillet"], step_w * FLANK_FILLET_MAX_SHARE)
        tool = _tool_with_floor_fillet(band, d["fw"], far, r_flank)
        rod -= tool
        rod -= tool.mirror(Plane.YZ)
    # I-beam pockets
    taper = (hs - hb) / L
    z0, z1 = d["be_R"] + df["pocket_end_gap"], L - d["se_R"] - df["pocket_end_gap"]
    h0, h1 = hb + taper * z0 - d["ft"], hb + taper * z1 - d["ft"]
    pk = fillet(Polygon((-h0, z0), (h0, z0), (h1, z1), (-h1, z1), align=None).vertices(), df["pocket_corner"])
    tool = _tool_with_floor_fillet(pk, d["tw"], far, df["web_fillet"])
    rod -= tool
    rod -= tool.mirror(Plane.YZ)
    rod -= cyl_x(d["be_r"], -2 * W, 2 * W)
    rod -= cyl_x(d["se_r"], -2 * W, 2 * W, z=L)
    if quarter:
        rod = rod & _halfspace(x=True, y=True)
    return rod


# ---------------------------------------------------------------------------------------------
# crankshaft single throw
# ---------------------------------------------------------------------------------------------
def throw_dims(p: Params) -> dict:
    c = p.section("crankshaft")
    mw = p.L(c["main_width"]) / 2
    return {"half_span": p.spacing / 2, "main_r": p.L(c["main_diameter"]) / 2, "pin_r": p.L(c["pin_diameter"]) / 2,
            "pw": p.L(c["pin_width"]) / 2, "web_x0": p.spacing / 2 - mw, "r_pin_f": p.L(c["pin_fillet_radius"]),
            "r_main_f": p.L(c["main_fillet_radius"]), "crank_r": p.crank_r,
            "web_pin_R": p.L(c["pin_diameter"]) / 2 + p.L(c["web_pin_margin"]),
            "web_main_R": p.L(c["main_diameter"]) / 2 + p.L(c["web_main_margin"])}


def crank_throw_solid(p: Params):
    """Throw of cylinder 1 (pin at +z) centred on x = 0, from main-bearing centre plane to centre plane."""
    d = throw_dims(p)
    S, wx, pw = d["half_span"], d["web_x0"], d["pw"]
    shaft = cyl_x(d["main_r"], -S, -wx) + cyl_x(d["main_r"], wx, S)
    shaft += cyl_x(d["pin_r"], -pw, pw, 0.0, d["crank_r"])
    shaft += web(p, 0, -wx, -pw, True) + web(p, 0, pw, wx, True)
    # fillets: circular junction edges of the pin and journals with the web faces
    edges = []
    for e in shaft.edges():
        c = e.center()
        if e.geom_type.name != "CIRCLE":
            continue
        if abs(abs(c.X) - pw) < TOL and abs(e.radius - d["pin_r"]) < TOL:
            edges.append(("pin", e))
        elif abs(abs(c.X) - wx) < TOL and abs(e.radius - d["main_r"]) < TOL:
            edges.append(("main", e))
    shaft = fillet([e for k, e in edges if k == "pin"], d["r_pin_f"])
    main_edges = [e for e in shaft.edges() if e.geom_type.name == "CIRCLE"
                  and abs(abs(e.center().X) - wx) < TOL and abs(e.radius - d["main_r"]) < TOL]
    shaft = fillet(main_edges, d["r_main_f"])
    return shaft


# ---------------------------------------------------------------------------------------------
# piston (quarter)
# ---------------------------------------------------------------------------------------------
def piston_solid(p: Params, fd: dict, quarter: bool = True):
    """Defeatured piston (see spec/fea.yaml piston), pin axis along x at z = 0, crown up. The body is
    cad.piston.piston_blank, the same solid the CAD piston is cut from (ring grooves, circlip grooves and valve
    pockets are the CAD-only features; see spec/fea.yaml piston for the defeaturing)."""
    body = piston_blank(p, {k: p.L(v) for k, v in fd["defeature"].items()})
    if quarter:
        body = body & _halfspace(x=True, y=True)
    return body


# ---------------------------------------------------------------------------------------------
# main bearing cap (half: y >= 0)
# ---------------------------------------------------------------------------------------------
# geometry lives in cad/maincap.py (one definition for the CAD assembly and the FEA)
from cad.maincap import main_cap_dims, main_cap_solid  # noqa: E402,F401


# ---------------------------------------------------------------------------------------------
# analytic anchors
# ---------------------------------------------------------------------------------------------
def cantilever(L, b, h):
    return box(0.0, L, -b / 2, b / 2, -h / 2, h / 2)


def lame_quarter(ri, ro, t):
    """Quarter of a thick cylinder (x, y >= 0), axis along z, length t."""
    return (cyl_z(ro, 0.0, t) - cyl_z(ri, 0.0, t)) & _halfspace(x=True, y=True)


def plate_hole_quarter(W, Hh, r, t):
    """Quarter of a plate 2W x 2H with a central hole of radius r, thickness t (x, y >= 0)."""
    return box(0.0, W, 0.0, Hh, 0.0, t) - cyl_z(r, -t, 2 * t)


def slab(L, a):
    return box(0.0, L, 0.0, a, 0.0, a)


def export(shape, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    export_step(shape, str(path), unit=Unit.MM)
    return path
