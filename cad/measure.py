"""Measurements taken from the CAD solids themselves (used by cad.build and the M1 tests)."""
from __future__ import annotations

import math

import numpy as np
from build123d import GeomType

from cad.conrod import big_end_bore_r, small_end_bore_r
from cad.crankshaft import journal_spans, pin_yz, web
from build123d import Pos

from cad.head import valve_cuts
from cad.massprops import MassProps, mass_props
from cad.params import Params
from cad.primitives import cyl_x, cyl_z

RADIUS_MATCH_TOL = 1e-6  # CAD units: picks circle edges of a given nominal radius


def common_volume(a, b) -> float:
    """Volume shared by two solids (CAD units^3). Sums the solids of the boolean common so an
    empty (null) result reads as zero."""
    return float(sum(s.volume for s in (a & b).solids()))


def circle_centres(shape, radius: float) -> list[np.ndarray]:
    return [np.array(tuple(e.arc_center)) for e in shape.edges()
            if e.geom_type == GeomType.CIRCLE and abs(e.radius - radius) < RADIUS_MATCH_TOL]


def conrod_centre_distance(p: Params, rod) -> float:
    """Distance between big-end and small-end bore axes (both along x), CAD units."""
    be = circle_centres(rod, big_end_bore_r(p))
    se = circle_centres(rod, small_end_bore_r(p))
    return float(np.hypot(*(np.mean(se, axis=0) - np.mean(be, axis=0))[1:]))


def crank_throws(p: Params, crank) -> list[float]:
    """Radial offset of every crankpin circle edge from the crank axis (x), CAD units."""
    r_pin = p.L(p.section("crankshaft")["pin_diameter"]) / 2
    return [float(math.hypot(c[1], c[2])) for c in circle_centres(crank, r_pin)]


def clearance_gas_volume(p: Params, head_solid, piston_solid, idx: int = 0) -> float:
    """CAD clearance volume of cylinder idx at TDC (CAD units^3), measured from the actual
    head and piston solids: the bore cylinder from half-way down the top land up to the roof
    ridge, minus the head, minus the piston placed at TDC (so crown reliefs count), minus the
    seat/throat cuts (closed valves fill them flush with the roof). Independent of the
    analytic roof-angle solve in Params. The thin top-land crevice (skirt clearance) inside
    the region is included; it is ~1e-4 of the clearance volume."""
    x = p.cyl_x(idx)
    z0 = p.crown_tdc_z - p.L(p.section("piston")["top_land"]) / 2
    region = cyl_z(p.bore_r, z0, p.head_face_z + p.roof_height, x)
    region -= head_solid
    region -= Pos(x, 0.0, p.crown_tdc_z - p.compression_height) * piston_solid
    region -= valve_cuts(p, idx)
    return region.volume


def throw_solids(p: Params, idx: int, with_counterweights: bool):
    """Crankpin plus its two webs for cylinder idx (crank frame, angle 0)."""
    c = p.section("crankshaft")
    mains = journal_spans(p)
    pw = p.L(c["pin_width"]) / 2
    xc = p.cyl_x(idx)
    py, pz = pin_yz(p, idx)
    solid = cyl_x(p.L(c["pin_diameter"]) / 2, xc - pw, xc + pw, py, pz)
    solid += web(p, idx, mains[idx][1], xc - pw, with_counterweights)
    solid += web(p, idx, xc + pw, mains[idx + 1][0], with_counterweights)
    return solid


def counterweight_moments(p: Params, crank_mp: MassProps, crank_density: float) -> dict:
    """Static balance of the crankshaft, first moments in kg m (radial y-z part).

    Per throw: the counterweight first moment is (throw with counterweights) minus (throw
    without), i.e. the net moment the sectors add beyond the web hull. It is compared with the
    throw's rotating unbalance: pin + web hull + the conrod rotating share (1 - small-end
    fraction) at the crank radius. balance_ratio = 1 means the throw is fully balanced
    locally; the flat-plane layout cancels throw unbalance globally regardless.
    """
    m = p.engine["masses"]
    rod_rot = m["conrod"] * (1 - m["conrod_small_end_fraction"])
    r = p.geom.crank_radius
    throws = {}
    cw_total = np.zeros(3)
    for idx in range(p.cylinders):
        with_cw = mass_props(throw_solids(p, idx, True), crank_density, p.k).first_moment()
        bare = mass_props(throw_solids(p, idx, False), crank_density, p.k).first_moment()
        py, pz = pin_yz(p, idx)
        u = np.array([0.0, py, pz]) / p.crank_r
        unbalance = bare + rod_rot * r * u
        cw = with_cw - bare
        cw_total += cw
        throws[f"cylinder_{idx + 1}"] = {
            "counterweight_first_moment_kgm": cw[1:].tolist(),
            "unbalance_pin_web_rod_rotating_kgm": unbalance[1:].tolist(),
            "balance_ratio": float(-np.dot(cw[1:], unbalance[1:]) / np.dot(unbalance[1:], unbalance[1:])),
        }
    total = crank_mp.first_moment()
    return {
        "note": "first moments are (y, z) in kg m about the crank axis, crank at 0 deg",
        "per_throw": throws,
        "counterweights_total_first_moment_kgm": cw_total[1:].tolist(),
        "crank_first_moment_kgm": total[1:].tolist(),
        "crank_static_imbalance_kgm": float(np.hypot(total[1], total[2])),
    }
