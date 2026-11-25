"""CAD audit gate (cad/audit.py, spec/cad_fits.yaml). All geometry is MODELLED.

The full audit (exact BRep clash check of every part pair, 720-deg / valve-lift sweeps of every moving pair, support
graph, tessellation quality) takes about an hour on 4 cores, so it is run explicitly (`python -m cad.audit --tag final`
then `python -m cad.audit --docs`) and its result is checked here: out/cad/audit_final.json when present, else the
committed compact copy docs/cad-audit.json. The JSON carries the hash of the CAD sources and specs it was produced
from (taken when the run starts), so a stale audit fails. A few live checks re-measure the essentials directly.

Stated tolerances (never loosened silently): contact <= 0.05 mm and overlap <= 0.01 mm^3 (spec/cad_fits.yaml);
tessellation chord <= export.max_chord_mm and facet angle <= export.max_facet_angle_rad (spec/cad.yaml); CAD
piston-to-valve minimum inside the analytic bracket [min(M3 gap, pocket wall), M3 gap] within PV_TOL_MM; exhaust
header duct lengths equal the gasdyn lengths within EXHAUST_LEN_TOL_MM.
"""
from __future__ import annotations

import json

import pytest

from cad.audit import DOCS_JSON, OUT_DIR, Fits, common_volume, distance, source_hash
from cad.params import params

AUDIT = OUT_DIR / "audit_final.json"
PV_TOL_MM = 0.05    # sweep resolution on the piston-to-valve bracket (audit resolution_mm scale)
EXHAUST_LEN_TOL_MM = 1e-3   # header duct centre lines vs gasdyn lengths (solver bisection converges far below this)


@pytest.fixture(scope="module")
def audit():
    found = [json.loads(p.read_text()) for p in (AUDIT, DOCS_JSON) if p.exists()]
    current = [r for r in found if r.get("source_hash") == source_hash()]
    if current or found:
        return (current or found)[0]
    pytest.fail("no audit result: run `python -m cad.audit --tag final && python -m cad.audit --docs`")


def test_audit_is_current(audit):
    assert audit["source_hash"] == source_hash(), "CAD sources changed since the audit: rerun python -m cad.audit"


def test_no_unexplained_overlaps_or_contacts(audit):
    s = audit["summary"]
    assert s["unexplained_overlaps"] == 0, audit["defects"]
    assert s["unexplained_contacts"] == 0, audit["defects"]
    assert s["missing_required_contacts"] == 0, audit["missing_required_contacts"]


def test_nothing_floats(audit):
    assert audit["floating"] == [], audit["floating"]
    assert audit["floating_without_unlisted_contacts"] == []


def test_moving_parts_clear(audit):
    tol = audit["tolerances"]["contact_mm"]
    assert audit["summary"]["moving_clashes"] == 0
    assert audit["summary"]["moving_pairs"] > 0
    for w in audit["sweep"]:
        if w["mode"]:   # allowed running contact: no penetration
            assert w["max_overlap_mm3"] <= audit["tolerances"]["overlap_mm3"] or not w["clash"], w
        else:
            assert w["min_mm"] > tol, w


def test_springs_do_not_bind(audit):
    for kind, v in audit["valve_checks"]["springs"].items():
        assert v["coil_bind_margin_mm"] >= v["min_margin_mm"], (kind, v)
        assert v["ok"], (kind, v)


def test_piston_to_valve_agrees_with_m3(audit):
    """The pocket floor is the M3 plane and the pocket wall is valve_relief_margin/2 from the head all round, so the
    3D minimum over the sweep is bounded analytically: min(M3, wall) <= CAD <= M3, within PV_TOL_MM."""
    wall = params().L(params().section("piston")["valve_relief_margin"]) / 2
    for key, v in audit["valve_checks"]["piston_valve"].items():
        assert v["pairs"] > 0, key
        assert v["cad_min_mm"] > audit["tolerances"]["contact_mm"], (key, v)
        assert v["cad_min_mm"] >= min(v["m3_min_mm"], wall) - PV_TOL_MM, (key, v)
        assert v["cad_min_mm"] <= v["m3_min_mm"] + PV_TOL_MM, (key, v)


def test_key_clearances_reported_and_positive(audit):
    """Every report_clearances pair (exhaust header to pan, flywheel, block, intake; intake to cams and belt) is
    measured and clear."""
    listed = {tuple(e["parts"]) for e in Fits().report}
    got = {tuple(k["parts"]): k for k in audit["key_clearances"]}
    assert listed and listed <= set(got), listed - set(got)
    for k in got.values():
        assert k["min_mm"] > audit["tolerances"]["contact_mm"], k


def test_intake_runner_length_matches_gasdyn():
    """Head port + flange + CAD runner centre line (flange face to bellmouth lip) = spec/gasdyn.yaml runner length,
    and the plenum's inner volume = the gasdyn plenum volume."""
    from cad import manifolds as M
    from cad.head import port_length
    p = params()
    total = port_length(p, "intake") + p.L(p.section("manifolds")["flange_thickness"]) + M.runner_length(p)
    assert total == pytest.approx(p.L(p.gas["intake"]["runner"]["length"]), abs=1e-6)
    g = M.plenum_geom(p)
    vol = M._rbox(g["x0"], g["x1"], g["y0"], g["y1"], g["z0"], g["z1"], g["rf"]).volume
    assert vol == pytest.approx(p.gas["intake"]["plenum"]["volume"] * p.k**3, rel=1e-6)


def test_exhaust_header_matches_gasdyn_and_packaging():
    """Every primary (head flange face to merge-cone mouth) and secondary (cone outlet to collector mouth) centre line
    measured on the CAD wire equals the spec/gasdyn.yaml length (EXHAUST_LEN_TOL_MM; the closed-form arc length used
    by the solver agrees), every bend radius is >= min_bend_to_od x the tube OD, the header is one solid, and its
    lowest point is within max_drop_below_pan_floor of the sump floor (no header hanging far below the engine)."""
    from cad import manifolds as M
    from cad.block import pan_floor_z
    p = params()
    g, me, w = p.gas["exhaust"], p.section("manifolds")["exhaust"], p.L(p.section("manifolds")["wall"])
    R, r = M.bend_radii(p), M.exhaust_radii(p)
    for i in range(p.cylinders):
        L = M.primary_length(p, i)
        assert L == pytest.approx(p.L(g["primary"]["length"]), abs=EXHAUST_LEN_TOL_MM), (i, L)
        assert M.poly_length(M.primary_points(p, i), R["primary"]) == pytest.approx(L, abs=EXHAUST_LEN_TOL_MM)
    for pair in (0, 1):
        L = M.secondary_length(p, pair)
        assert L == pytest.approx(p.L(g["secondary"]["length"]), abs=EXHAUST_LEN_TOL_MM), (pair, L)
        assert M.poly_length(M.secondary_points(p, pair), R["secondary"]) == pytest.approx(L, abs=EXHAUST_LEN_TOL_MM)
    for k in ("primary", "secondary"):
        assert R[k] >= float(me["min_bend_to_od"]) * 2.0 * (r[k] + w), k
    h = M.exhaust_header(p)
    assert len(h.solids()) == 1
    assert h.bounding_box().min.Z >= pan_floor_z(p) - p.L(me["max_drop_below_pan_floor"])


def test_tessellation_within_limits(audit):
    t = audit["tessellation"]
    assert t["max_chord_mm"] <= t["limits"]["max_chord_mm"]
    assert t["max_facet_angle_rad"] <= t["limits"]["max_facet_angle_rad"]
    assert audit["summary"]["tessellation_fail"] == 0


def test_fits_reasons_present():
    f = Fits()
    assert f.entries
    for e in f.entries:
        assert e.get("reason"), e
        assert e["type"] in ("contact", "overlap")


def test_live_tessellation_of_curved_part():
    """Fresh mesh of the crankshaft at the export settings: exact chordal deviation and facet angle in limits."""
    from cad.crankshaft import crankshaft
    from cad.tess import triangulate
    p = params()
    ex = p.section("export")
    _m, q = triangulate(crankshaft(p), p.L(ex["linear_deflection"]), float(ex["angular_deflection_rad"]), measure=True)
    assert q["chord_max_mm"] <= float(ex["max_chord_mm"])
    assert q["facet_angle_max_rad"] <= float(ex["max_facet_angle_rad"])


def test_live_crank_clears_block_and_caps():
    """Baseline defect (crank web faces flush with the bulkheads) stays fixed: real axial clearance at the mains."""
    from cad.block import block
    from cad.crankshaft import crankshaft
    from cad.maincap import main_cap_solid
    from build123d import Pos
    p = params()
    crank, blk = crankshaft(p), block(p)
    assert common_volume(crank.wrapped, blk.wrapped) == 0.0
    d = distance(crank.wrapped, blk.wrapped)[0]
    assert d > Fits().tol["contact_mm"]
    cap = Pos(p.main_x(2), 0, 0) * main_cap_solid(p, half=False)
    assert distance(crank.wrapped, cap.wrapped)[0] > Fits().tol["contact_mm"]
