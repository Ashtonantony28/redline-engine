"""M1 exit criteria: parametric CAD mass properties, balance, kinematic dimensions and
clearances (PLAN.md sec. 5). All values are modelled.

Tolerances (stated, never loosened silently): masses and crank inertia within 3% of the
spec budget; clearance volume within 5% of the geometric value from the spec CR; conrod
centre distance and crank throw exact to CAD float precision (EXACT_TOL_M).
"""
from __future__ import annotations

import json
import math

import pytest

from cad.assembly import cylinder_pose, parts, placed_parts
from cad.build import main as build_main
from cad.build import mass_report
from cad.measure import common_volume
from cad.params import params
from kin.slider_crank import piston_position

BUDGET_TOL = 0.03
CLEARANCE_VOL_TOL = 0.05
EXACT_TOL_M = 1e-9
STATIC_IMBALANCE_TOL_KGM = 1e-9   # flat-plane symmetry: zero up to float noise
INTERFERENCE_TOL_MM3 = 1e-6       # boolean common volume treated as zero below this
DISTANCE_REL_TOL = 1e-3          # BRepExtrema distance solver tolerance on the crown gap
CHECK_ANGLES_DEG = (0.0, 45.0, 90.0, 135.0, 180.0, 270.0)


@pytest.fixture(scope="module")
def p():
    return params()


@pytest.fixture(scope="module")
def report(p):
    return mass_report(p)


def test_piston_assembly_mass_within_budget(report):
    pa = report["piston_assembly"]
    assert abs(pa["rel_err"]) <= BUDGET_TOL, pa


def test_conrod_mass_within_budget(report):
    rod = report["conrod"]
    assert abs(rod["rel_err"]) <= BUDGET_TOL, rod


def test_reciprocating_mass_within_budget(report):
    assert abs(report["reciprocating_mass"]["rel_err"]) <= BUDGET_TOL


def test_crank_inertia_within_budget(report):
    cr = report["crankshaft"]
    assert abs(cr["rel_err"]) <= BUDGET_TOL, cr["inertia_about_crank_axis_kgm2"]


def test_crank_static_balance_and_counterweight_moment_reported(report):
    bal = report["crankshaft"]["balance"]
    assert bal["crank_static_imbalance_kgm"] < STATIC_IMBALANCE_TOL_KGM
    assert len(bal["per_throw"]) == params().cylinders
    for throw in bal["per_throw"].values():
        cw = math.hypot(*throw["counterweight_first_moment_kgm"])
        assert cw > 0.0
        assert 0.0 < throw["balance_ratio"] <= 1.0  # counterweights oppose the throw unbalance


def test_conrod_centre_distance_exact(p, report):
    rod = report["conrod"]
    assert rod["centre_distance_m"] == pytest.approx(p.geom.conrod_length, abs=EXACT_TOL_M)


def test_crank_throw_exact(p, report):
    t = report["crankshaft"]["throw_m"]
    assert t["min"] == pytest.approx(p.geom.crank_radius, abs=EXACT_TOL_M)
    assert t["max"] == pytest.approx(p.geom.crank_radius, abs=EXACT_TOL_M)


def test_flat_plane_pin_phasing(p):
    """Crankpins at 0/180/180/0 deg for firing order 1-3-4-2 (ADR 0001)."""
    assert [p.pin_angle_deg(i) for i in range(p.cylinders)] == [0.0, 180.0, 180.0, 0.0]


def test_clearance_volume_matches_compression_ratio(report):
    ch = report["combustion_chamber"]
    assert ch["clearance_volume_cad_m3"] > 0.0
    assert abs(ch["rel_err"]) <= CLEARANCE_VOL_TOL, ch


def test_no_piston_head_interference_at_tdc(p):
    """Every piston at its TDC: zero common volume with the head and a positive gap."""
    for idx in range(p.cylinders):
        crank_deg = p.firing_phase_deg(idx)
        pp = placed_parts(crank_deg, p)
        head, piston = pp["head"][0], pp[f"piston_{idx + 1}"][0]
        assert common_volume(head, piston) < INTERFERENCE_TOL_MM3
        gap = head.distance_to(piston)
        assert gap > 0.0
        assert gap == pytest.approx(p.head_face_z - p.crown_tdc_z, rel=DISTANCE_REL_TOL)


def test_valve_seats_fit_inside_bore(p):
    for _kind, _sign, x, y, d in p.valve_layout():
        assert math.hypot(x, y) + d / 2 < p.bore_r
        assert abs(y) > d / 2 * math.cos(math.radians(p.roof_angle_deg))  # clear of the ridge


def test_piston_follows_slider_crank(p):
    for idx in range(p.cylinders):
        for deg in CHECK_ANGLES_DEG:
            _py, _pz, _rod, wz = cylinder_pose(p, idx, deg)
            a = math.radians(p.pin_angle_deg(idx, deg))
            assert wz == pytest.approx(float(piston_position(a, p.crank_r, p.rod_len)))


@pytest.mark.parametrize("crank_deg", CHECK_ANGLES_DEG)
def test_moving_parts_clear_block(p, crank_deg):
    """Crank, conrods and pistons have no common volume with the block over a revolution."""
    pp = placed_parts(crank_deg, p)
    block = pp["block"][0]
    movers = ["crankshaft"] + [f"{n}_{i + 1}" for i in range(p.cylinders) for n in ("conrod", "piston")]
    for name in movers:
        assert common_volume(block, pp[name][0]) < INTERFERENCE_TOL_MM3, name


@pytest.mark.parametrize("crank_deg", CHECK_ANGLES_DEG)
def test_conrods_and_pistons_clear_crank(p, crank_deg):
    """Big ends clear the webs; counterweights clear the piston skirts near BDC."""
    pp = placed_parts(crank_deg, p)
    crank = pp["crankshaft"][0]
    for i in range(p.cylinders):
        for n in ("conrod", "piston"):
            assert common_volume(crank, pp[f"{n}_{i + 1}"][0]) < INTERFERENCE_TOL_MM3, (n, i)


def test_build_writes_exports_and_mass_properties(tmp_path):
    assert build_main(["--out", str(tmp_path)]) == 0
    for name in parts():
        assert (tmp_path / f"{name}.step").stat().st_size > 0
        assert (tmp_path / f"{name}.stl").stat().st_size > 0
    assert (tmp_path / "assembly_0deg.step").exists()
    data = json.loads((tmp_path / "mass_properties.json").read_text())
    assert "inertia_about_crank_axis_kgm2" in data["parts"]["crankshaft"]
    assert len(data["parts"]["crankshaft"]["inertia_about_crank_axis_kgm2"]) == 3
