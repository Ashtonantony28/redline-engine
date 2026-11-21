"""M6: structural FEA (gmsh C3D10 + CalculiX). Analytic anchors first (CLAUDE.md rule 2), then the
load-case and fatigue helpers, then the M6 exit criteria (PLAN.md sec. 5). All values are modelled.

Tolerances (stated, never loosened silently):
- cantilever tip deflection and bending stress vs Euler-Bernoulli: 2% (L/h = 20, quadratic tets)
- Lame thick cylinder radial / hoop stress at every node: 2% of the inner-wall hoop stress
- plate with a hole Kt vs the Heywood finite-width value (3.03 at d/W = 0.1): 3%
- 1-D conduction slab: 0.01 K at every node (quadratic elements reproduce a linear field exactly)
- overrun tension vs m_rec r w^2 (1 + lambda): 1e-9 relative (closed form at TDC with zero gas force)
- cooling-jet hydraulics, pump FMEP and Stevens & Webb coefficient vs closed forms: 1e-12 relative
- piston thermal consistency: FEA area-mean crown temperature vs thermo.yaml T_piston within
  spec/fea.yaml piston.thermal.consistency_tol_K (3 K); thermal-solve heat balance within 1%
- exit criteria: every part meets its SF target (conrod buckling >= 3, fatigue >= 1.5) on the finest
  mesh and the peak stress changes < 5% between the two finest meshes.

The exit-criteria tests read out/fea/results.json written by `python -m fea.run --all`; they are skipped
when that file is missing or stale (spec hash differs) unless REDLINE_RUN_SLOW=1, which runs the full
study (~15-25 min). The anchors and unit tests run in well under a minute.
"""
from __future__ import annotations

import json
import math
import os
import shutil

import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
if shutil.which("ccx") is None and not os.environ.get("REDLINE_CCX"):
    pytest.skip("CalculiX (ccx) not installed", allow_module_level=True)

from fea import anchors  # noqa: E402
from fea.config import base_params, fea_spec, spec_hash  # noqa: E402
from fea.fatigue import goodman_sf, johnson_euler, peterson_kt_bending, richardson, signed_von_mises  # noqa: E402
from fea.geometry import cantilever, export  # noqa: E402
from fea.mesh import mesh_step  # noqa: E402
from redline.spec import ROOT  # noqa: E402

BEAM_TOL = 0.02
LAME_TOL = 0.02
KT_TOL = 0.03
SLAB_TOL_K = 0.01
EXACT_REL = 1e-9
CONV_TOL = 0.05


@pytest.fixture(scope="module")
def wd(tmp_path_factory):
    return tmp_path_factory.mktemp("fea_anchors")


# ---------------------------------------------------------------------------------------------
# analytic anchors
# ---------------------------------------------------------------------------------------------
def test_cantilever_matches_euler_bernoulli(wd):
    r = anchors.cantilever_case(wd)
    assert abs(r["tip_rel_err"]) < BEAM_TOL, r
    assert abs(r["stress_rel_err"]) < BEAM_TOL, r


def test_lame_thick_cylinder(wd):
    r = anchors.lame_case(wd)
    assert r["hoop_max_err_rel"] < LAME_TOL, r
    assert r["radial_max_err_rel"] < LAME_TOL, r


def test_plate_with_hole_kt_about_three(wd):
    r = anchors.plate_hole_case(wd)
    assert abs(r["rel_err"]) < KT_TOL, r
    assert 2.9 < r["kt_fea"] < 3.2


def test_slab_conduction_linear_profile(wd):
    r = anchors.slab_case(wd)
    assert r["max_err_K"] < SLAB_TOL_K, r


# ---------------------------------------------------------------------------------------------
# mesh bookkeeping
# ---------------------------------------------------------------------------------------------
def test_mesh_mapping_and_outward_normals(wd):
    step = export(cantilever(40.0, 10.0, 10.0), wd / "box.step")
    m = mesh_step(step, 1.0e-3, 4.0, 0.4, 0)
    assert m.check_midside() < 1e-9               # straight box: every mid-side node exactly mid-edge
    _e, _f, _cn, _six, c, n, a = m.boundary_faces()
    # divergence theorem: closed surface, outward normals -> sum(n a) = 0 and sum(c.n a)/3 = V
    assert np.abs((n * a[:, None]).sum(axis=0)).max() < 1e-12
    assert (np.einsum("ij,ij->i", c, n) * a).sum() / 3 == pytest.approx(m.volume(), rel=1e-9)
    assert m.volume() == pytest.approx(40e-3 * 10e-3 * 10e-3, rel=1e-9)


# ---------------------------------------------------------------------------------------------
# fatigue / buckling helpers
# ---------------------------------------------------------------------------------------------
def test_goodman_closed_form():
    se, sut = 400.0e+6, 1000.0e+6
    assert goodman_sf(100.0e+6, 250.0e+6, se, sut) == pytest.approx(1.0 / (0.25 + 0.25))
    assert goodman_sf(100.0e+6, -300.0e+6, se, sut) == pytest.approx(4.0)   # compressive mean not credited


def test_signed_von_mises_uniaxial():
    s = np.array([[200.0, 0, 0, 0, 0, 0], [-150.0, 0, 0, 0, 0, 0]])
    assert signed_von_mises(s) == pytest.approx([200.0, -150.0])


def test_johnson_meets_euler_at_transition():
    E, sy, A, I = 210.0e+9, 850.0e+6, 2.0e-4, 3.0e-9
    rho = math.sqrt(I / A)
    cc = math.sqrt(2 * math.pi**2 * E / sy)
    Le = cc * rho
    pj, reg_j, *_ = johnson_euler(A, I, Le * (1 - 1e-9), E, sy)
    pe, reg_e, *_ = johnson_euler(A, I, Le * (1 + 1e-9), E, sy)
    assert (reg_j, reg_e) == ("Johnson", "Euler")
    assert pj == pytest.approx(pe, rel=1e-6)
    assert pj == pytest.approx(A * sy / 2, rel=1e-6)


def test_richardson_recovers_second_order():
    h = np.array([1.0, 0.5, 0.25])
    f = 10.0 + 3.0 * h**2
    p, fe, _ = richardson(f, 2.0)
    assert p == pytest.approx(2.0)
    assert fe == pytest.approx(10.0)


def test_peterson_fit_sanity():
    # stepped round bar in bending, D/d = 1.2, r/d = 0.1: chart value ~1.6 (Shigley Fig. A-15-9)
    assert peterson_kt_bending(1.2, 1.0, 0.1) == pytest.approx(1.62, abs=0.06)


# ---------------------------------------------------------------------------------------------
# piston alloy data (spec/engine.yaml al_2618_t61) reproduces its cited sources
# ---------------------------------------------------------------------------------------------
BASQUIN_RT = (5.3e+33, 12.7)        # Han et al. 2026, Table 4: 1000 h / 230 C overaged, tested at RT
BASQUIN_230C = (6.2e+20, 7.3)       # same, tested at 230 C, 20 Hz
T_230C_K = 503.15
TABLE_TOL_PA = 0.1e+6               # tables are rounded to 0.1 MPa


def test_piston_alloy_tables_match_citations():
    from fea.config import material, strength_at
    p = base_params()
    m = material(p.engine, p.cad["materials"]["piston"])
    n = m["fatigue_life_cycles"]
    se_rt, se_230 = ((c / n) ** (1.0 / k) * 1.0e+6 for c, k in (BASQUIN_RT, BASQUIN_230C))   # MPa -> Pa
    sut, se = strength_at([293.15, T_230C_K], m)
    assert se[0] == pytest.approx(se_rt, abs=TABLE_TOL_PA)
    assert se[1] == pytest.approx(se_230, abs=TABLE_TOL_PA)
    assert sut[0] == pytest.approx(351.0e+6)                        # Han Table 3, overaged, RT
    T = m["T_K"]
    assert np.all(np.diff(T) > 0) and np.all(np.diff(m["sut_T"]) < 0) and np.all(np.diff(m["se_T"]) < 0)
    # above 230 C the fatigue strength follows the ASM tensile ratio (documented assumption)
    hot = T > T_230C_K
    ratio = m["sut_T"][hot] / m["sut_T"][T == T_230C_K]
    assert m["se_T"][hot] == pytest.approx(se[1] * ratio, abs=TABLE_TOL_PA)


# ---------------------------------------------------------------------------------------------
# piston cooling jets (thermo/oiljet.py, ADR 0006): closed forms before the FEA uses them
# ---------------------------------------------------------------------------------------------
JET_REL = 1e-12


def test_jet_hydraulics_and_pump_fmep_closed_form():
    from redline.spec import geometry
    from thermo.oiljet import fmep_jets, jet_spec, jet_state, pump_power
    js, st = jet_spec(), jet_state()
    rho, dp, d = js["oil"]["rho"], js["supply_pressure"], js["orifice_diameter"]
    q = js["discharge_coefficient"] * math.pi * d**2 / 4 * math.sqrt(2 * dp / rho)
    assert st.q_orifice == pytest.approx(q, rel=JET_REL)
    assert st.re == pytest.approx(rho * st.u * d / st.mu, rel=JET_REL)
    g = geometry()
    n = js["orifices_per_cylinder"] * g.cylinders
    rpm = 6000.0
    p_w = n * q * dp / js["pump_efficiency"]
    assert pump_power(rpm) == pytest.approx(p_w, rel=JET_REL)
    assert fmep_jets(rpm) == pytest.approx(p_w / (g.displacement * rpm / 120.0), rel=JET_REL)   # 4-stroke
    assert fmep_jets(0.5 * js["on_above_rpm"]) == 0.0                                            # check valve shut


def test_stevens_webb_jet_coefficient():
    """Stevens & Webb (1991) as restated by Varghese & Agarwal (SAE 2004-28-0061): stagnation value by hand,
    within the published Re range; the radial factor is ~1 at the stagnation point and decays outward."""
    from thermo.oiljet import jet_state, local_h, stagnation_nusselt
    from thermo.oiljet import jet_spec
    corr = fea_spec()["piston"]["thermal"]["jets"]["correlation"]
    st, z = jet_state(), jet_spec()["nozzle_to_crown_mean"]
    assert corr["re_range"][0] <= st.re <= corr["re_range"][1]
    nu0 = 2.67 * st.re**0.567 * st.pr**0.4 * (z / st.d) ** -0.0336 * (st.u / st.d) ** -0.237
    assert stagnation_nusselt(st, z, corr) == pytest.approx(nu0, rel=JET_REL)
    h0 = nu0 * st.k / st.d
    h = local_h(np.array([0.0, 2.0, 5.0, 10.0]) * st.d, st, z, corr)
    assert h[0] == pytest.approx(h0 * (1 + 1.13**-9) ** (-1 / 9), rel=JET_REL)
    assert np.all(np.diff(h) < 0)
    assert h[-1] == pytest.approx(h0 * 1.13 * math.exp(-2.3), rel=1e-3)      # far field ~ f(r/d)


def test_cooling_jet_clears_conrod_small_end():
    from thermo.oiljet import jet_spec
    p = base_params()
    js = jet_spec()
    r_se = p.cad["wrist_pin"]["outer_diameter"] / 2 + p.cad["conrod"]["small_end_wall"]
    assert js["target_offset"] - js["orifice_diameter"] / 2 - r_se >= js["min_small_end_clearance"]


# ---------------------------------------------------------------------------------------------
# load cases come from thermo / dyn / CAD
# ---------------------------------------------------------------------------------------------
def test_overrun_tension_is_pure_inertia_at_tdc():
    from fea.loads import conrod_loads
    p = base_params()
    ld = conrod_loads(p)
    m = ld["masses"]["reciprocating"]
    r, l = p.geom.crank_radius, p.geom.conrod_length
    w = ld["redline_rpm"] * 2 * math.pi / 60
    assert ld["tension_theta_deg"] in (0.0, 360.0)
    assert ld["tension_N"] == pytest.approx(m * r * w**2 * (1 + r / l), rel=EXACT_REL)
    assert ld["compression_N"] > ld["tension_N"] > 0.0


def test_peak_firing_load_from_thermo_sweep():
    from fea.loads import conrod_loads, piston_loads, wot_sweep
    pts = wot_sweep()
    assert all(pt.metrics["converged"] for pt in pts)
    pl = piston_loads()
    assert pl["p_max_Pa"] == max(pt.p_max for pt in pts)
    cc = conrod_loads(base_params())["compression_case"]
    assert cc["F_rod_max_N"] == pytest.approx(cc["F_gas_N"] + cc["F_inertia_N"], rel=0.01)  # cos(beta) ~ 1


# ---------------------------------------------------------------------------------------------
# M6 exit criteria
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def results():
    path = ROOT / fea_spec()["solver"]["outdir"] / "results.json"
    parts = {"conrod", "crankshaft", "piston", "main_cap"}
    fresh = False
    if path.exists():
        d = json.loads(path.read_text())
        fresh = d.get("spec_hash") == spec_hash() and set(d.get("parts", {})) == parts
    if not fresh:
        if os.environ.get("REDLINE_RUN_SLOW") != "1":
            pytest.skip("out/fea/results.json missing or stale: run `python -m fea.run --all` "
                        "(or set REDLINE_RUN_SLOW=1)")
        from fea.run import main
        main(["--all"])
    data = json.loads(path.read_text())
    assert set(data["parts"]) == parts
    return data


@pytest.mark.slow
@pytest.mark.parametrize("part", ["conrod", "crankshaft", "piston", "main_cap"])
def test_part_meets_sf_targets(results, part):
    r = results["parts"][part]
    for name, target in r["targets"].items():
        assert r["sf"][name] >= target, (part, name, r["sf"][name], target)


def test_conrod_buckling_target_named(results):
    assert results["parts"]["conrod"]["targets"]["buckling"] == fea_spec()["targets"]["buckling_sf"]


@pytest.mark.slow
@pytest.mark.parametrize("part", ["conrod", "crankshaft", "piston", "main_cap"])
def test_mesh_convergence(results, part):
    c = results["parts"][part]["convergence"]
    assert len(c["peaks_Pa"]) == fea_spec()["convergence"]["levels"] == 3
    assert c["nodes"][0] < c["nodes"][1] < c["nodes"][2]
    assert c["change_finest_two"] < CONV_TOL, c


@pytest.mark.slow
def test_mass_budget_kept_by_iterations(results):
    for part, r in results["parts"].items():
        if r["mass_final"] is not None:
            assert r["mass_ok"], (part, r["mass_final"])


@pytest.mark.slow
def test_piston_thermal_state_consistent_with_thermo(results):
    """The FEA area-mean crown temperature equals thermo.yaml T_piston (the wall temperature behind the
    Woschni flux and the cycle-model heat loss) within the stated tolerance (ADR 0006)."""
    th = results["parts"]["piston"]["data"]["thermal"]
    assert th["consistent"], th
    assert abs(th["crown_mean_T_K"] - th["T_piston_spec_K"]) <= fea_spec()["piston"]["thermal"]["consistency_tol_K"]
    assert abs(th["heat_W_quarter"]["balance_rel"]) < 1e-2          # steady heat balance of the thermal solve
