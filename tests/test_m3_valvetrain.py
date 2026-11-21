"""M3: cam design, VVT, piston-to-valve clearance, valvetrain dynamics, Hertz stress.

Analytic anchors first (CLAUDE.md rule 2), then the M3 exit criteria. All results modelled.
Verify: /home/user/redline-venv/bin/python -m pytest -q tests/test_m3_valvetrain.py
"""
import dataclasses
import math

import numpy as np
import pytest

from kin import cam as cm
from kin import valve_dynamics as vd
from kin import valvetrain as vt
from redline.spec import geometry

ENGINE, VT = vt.load_specs()
CAMS = cm.design_cams(ENGINE, VT)
LIM = VT["limits"]
FD_STEP = 1e-5          # rad, central-difference step for derivative checks


@pytest.fixture(scope="module")
def results():
    return vt.analyse(ENGINE, VT)


# ============================================================== analytic anchors: polynomials
def test_345_polynomial_coefficients_exact():
    powers, c = cm.poly345_coefficients()
    np.testing.assert_allclose(c, [0, 0, 0, 10, -15, 6], atol=1e-12)


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_profile_boundary_conditions_exact(name):
    p = CAMS[name]
    bm, br, ba = p.main_half_angle, p.ramp_angle, p.ramp_accel_angle
    tol = 1e-12
    # base circle: lift and three derivatives vanish at the start of the event
    d0 = p.derivatives(np.array([-(bm + br), bm + br]))
    np.testing.assert_allclose(d0, 0.0, atol=tol)
    # nose: peak cam lift, zero velocity and jerk (even polynomial)
    dn = p.derivatives(np.array(0.0))
    assert dn[0] == pytest.approx(p.cam_lift, abs=tol)
    assert dn[1] == pytest.approx(0.0, abs=tol) and dn[3] == pytest.approx(0.0, abs=tol)
    # ramp / main junction from both sides: h_r, v_r, 0, 0
    eps = 1e-12
    for a in (bm - eps, bm + eps):
        d = p.derivatives(np.array(a))
        assert d[0] == pytest.approx(p.ramp_height, abs=1e-12)
        assert d[1] == pytest.approx(-p.ramp_velocity, abs=1e-9)
        assert abs(d[2]) < 1e-8 and abs(d[3]) < 1e-6
    # end of the ramp polynomial section: h_a, v_r, 0, 0
    d = p.derivatives(np.array(bm + br - ba))
    assert d[0] == pytest.approx(p.ramp_accel_height, abs=tol)
    assert d[1] == pytest.approx(-p.ramp_velocity, abs=1e-12)
    # valve lift equals the timing check lift exactly at the spec open / close angles
    v = ENGINE["valves"]
    o, c = (v["ivo_deg"], v["ivc_deg"]) if name == "intake" else (v["evo_deg"], v["evc_deg"])
    lift = p.valve_lift_crank(np.radians([o, c]))
    np.testing.assert_allclose(lift, VT["timing"]["check_lift"], atol=1e-10)
    assert p.valve_max_lift == pytest.approx(v[f"{name}_max_lift"], abs=1e-15)


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_profile_derivatives_consistent(name):
    p = CAMS[name]
    phi = np.linspace(-p.half_event * 0.999, p.half_event * 0.999, 4001)
    d = p.derivatives(phi)
    dp, dm = p.derivatives(phi + FD_STEP), p.derivatives(phi - FD_STEP)
    for n in range(3):
        num = (dp[n] - dm[n]) / (2 * FD_STEP)
        scale = np.abs(d[n + 1]).max()
        np.testing.assert_allclose(num, d[n + 1], atol=1e-5 * scale)


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_profile_is_c3_with_bounded_jerk(name):
    """Lift, velocity, acceleration and jerk continuous at every segment join (C2 lift with
    continuous, hence bounded, jerk)."""
    p = CAMS[name]
    joins = [p.main_half_angle, p.main_half_angle + p.ramp_angle - p.ramp_accel_angle, p.half_event]
    eps = 1e-9
    for a in joins:
        for s in (-1.0, 1.0):
            left, right = p.derivatives(np.array(s * (a - eps))), p.derivatives(np.array(s * (a + eps)))
            scale = np.abs(p.derivatives(np.linspace(-p.half_event, p.half_event, 2001))).max(axis=1)
            np.testing.assert_allclose(left, right, atol=1e-6 * scale.max(), rtol=0)
    jerk = p.derivatives(np.linspace(-p.half_event, p.half_event, 20001))[3]
    assert np.isfinite(jerk).all()


# ============================================================== analytic anchors: harmonic cam
HARM = cm.HarmonicCam(lift=0.008, half_width=math.radians(90), base_radius=0.02)  # rho > 0 everywhere


def test_harmonic_cam_closed_form():
    phi = np.linspace(-HARM.half_width, HARM.half_width, 2001)
    d = HARM.derivatives(phi)
    k = math.pi / HARM.half_width
    assert d[0].max() == pytest.approx(HARM.lift, rel=1e-15)
    # simple harmonic motion: y'' = -k^2 (y - h/2) exactly
    np.testing.assert_allclose(d[2], -k**2 * (d[0] - HARM.lift / 2), atol=1e-15)
    num = np.gradient(d[0], phi)
    np.testing.assert_allclose(num[1:-1], d[1][1:-1], atol=1e-6)


@pytest.mark.parametrize("which", ("harmonic", "intake", "exhaust"))
def test_flat_follower_curvature_matches_contour_geometry(which):
    """rho = r_b + y + y'' (formula) vs curvature of the generated cam contour (geometry)."""
    if which == "harmonic":
        rb, half, der = HARM.base_radius, HARM.half_width * 0.95, HARM.derivatives
    else:
        p = CAMS[which]
        rb, half, der = p.base_radius, p.half_event * 0.95, p.derivatives
    phi = np.linspace(-half, half, 40001)
    d = der(phi)
    x, z = cm.flat_follower_contour(rb, phi, d[0], d[1])
    geo = cm.contour_curvature_radius(x, z, phi)
    formula = cm.flat_follower_radius(rb, d[0], d[2])
    np.testing.assert_allclose(geo[5:-5], formula[5:-5], rtol=1e-4)


# ============================================================== analytic anchors: dynamics
MODEL = vd.model_from_spec("intake", VT)


def test_single_dof_period_matches_sqrt_k_over_m():
    m, k = MODEL.single_mass, MODEL.spring_rate
    t = vd.free_vibration_period(m, k, x0=1e-3, steps_per_period=400, periods=20)
    assert t == pytest.approx(2 * math.pi * math.sqrt(m / k), rel=1e-4)


def test_multi_mass_reduces_to_single_dof():
    """Massless spring, rigid contact/stem: first mode -> sqrt(k / (m_t + m_v)) / 2 pi."""
    rigid = MODEL.k_stem * 1e4
    m = dataclasses.replace(MODEL, spring_mass=MODEL.spring_mass * 1e-6, k_contact=0.0, k_stem=rigid)
    f1 = m.natural_frequencies(contact=False, seat=False)[0]
    expected = math.sqrt(MODEL.spring_rate / (MODEL.m_tappet + MODEL.m_valve)) / (2 * math.pi)
    assert f1 == pytest.approx(expected, rel=1e-4)


def test_spring_chain_matches_discrete_and_continuous_surge():
    fd = vd.surge_frequencies_discrete(MODEL.spring_rate, MODEL.spring_mass, MODEL.n_spring)
    np.testing.assert_allclose(MODEL.spring_chain_frequencies()[:3], fd, rtol=1e-10)
    fine = dataclasses.replace(MODEL, n_spring=60)
    fc = vd.surge_frequencies_continuous(MODEL.spring_rate, MODEL.spring_mass)
    np.testing.assert_allclose(fine.spring_chain_frequencies()[:3], fc, rtol=2e-3)


def test_static_spring_force_is_hooke():
    for x in (0.0, 0.002, 0.0105):
        assert MODEL.static_valve_force(x) == pytest.approx(
            MODEL.preload + MODEL.spring_rate * x, rel=1e-12)
        assert vd.spring_force(x, MODEL.preload, MODEL.spring_rate) == MODEL.preload + MODEL.spring_rate * x


def test_dynamic_model_tracks_quasi_static_at_low_speed():
    p = CAMS["intake"]
    rpm = ENGINE["engine"]["idle_rpm"]
    r = vd.simulate(MODEL, p, [rpm], revolutions=1, step_safety=VT["simulation"]["step_safety"],
                    record_points=1441)
    mid = np.abs(r.phi) < 0.8 * p.main_half_angle
    qs = vd.quasi_static_contact_force(MODEL, p, r.phi[mid], rpm)
    dyn = r.contact[mid, 0]
    assert np.abs(dyn - qs).max() < 0.03 * MODEL.preload


def test_hertz_line_contact_identities():
    e, nu = 210e9, 0.3
    assert vt.reduced_modulus(e, nu, e, nu) == pytest.approx(e / (2 * (1 - nu**2)), rel=1e-15)
    f, r, w = 600.0, 0.008, 0.014
    es = vt.reduced_modulus(170e9, 0.27, e, nu)
    p0 = vt.hertz_line_pressure(f, r, w, es)
    b = vt.hertz_contact_halfwidth(f, r, w, es)
    assert p0 == pytest.approx(2 * f / w / (math.pi * b), rel=1e-12)   # p0 = 2 P' / (pi b)


# ============================================================== exit criteria
@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_cam_shape_and_tappet_geometry(name, results):
    r = results[name]
    p = r.profile
    phi = np.linspace(-p.half_event, 0.0, 20001)
    assert (p.derivatives(phi)[1] >= -1e-15).all()                    # monotone rise
    assert r.min_rho > LIM["min_radius_of_curvature"] > 0              # no undercut
    assert r.required_tappet_diameter <= VT["tappet"]["diameter"]      # contact stays on face
    assert r.coil_bind_margin >= VT["spring"]["min_coil_bind_margin"]


def test_vvt_phasing_shifts_intake_events():
    p = CAMS["intake"]
    theta = np.linspace(0, 4 * math.pi, 7201)
    adv = vt.vvt_advances(ENGINE, VT)
    assert math.degrees(adv[-1]) == pytest.approx(ENGINE["valves"]["vvt_intake_range_deg"])
    for a in adv[[0, len(adv) // 2, -1]]:
        np.testing.assert_allclose(p.valve_lift_crank(theta, a), p.valve_lift_crank(theta + a), atol=1e-15)


def test_chamber_geometry_is_consistent_with_compression_ratio():
    g = geometry(ENGINE)
    assert VT["combustion_chamber"]["valve_to_crown_at_tdc"] < g.clearance_volume / g.piston_area


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_piston_to_valve_clearance_over_full_vvt_range(name, results):
    r = results[name]
    if name == "intake":
        assert math.degrees(r.advances.max()) == pytest.approx(ENGINE["valves"]["vvt_intake_range_deg"])
    assert r.min_clearance > LIM["min_piston_valve_clearance"]


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_no_float_below_limit(name, results):
    r = results[name]
    assert r.qs_float_rpm > LIM["float_free_rpm"]
    assert r.float_rpm > LIM["float_free_rpm"]
    below = r.sweep.rpm <= LIM["float_free_rpm"]
    assert not r.sweep.floats[below].any()


@pytest.mark.parametrize("name", ("intake", "exhaust"))
def test_hertz_stress_under_material_limit(name, results):
    assert results[name].peak_hertz < LIM["hertz_limit"]


def test_report_and_outputs(results, tmp_path):
    text = vt.report(results, ENGINE, VT)
    assert "MODELLED" in text
    vt.write_outputs(results, ENGINE, VT, tmp_path)
    for f in ("cam_intake.csv", "cam_exhaust.csv", "sweep_intake.csv", "sweep_exhaust.csv",
              "clearance_intake_vvt.csv", "cam_profiles.png", "vvt_clearance.png", "float_sweep.png"):
        assert (tmp_path / f).stat().st_size > 0
