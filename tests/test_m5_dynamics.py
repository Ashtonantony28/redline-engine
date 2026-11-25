"""M5 dynamics/NVH tests: analytic anchors first, then the PLAN.md M5 exit criteria.

All results are modelled. Tolerances are stated next to each assertion.
"""
import math

import numpy as np
import pytest

from dyn import balance, bearings, flywheel, forces, params, torsional
from dyn.pressure import PressureTrace, placeholder_wiebe_trace

GRID = 2880


@pytest.fixture(scope="module")
def trace():
    return placeholder_wiebe_trace("full_load")


@pytest.fixture(scope="module")
def model():
    return torsional.build_model()


# ------------------------------------------------------------------ analytic anchors

def test_two_mass_torsional_closed_form():
    J1, J2, k = 0.01, 0.09, 5.0e5
    f_cf = torsional.two_mass_frequency(J1, J2, k)
    f_eig = torsional.eigen_frequencies(np.array([J1, J2]), np.array([k]))[1]
    f_hol = torsional.holzer_frequencies(np.array([J1, J2]), np.array([k]), 2 * f_cf)[1]
    assert f_eig == pytest.approx(f_cf, rel=1e-12)
    assert f_hol == pytest.approx(f_cf, rel=1e-10)


def test_holzer_matches_eigen_on_crank_model(model):
    f_eig = torsional.eigen_frequencies(model.J, model.k)
    f_hol = torsional.holzer_frequencies(model.J, model.k, f_eig[-1] * 1.05)
    assert len(f_hol) == len(f_eig)
    assert np.allclose(f_hol[1:], f_eig[1:], rtol=1e-9)


def test_shaking_primary_and_moments_cancel():
    th = np.linspace(0, params.REV, 1024, endpoint=False)
    w = params.omega_from_rpm(7000)
    s = balance.shaking(th, w)
    g = params.geometry()
    ref = params.masses().reciprocating * g.crank_radius * w**2       # one cylinder's primary amplitude
    tol = 1e-12 * ref                                                    # round-off only
    assert np.abs(s.Fz_primary).max() < tol
    assert np.abs(s.My_primary).max() < tol
    assert np.abs(s.My_secondary).max() < tol
    assert np.abs(s.My).max() < tol and np.abs(s.Mz).max() < tol
    assert np.hypot(s.Fz_rot, s.Fy).max() < tol


def test_mean_torque_identity(trace):
    th = forces.cycle_grid(GRID)
    for rpm in (800, 4500, 7000):
        _, T = forces.multi_cylinder_torque(th, trace, params.omega_from_rpm(rpm))
        # trapezoid on V vs uniform-theta mean: O(h^2) discretisation, 1e-5 relative
        assert T.mean() == pytest.approx(forces.mean_torque_identity(trace), rel=1e-5)


def test_inertia_torque_zero_mean(trace):
    th = forces.cycle_grid(GRID)
    w = params.omega_from_rpm(7000)
    f = forces.cylinder_forces(th, trace, w)
    assert abs(f.T_inertia.mean()) < 1e-9 * np.abs(f.T_inertia).max()


def test_torque_equals_tangential_force_times_radius(trace):
    th = forces.cycle_grid(GRID)
    f = forces.cylinder_forces(th, trace, params.omega_from_rpm(4500))
    r = params.geometry().crank_radius
    assert np.allclose(f.T, f.F_tangential * r, rtol=1e-10, atol=1e-9)


def test_ocvirk_static_short_bearing():
    b = bearings.bearing("main")
    film = bearings.film_for()
    w = params.omega_from_rpm(3000)
    for W in (1.0e3, 1.0e4, 4.0e4):
        eps, att = bearings.static_eccentricity(W, b, w, film)
        # quadrature of the film integrals on the psi grid: 1e-3 relative
        assert bearings.ocvirk_static_load(eps, b, w) == pytest.approx(W, rel=1e-3)
        assert att == pytest.approx(float(bearings.ocvirk_attitude(eps)), abs=1e-3)


def test_flywheel_energy_method_vs_time_integration(trace):
    th = forces.cycle_grid(GRID)
    w = params.omega_from_rpm(4500)
    _, T = forces.multi_cylinder_torque(th, trace, w)
    I = flywheel.rotating_inertia()
    cs = flywheel.cyclic_irregularity(flywheel.fluctuation_energy(th, T), I, w)
    _, ws = flywheel.simulate_speed(th, T, I, w)
    # small-C_s linearisation: agreement to 2 % at C_s ~ 0.005
    assert (ws.max() - ws.min()) / ws.mean() == pytest.approx(cs, rel=0.02)


def test_pressure_interface_validation():
    th = np.linspace(0, 720, 721)
    with pytest.raises(ValueError):
        PressureTrace.from_degrees(th[:-10], np.ones(711) * 1e5)        # does not cover 720 deg
    with pytest.raises(ValueError):
        PressureTrace.from_degrees(th, -np.ones(721))                     # not absolute
    t = PressureTrace.from_degrees(th, 1e5 + th)
    assert t.at(math.radians(730.0)) == pytest.approx(1e5 + 10.0)        # periodic over 720 deg


# ------------------------------------------------------------------ M5 exit criteria

def test_exit_secondary_force_matches_analytic():
    th = np.linspace(0, params.REV, 4096, endpoint=False)
    w = params.omega_from_rpm(7000)
    s = balance.shaking(th, w)
    lam = params.geometry().lam
    exact2 = abs(balance.harmonic(s.Fz, th, 2))
    formula = abs(balance.harmonic(balance.secondary_force_analytic(th, w), th, 2))
    # leading order 4 m r w^2 lam: next term is lam^2/4 relative (4A2 = lam + lam^3/4 + ...)
    assert exact2 / formula - 1 == pytest.approx(lam**2 / 4, abs=lam**4)
    # series to lam^5 (kin.fourier_coefficients) vs exact: truncation O(lam^7), 1e-4
    assert abs(balance.harmonic(s.Fz_secondary, th, 2)) == pytest.approx(exact2, rel=1e-4)
    # the exact force minus the 2nd-order series is bounded by the 4th and higher harmonics
    # (plus the 1e-4 2nd-order truncation): no 1st/3rd/odd content is left
    higher = sum(abs(balance.harmonic(s.Fz, th, q)) for q in range(4, 40, 2))
    assert np.abs(s.Fz - s.Fz_secondary).max() <= higher + 1e-4 * exact2
    for q in (1, 3, 5):
        assert abs(balance.harmonic(s.Fz, th, q)) < 1e-9 * exact2


def test_exit_lanchester_cancels_secondary():
    th = np.linspace(0, params.REV, 4096, endpoint=False)
    w = params.omega_from_rpm(7000)
    s = balance.shaking(th, w)
    pair = balance.lanchester_sizing()
    fz, fy = balance.lanchester_force(th, w, pair)
    exact2 = abs(balance.harmonic(s.Fz, th, 2))
    assert abs(balance.harmonic(s.Fz + fz, th, 2)) < 1e-9 * exact2
    assert np.abs(fy).max() < 1e-9 * exact2


def test_exit_torsional_modes_agree_within_one_percent(model):
    d = torsional.den_hartog_damper(model, params.dyn_spec()["torsional"]["damper"]["ring_inertia"])
    for mdl in (model, torsional.with_damper(model, d)):
        f_eig = torsional.eigen_frequencies(mdl.J, mdl.k)
        f_hol = torsional.holzer_frequencies(mdl.J, mdl.k, f_eig[-1] * 1.05)
        assert np.max(np.abs(f_hol[1:] / f_eig[1:] - 1)) < 0.01


def test_tuned_damper_lowers_resonant_amplitude(model, trace):
    ts = params.dyn_spec()["torsional"]
    d = torsional.den_hartog_damper(model, ts["damper"]["ring_inertia"], ts["damper"]["tuned_mode"])
    md = torsional.with_damper(model, d)
    orders, Cg, Ci = torsional.torque_orders(trace, ts["max_order"])
    f1 = torsional.eigen_frequencies(model.J, model.k)[1]
    q = 6.0                                     # I4 major order with resonance below redline
    assert 60 * f1 / q < params.engine_spec()["engine"]["redline_rpm"]
    i = int(round(2 * q)) - 1
    rpm = np.linspace(60 * f1 / q * 0.6, 60 * f1 / q * 1.3, 200)
    r0 = torsional.forced_response(model, rpm, orders[i:i + 1], Cg[i:i + 1], Ci[i:i + 1])
    r1 = torsional.forced_response(md, rpm, orders[i:i + 1], Cg[i:i + 1], Ci[i:i + 1])
    tw0 = np.abs(r0[0, :, model.hub] - r0[0, :, -1]).max()
    tw1 = np.abs(r1[0, :, md.hub] - r1[0, :, -1]).max()
    assert tw1 < 0.5 * tw0


def test_bearing_loads_consistent(trace):
    th = forces.cycle_grid(720)
    w = params.omega_from_rpm(7000)
    be = bearings.big_end_loads(th, trace, w)
    # frame changes are rotations: magnitudes preserved
    m = np.hypot(*be.F_pin)
    assert np.allclose(np.hypot(*be.pin_frame), m, rtol=1e-12)
    assert np.allclose(np.hypot(*be.rod_frame), m, rtol=1e-12)
    # tangential pin force + big-end centrifugal radial term reproduce forces module
    f = forces.cylinder_forces(th, trace, w)
    assert np.allclose(be.pin_frame[1], f.F_tangential, rtol=1e-10, atol=1e-6)
    # mains carry the sum of all throw forces (pin forces + net throw unbalance)
    ml = bearings.main_bearing_loads(th, trace, w)
    total = np.zeros_like(be.F_pin)
    for ph in params.firing_offsets():
        total += bearings.big_end_loads(th - ph, trace, w).F_pin
    # net unbalance of the 4 flat-plane throws cancels, so the mains sum to the pin forces
    assert np.allclose(ml.sum(0), total, rtol=1e-10, atol=1e-6)


def test_big_end_orbit_film_positive(trace):
    """Smoke test of the mobility orbit (coarse: 2 cycles, 720 steps): film stays inside (0, c)."""
    th = forces.cycle_grid(720)
    w = params.omega_from_rpm(7000)
    be = bearings.big_end_loads(th, trace, w)
    b = bearings.bearing("big_end")
    o = bearings.journal_orbit(th, be.F_pin, w, be.omega_rod, b, w, bearings.film_for(), 2, 720)
    assert 0.0 < o.min_film < b.radial_clearance


# ------------------------------------------------------------------ M10: REAL pressure sources

@pytest.fixture(scope="module")
def m4_full_load():
    from dyn.pressure import m4_trace
    return m4_trace(params.engine_spec()["targets"]["torque_rpm"])


def test_m4_trace_reproduces_m4_imep(m4_full_load):
    """The stitched M4 trace (network open part + two-zone closed part) integrates to the M4 run's own
    net IMEP within 0.5 % (the trace is resampled at 0.25 deg; M4 integrates work on its own grid)."""
    tr = m4_full_load
    assert not tr.meta["placeholder"] and tr.meta["converged"]
    assert forces.imep(tr) / 1e5 == pytest.approx(tr.meta["imep_net_bar"], rel=5e-3)
    assert tr.peak / 1e5 == pytest.approx(tr.meta["p_max_bar_m4"], rel=1e-2)
    # firing peak after TDC, and the mean-torque identity still holds on the real trace
    assert 0.0 < tr.peak_angle_deg < 40.0
    th = forces.cycle_grid(GRID)
    _, T = forces.multi_cylinder_torque(th, tr, params.omega_from_rpm(tr.meta["rpm"]))
    assert T.mean() == pytest.approx(forces.mean_torque_identity(tr), rel=1e-3)


def test_idle_trace_zero_brake_torque():
    """Idle trace: net IMEP equals the Chen-Flynn FMEP (unloaded idle) to 1e-3 bar."""
    from dyn.pressure import idle_trace
    tr = idle_trace()
    assert forces.imep(tr) / 1e5 == pytest.approx(tr.meta["fmep_bar"], abs=1e-3)
    lo, hi = params.dyn_spec()["pressure_m4"]["idle"]["p_intake_bracket"]
    assert float(lo) < tr.meta["p_intake"] < float(hi)


def test_provider_placeholder_ignores_speed():
    from dyn.pressure import PressureProvider
    pv = PressureProvider("placeholder")
    assert pv("full_load", 1500) is pv("full_load", 7000)


def test_bearing_diameters_match_cad():
    """spec/dynamics.yaml journal diameters are the CAD crankshaft's (spec/cad.yaml), so the film model and the
    3D bearing shells describe the same bearings."""
    import yaml
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    dyn = yaml.safe_load((root / "spec" / "dynamics.yaml").read_text())["bearings"]
    crank = yaml.safe_load((root / "spec" / "cad.yaml").read_text())["crankshaft"]
    assert dyn["main"]["diameter"] == pytest.approx(crank["main_diameter"], abs=1e-9)
    assert dyn["big_end"]["diameter"] == pytest.approx(crank["pin_diameter"], abs=1e-9)
