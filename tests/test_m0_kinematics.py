"""M0: kinematics match closed forms; Otto cycle matches 1 - r^(1-gamma)."""
import numpy as np
import pytest

from redline.spec import geometry, load
from kin import slider_crank as sc
from thermo.ideal import otto_efficiency, otto_cycle_numeric

G = geometry()
R, L = G.crank_radius, G.conrod_length
TH = np.linspace(0, 4 * np.pi, 7201)


def test_spec_loads_and_displacement():
    assert load()["engine"]["cylinders"] == 4
    assert G.displacement == pytest.approx(1.998e-3, rel=1e-3)


def test_position_endpoints():
    assert sc.displacement_from_tdc(0.0, R, L) == pytest.approx(0.0, abs=1e-15)
    assert sc.displacement_from_tdc(np.pi, R, L) == pytest.approx(G.stroke, rel=1e-12)


def test_volume_ratio_is_compression_ratio():
    vmax = sc.cylinder_volume(np.pi, G)
    vmin = sc.cylinder_volume(0.0, G)
    assert vmax / vmin == pytest.approx(G.compression_ratio, rel=1e-12)


def test_velocity_is_derivative_of_position():
    h = 1e-6
    num = (sc.displacement_from_tdc(TH + h, R, L) - sc.displacement_from_tdc(TH - h, R, L)) / (2 * h)
    np.testing.assert_allclose(sc.piston_velocity(TH, 1.0, R, L), num, rtol=0, atol=1e-9)


def test_acceleration_is_derivative_of_velocity():
    h = 1e-6
    num = (sc.piston_velocity(TH + h, 1.0, R, L) - sc.piston_velocity(TH - h, 1.0, R, L)) / (2 * h)
    np.testing.assert_allclose(sc.piston_acceleration(TH, 1.0, R, L), num, rtol=0, atol=1e-8)


def test_tdc_bdc_acceleration_closed_form():
    w = 7000 * 2 * np.pi / 60
    lam = R / L
    assert sc.piston_acceleration(0.0, w, R, L) == pytest.approx(w**2 * R * (1 + lam), rel=1e-12)
    assert sc.piston_acceleration(np.pi, w, R, L) == pytest.approx(-w**2 * R * (1 - lam), rel=1e-12)


@pytest.mark.parametrize("order,tol", [(2, 2e-2), (4, 1e-3), (6, 1e-4)])
def test_series_converges(order, tol):
    exact = sc.piston_acceleration(TH, 1.0, R, L)
    approx = sc.acceleration_series(TH, 1.0, R, L, order)
    assert np.max(np.abs(exact - approx)) / R < tol


def test_mean_piston_speed():
    w = 6500 * 2 * np.pi / 60
    v = np.abs(sc.piston_velocity(np.linspace(0, 2 * np.pi, 200001)[:-1], w, R, L))
    assert v.mean() == pytest.approx(2 * G.stroke * 6500 / 60, rel=1e-6)


def test_otto_numeric_matches_closed_form():
    eta, imep, _ = otto_cycle_numeric(G)
    assert eta == pytest.approx(otto_efficiency(G.compression_ratio), rel=1e-6)
    assert imep > 0
