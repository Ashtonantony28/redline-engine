"""M10 stretch: turbocharged variant Redline I4-2.0T (synthetic maps, wastegate PI, intercooler, KLSA).

Verify: /home/user/redline-venv/bin/python -m pytest -q tests/test_turbo.py tests/test_m4_gasdyn.py
        /home/user/redline-venv/bin/python -m sim.wot_turbo --sweep 1500:6500:500
Analytic anchors first (CLAUDE.md rule 2): compressor power closed form, turbine isentropic relation and
nozzle flow, shaft power balance, map surge / choke monotonicity, NA engine unchanged by the hooks; then the
engine-level checks on two steady WOT points and a tip-in. All numbers are MODELLED (no physical engine);
the compressor and turbine maps are SYNTHETIC (turbo/maps.py).
"""
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pytest

from gasdyn import fv
from gasdyn.elements import Intercooler, plenum_inject, plenum_remove_energy
from gasdyn.engine import EngineModel
from gasdyn.exact import critical_pressure_ratio, orifice_mass_flow
from gasdyn.network import Network
from redline.spec import ROOT, load
from turbo.control import BoostPI, TurboShaft, rpm_to_rad
from turbo.maps import CompressorMap, TurbineMap
from turbo.spec import engine_spec, turbo_spec

# Stated tolerances (CLAUDE.md: never loosened silently)
CLOSED_FORM_REL = 1.0e-12      # map relations vs their closed forms (same arithmetic up to round-off)
ENERGY_REL = 1.0e-12           # rotor kinetic-energy update vs net power x dt
FIXTURE_REL = 1.0e-5           # NA engine vs the committed M4 fixture (written with 6 significant digits)
SURGE_SLOPE_REL = 1.0e-3       # |dPR/dm| at the surge point relative to the mean line slope
LAG_BAND_S = (0.2, 3.0)        # plausibility band of the 2000 rpm time to 90 % boost rise [s]
NA_FIXTURE_RPM = 4500.0

TS = turbo_spec()
GAMMA_AIR, R_AIR = 1.4, 287.0
GAMMA_EX, R_EX = 1.33, 290.0
P1, T1 = 1.0e5, 298.15


@pytest.fixture(scope="module")
def cmap():
    return CompressorMap(TS["compressor"], GAMMA_AIR, R_AIR)


@pytest.fixture(scope="module")
def tmap():
    return TurbineMap(TS["turbine"], GAMMA_EX, R_EX)


# ---------------------------------------------------------------------------------------------
# analytic anchors
# ---------------------------------------------------------------------------------------------
def test_compressor_power_matches_isentropic_closed_form(cmap):
    """P_C = m cp T1 (PR^((g-1)/g) - 1) / eta over the whole map (Euler work = isentropic work / eta)."""
    cp = GAMMA_AIR * R_AIR / (GAMMA_AIR - 1.0)
    k = (GAMMA_AIR - 1.0) / GAMMA_AIR
    n = 0
    for krpm in (50.0, 80.0, 110.0, 140.0, 170.0):
        w = rpm_to_rad(krpm * 1.0e3)
        m_s, _, m_ch = cmap.limits(w, P1, T1)
        for m in np.linspace(m_s, 0.98 * m_ch, 7):
            pt = cmap.line(w, P1, T1, m)
            P = cmap.power(w, P1, T1, m)
            closed = m * cp * T1 * (pt["pr"] ** k - 1.0) / pt["eta"]
            assert P == pytest.approx(closed, rel=CLOSED_FORM_REL)
            assert P == pytest.approx(m * pt["lam"] * pt["U2"] ** 2, rel=CLOSED_FORM_REL)
            n += 1
    assert n == 35


def test_compressor_outlet_temperature_from_power(cmap):
    w = rpm_to_rad(120.0e3)
    m_s, _, m_ch = cmap.limits(w, P1, T1)
    m_target = 0.5 * (m_s + m_ch)
    pr_target = cmap.line(w, P1, T1, m_target)["pr"]
    op = cmap.operating_point(w, P1, T1, lambda m: pr_target)
    assert op["mdot"] == pytest.approx(m_target, rel=1.0e-8)
    cp = GAMMA_AIR * R_AIR / (GAMMA_AIR - 1.0)
    T2s = T1 * op["pr"] ** ((GAMMA_AIR - 1.0) / GAMMA_AIR)
    assert op["T2"] == pytest.approx(T1 + (T2s - T1) / op["eta"], rel=1.0e-10)
    assert op["T2"] == pytest.approx(T1 + op["power"] / (op["mdot"] * cp), rel=CLOSED_FORM_REL)


def test_turbine_isentropic_relation_and_nozzle_flow(tmap):
    cp = GAMMA_EX * R_EX / (GAMMA_EX - 1.0)
    k = (GAMMA_EX - 1.0) / GAMMA_EX
    T03, p4 = 1100.0, 1.05e5
    for p03 in (1.2e5, 1.6e5, 2.2e5, 3.0e5):
        m = tmap.flow(p03, T03, p4)
        assert m == pytest.approx(orifice_mass_flow(tmap.A, p03, T03, p4, GAMMA_EX, R_EX), rel=CLOSED_FORM_REL)
        c_is = math.sqrt(2.0 * cp * T03 * (1.0 - (p4 / p03) ** k))
        w_opt = tmap.bsr_opt * c_is / (tmap.D / 2.0)
        r = tmap.power(m, p03, T03, p4, w_opt)
        assert r["bsr"] == pytest.approx(tmap.bsr_opt, rel=1.0e-12)
        assert r["eta"] == pytest.approx(tmap.eta_max, rel=1.0e-12)
        assert r["power"] == pytest.approx(m * cp * T03 * r["eta"] * (1.0 - (p4 / p03) ** k), rel=CLOSED_FORM_REL)
        # off-design blade-speed ratio always costs efficiency
        assert tmap.power(m, p03, T03, p4, 0.7 * w_opt)["eta"] < r["eta"]
        assert tmap.power(m, p03, T03, p4, 1.3 * w_opt)["eta"] < r["eta"]
    # choked nozzle: flow independent of the back pressure below the critical ratio
    crit = critical_pressure_ratio(GAMMA_EX)
    assert tmap.flow(3.0e5, T03, 0.9 * crit * 3.0e5) == pytest.approx(tmap.flow(3.0e5, T03, 0.5 * crit * 3.0e5),
                                                                       rel=CLOSED_FORM_REL)


def test_steady_shaft_power_balance_and_energy_update():
    sh = TurboShaft(TS["shaft"])
    w0 = rpm_to_rad(100.0e3)
    sh.omega = w0
    P_C = 8.0e3
    P_T = P_C / sh.eta_m                     # steady balance: eta_m P_T = P_C
    for _ in range(1000):
        sh.step(P_T, P_C, 1.0e-4)
    assert sh.omega == pytest.approx(w0, rel=ENERGY_REL)
    dt, surplus = 1.0e-3, 500.0
    E0 = 0.5 * sh.J * sh.omega ** 2
    sh.step((P_C + surplus) / sh.eta_m, P_C, dt)
    assert 0.5 * sh.J * sh.omega ** 2 - E0 == pytest.approx(surplus * dt, rel=1.0e-9)


def test_compressor_map_surge_and_choke_lines_are_monotonic(cmap):
    speeds = [rpm_to_rad(k * 1.0e3) for k in range(40, 191, 10)]
    lim = np.array([cmap.limits(w, P1, T1) for w in speeds])
    m_s, pr_s, m_ch = lim[:, 0], lim[:, 1], lim[:, 2]
    assert np.all(np.diff(m_s) > 0.0) and np.all(np.diff(pr_s) > 0.0) and np.all(np.diff(m_ch) > 0.0)
    assert np.all(m_ch > m_s)
    for w in speeds:
        m, pr, eta = cmap.speed_line(w, P1, T1, n=200)
        assert np.all(np.diff(pr) < 0.0)                     # right of surge: PR falls with flow
        assert np.all(eta <= TS["compressor"]["eta_peak"] + 1e-12) and np.all(eta >= 0.0)
        # surge = zero slope: one-sided slope at the surge point is ~0 relative to the mean slope
        ms, prs, _ = cmap.limits(w, P1, T1)
        h = 1.0e-5 * ms
        slope = (cmap.line(w, P1, T1, ms + h)["pr"] - prs) / h
        assert abs(slope) < SURGE_SLOPE_REL * abs((pr[-1] - pr[0]) / (m[-1] - m[0]))
    # demand above the surge pressure ratio is flagged and lowers the flow
    w = speeds[8]
    ms, prs, _ = cmap.limits(w, P1, T1)
    op = cmap.operating_point(w, P1, T1, lambda m: prs + 0.5 * TS["compressor"]["surge_band_pr"])
    assert op["surge"] and op["mdot"] < ms
    # passive restriction beyond choke and at a stopped wheel: PR < 1 gives positive flow
    op = cmap.operating_point(rpm_to_rad(TS["shaft"]["min_speed_rpm"]), P1, T1, lambda m: 0.95)
    assert op["mdot"] > 0.0 and op["pr"] == pytest.approx(0.95, rel=1.0e-8)


def test_intercooler_and_plenum_source_terms():
    ic = Intercooler(0.75, 300.0, 1.0e4, 0.2, 2.0)
    assert ic.T_out(400.0) == pytest.approx(325.0, rel=CLOSED_FORM_REL)
    assert ic.dp(0.1, 1.0) == pytest.approx(1.0e4 * 0.25 * 2.0, rel=CLOSED_FORM_REL)
    net = Network()
    iv = net.add_plenum(1.0e-3, 1.0e5, 300.0, GAMMA_AIR, R_AIR)
    net.add_pipe(0.1, 0.01, GAMMA_AIR, R_AIR, left=(iv, 0.0), n_cells=4)
    net.finalize()
    m0, E0 = net.S[iv, fv.S_M], net.S[iv, fv.S_E]
    cp = GAMMA_AIR * R_AIR / (GAMMA_AIR - 1.0)
    plenum_inject(net, iv, 1.0e-4, 330.0)
    assert net.S[iv, fv.S_M] - m0 == pytest.approx(1.0e-4, rel=CLOSED_FORM_REL)
    assert net.S[iv, fv.S_E] - E0 == pytest.approx(1.0e-4 * cp * 330.0, rel=1.0e-9)
    plenum_remove_energy(net, iv, 1.0)
    assert net.S[iv, fv.S_E] - E0 == pytest.approx(1.0e-4 * cp * 330.0 - 1.0, rel=1.0e-9)


def test_boost_pi_holds_and_saturates():
    pi = BoostPI(TS["boost_control"], TS["wastegate"])
    pi.target = 1.6e5
    for _ in range(2000):
        op = pi.step(1.2e5, 1.0e-3)              # far below target: wastegate driven shut
    assert op < 1.0e-3 and pi.i == 0.0
    for _ in range(2000):
        op = pi.step(2.0e5, 1.0e-3)              # far above target: driven open, no windup beyond 1
    e = (pi.target - 2.0e5) / pi.target
    assert op > 1.0 - 1.0e-3 and pi.i <= 1.0
    assert 1.0 <= pi.i - pi.kp * e <= 1.0 + pi.ki * abs(e) * 1.0e-3 + 1.0e-12   # frozen once saturated
    pi.target = None
    nd = TS["boost_control"]["no_demand_opening"]
    for _ in range(10):
        op = pi.step(1.0e5, 1.0)
    assert op == pytest.approx(nd, abs=1.0e-6) and pi.i == nd       # no boost demand: spring position


def test_variant_is_an_override_layer():
    base = load()
    cr = base["engine"]["compression_ratio"]
    sp = engine_spec(TS)
    assert sp["engine"]["compression_ratio"] == TS["variant"]["engine_overrides"]["engine.compression_ratio"]
    assert sp["engine"]["name"] == TS["variant"]["name"]
    assert load()["engine"]["compression_ratio"] == cr               # engine.yaml untouched


def test_na_engine_with_noop_hooks_reproduces_committed_m4_point():
    """Turbo hooks off (boost=None) plus no-op sub-step / cycle hooks: the M4 fixture point is reproduced."""
    import csv
    with open(ROOT / "fixtures" / "wot" / "wot_coarse.csv") as f:
        ref = next(r for r in csv.DictReader(f) if float(r["rpm"]) == NA_FIXTURE_RPM)
    eng = EngineModel()
    assert eng.boost is None and eng.v_charge is None and eng.turbine_pipe is None
    calls = []
    eng.substep_fn = lambda net, th, dt: calls.append(dt)
    eng.substep_deg = TS["steady"]["substep_deg"]
    eng.cycle_fn = lambda icyc, rec: None
    m = eng.run(NA_FIXTURE_RPM).metrics
    assert calls
    for k in ("torque_nm", "power_kw", "ve_trapped", "bsfc_g_per_kwh", "p_max_bar"):
        assert m[k] == pytest.approx(float(ref[k]), rel=FIXTURE_REL), k
    assert m["cycles_run"] == int(ref["cycles_run"])


# ---------------------------------------------------------------------------------------------
# engine level: steady WOT points (knock-limited spark) and the 2000 rpm tip-in
# ---------------------------------------------------------------------------------------------
def _steady(rpm):
    from sim.wot_turbo import _point
    return _point(rpm)


def _tipin(_):
    """Tip-in with the spark from TurboEngine.wot_spark_estimate (first knock-limited pass of the steady
    point) so that it runs in parallel with the steady points; `python -m sim.wot_turbo` uses the final
    steady spark instead."""
    from sim.wot_turbo import _engine, run_tipin
    tp = run_tipin(_engine().wot_spark_estimate(float(TS["tipin"]["rpm"])))
    tp.pop("log")
    return tp


@pytest.fixture(scope="module")
def engine_runs():
    from sim.wot import _warm_kernels
    _warm_kernels()
    rpms = [float(r) for r in TS["sweep"]["test_rpm"]]
    with ProcessPoolExecutor(max_workers=len(rpms) + 1) as pool:
        tip = pool.submit(_tipin, None)
        futs = {r: pool.submit(_steady, r) for r in rpms}
        steady = {r: f.result() for r, f in futs.items()}
        tip = tip.result()
    return steady, tip


def test_steady_points_converge_and_hold_boost_target(engine_runs):
    steady, _ = engine_runs
    tol = float(TS["boost_control"]["hold_tol"])
    for rpm, r in steady.items():
        assert r["converged"], rpm
        assert r["floor_hits"] == 0
        # either the target is held within +/-3 %, or the wastegate is shut and the turbo cannot reach it
        if abs(r["boost_error"]) > tol:
            assert r["boost_error"] < 0.0 and r["wg_opening"] < 1.0e-3, rpm
    held = [r for r in steady.values() if abs(r["boost_error"]) <= tol]
    assert held, "boost target held at no tested speed"
    for r in steady.values():
        assert abs(r["shaft_balance"]) < 0.05       # cycle-mean eta_m P_T = P_C at the steady state


def test_compressor_stays_right_of_surge_with_margin(engine_runs):
    steady, _ = engine_runs
    sm = float(TS["steady"]["surge_margin_min"])
    for rpm, r in steady.items():
        assert r["surge_margin"] >= sm, (rpm, r["surge_margin"])
        assert r["surge_time"] == 0.0, rpm
        assert r["mdot_c"] < r["m_surge"] * 10.0 and not r["overspeed"]


def test_knock_limited_spark_is_respected(engine_runs):
    steady, _ = engine_runs
    from thermo.m7spec import combustion_spec
    k = combustion_spec()["knock"]
    for rpm, r in steady.items():
        assert r["spark_deg"] <= r["mbt_deg"] + 1.0e-9
        if r["knock_limited"]:
            assert r["spark_deg"] == pytest.approx(r["klsa_deg"] - float(k["klsa_margin_deg"]), abs=1.0e-9)
        assert r["knock_margin"] < float(k["integral_threshold"]), rpm   # no knock at the operating spark


def test_tipin_turbo_lag_is_reported(engine_runs):
    _, tip = engine_runs
    assert tip["map_pre"] < 0.9 * tip["target"]
    assert np.isfinite(tip["t90_rise"]), "boost never reached 90 % of the rise"
    assert LAG_BAND_S[0] <= tip["t90_rise"] <= LAG_BAND_S[1]
    assert tip["torque_final"] > tip["torque_pre"]
    # cycle-mean compressor flow stays right of surge during the spool-up (short instantaneous excursions
    # within the exhaust-pulse period are reported as surge_time, not failed: deep surge is not modelled)
    assert tip["surge_margin_cycle_min"] > 0.0
