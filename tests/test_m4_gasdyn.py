"""M4 exit criteria: 1D finite-volume gas dynamics coupled to four cylinders, full-load curve.

Verify: /home/user/redline-venv/bin/python -m pytest -q tests/test_m4_gasdyn.py tests/test_m2_thermo.py
        /home/user/redline-venv/bin/python -m sim.wot --sweep 1000:7000:250
Analytic references first (CLAUDE.md rule 2), then the engine-level exit criteria of docs/PLAN.md sec. 5.
All engine numbers are MODELLED (no physical engine exists).
"""
import math

import numpy as np
import pytest

from gasdyn import fv
from gasdyn.engine import EngineModel, throttle_area
from gasdyn.exact import (critical_pressure_ratio, orifice_mass_flow, quarter_wave_frequency, riemann_exact,
                          runner_tuned_rpm)
from gasdyn.friction import fmep_chen_flynn
from gasdyn.network import WALL, Network
from gasdyn.spec import gasdyn_spec
from kin.slider_crank import cylinder_volume
from redline.spec import geometry, load
from sim.wot import peaks, run_sweep, sweep_speeds
from thermo.gasexchange import orifice

# Stated tolerances (CLAUDE.md: never loosened silently)
SOD_L1_DENSITY = {100: 8.0e-3, 200: 4.5e-3, 400: 2.5e-3}   # mean |rho - rho_exact| [kg/m^3], rho_left = 1
SOD_MIN_ORDER = 0.7             # observed L1 order between successive refinements (1 is the limit with shocks)
ORIFICE_REL_TOL = 1.0e-3        # steady pipe + orifice mass flow vs closed form
ORIFICE_FN_REL_TOL = 1.0e-12    # orifice law vs the independent closed form
RESONANCE_REL_TOL = 0.02        # organ-pipe fundamental vs c / (4 L)
MASS_ROUNDOFF = 1.0e-12         # relative mass drift of a closed network
PEAK_REL_TOL = 0.10             # WOT torque / power peak vs spec targets (PLAN.md M4)
PEAK_RPM_REL_TOL = 0.12         # speed of the peaks vs spec targets (one coarse 500 rpm step at 4500 rpm)
IMEP_CYCLE_TOL = 1.0e-3         # cycle-to-cycle IMEP change at convergence (0.1 %)
ENGINE_MASS_BALANCE = 5.0e-3    # intake vs exhaust mass per cycle of a converged point
VE_BAND = (0.5, 1.2)            # plausibility band for WOT volumetric efficiency over the sweep
TUNING_STEP_RPM = 250.0         # sweep step of the runner-tuning test = resolution of the VE peak
RUNNER_LENGTHS = (0.45, 0.65)   # m, either side of the spec runner
TUNING_SWEEP = (3000.0, 6500.0, TUNING_STEP_RPM)
JOBS = 4

GAMMA, R_AIR = 1.4, 287.0


def _sod(n):
    net = Network()
    net.add_pipe(1.0, 0.05, GAMMA, R_AIR, n_cells=n)
    net.finalize()
    x = net.x(0)
    left = x < 0.5
    p0 = np.where(left, 1.0e5, 1.0e4)
    rho0 = np.where(left, 1.0, 0.125)
    net.set_pipe(0, p0, p0 / (rho0 * R_AIR))
    t_end = 6.0e-4                      # waves still inside the tube
    net.advance(t_end)
    x, rho, u, p, _ = net.primitives(0)
    rho_e, u_e, p_e = riemann_exact((1.0, 0.0, 1.0e5), (0.125, 0.0, 1.0e4), GAMMA, (x - 0.5) / t_end)
    assert net.nbad == 0
    return float(np.mean(np.abs(rho - rho_e))), float(np.mean(np.abs(p - p_e))) / 1.0e5


# --------------------------------------------------------------------------------------------- (1)
def test_sod_shock_tube_converges_to_exact_riemann_solution():
    err = {n: _sod(n) for n in SOD_L1_DENSITY}
    ns = sorted(err)
    for n in ns:
        assert err[n][0] < SOD_L1_DENSITY[n]
    for a, b in zip(ns[:-1], ns[1:]):
        for k in range(2):
            assert math.log(err[a][k] / err[b][k]) / math.log(b / a) > SOD_MIN_ORDER


# --------------------------------------------------------------------------------------------- (2)
def test_orifice_law_matches_closed_form():
    A, p0, T0 = 3.0e-4, 2.0e5, 320.0
    for pr in (0.999, 0.9, 0.6, critical_pressure_ratio(GAMMA), 0.2):
        ref = orifice_mass_flow(A, p0, T0, pr * p0, GAMMA, R_AIR)
        assert abs(orifice(A, p0, T0, R_AIR, GAMMA, pr * p0, fv.NO_PRESSURE_LINEARISATION) / ref - 1.0) < ORIFICE_FN_REL_TOL
    # choked flow closed form: mdot = A p0 sqrt(g / (R T0)) (2 / (g + 1))^((g + 1) / (2 (g - 1)))
    choked = A * p0 * math.sqrt(GAMMA / (R_AIR * T0)) * (2.0 / (GAMMA + 1.0)) ** ((GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)))
    assert abs(orifice_mass_flow(A, p0, T0, 0.1 * p0, GAMMA, R_AIR) / choked - 1.0) < ORIFICE_FN_REL_TOL


@pytest.mark.parametrize("pr,area_ratio", [(0.8, 0.3), (0.3, 0.3), (0.9, 0.6), (0.95, 1.0)])
def test_steady_isentropic_nozzle_flow(pr, area_ratio):
    """Reservoir -> frictionless pipe -> orifice (or plain open end) -> back pressure: the steady mass flow
    is the isentropic orifice flow from the reservoir stagnation state (subsonic and choked cases)."""
    p0, T0, D = 2.0e5, 300.0, 0.05
    net = Network()
    res = net.add_ambient(p0, T0, GAMMA, R_AIR)
    back = net.add_ambient(pr * p0, T0, GAMMA, R_AIR)
    Ap = Network.area(D)
    net.add_pipe(0.5, D, GAMMA, R_AIR, left=(res, Ap), right=(back, area_ratio * Ap), p=pr * p0, T=T0)
    net.finalize()
    net.advance(0.3)                    # ~200 acoustic transits: steady
    _, rho, u, p, _ = net.primitives(0)
    mdot = rho * u * Ap
    ref = orifice_mass_flow(area_ratio * Ap, p0, T0, pr * p0, GAMMA, R_AIR)
    assert np.all(np.abs(mdot / ref - 1.0) < ORIFICE_REL_TOL)
    assert net.nbad == 0


# --------------------------------------------------------------------------------------------- (3)
def test_open_closed_pipe_resonates_at_quarter_wave_frequency():
    T0, p0, L = 300.0, 1.0e5, 1.0
    net = Network()
    amb = net.add_ambient(p0, T0, GAMMA, R_AIR)
    net.add_pipe(L, 0.05, GAMMA, R_AIR, left=WALL, right=(amb, Network.area(0.05)), n_cells=100, p=p0, T=T0)
    net.finalize()
    x = net.x(0)
    net.set_pipe(0, p0 * (1.0 + 1.0e-3 * np.cos(0.5 * np.pi * x / L)), T0)   # fundamental mode shape
    f_ref = quarter_wave_frequency(math.sqrt(GAMMA * R_AIR * T0), L)
    per_period, periods = 64, 12
    dt = 1.0 / (f_ref * per_period)
    sig = np.empty(per_period * periods)
    for k in range(sig.size):
        net.advance(dt)
        sig[k] = net.primitives(0)[3][0] - p0          # pressure at the closed end
    up = np.where((sig[:-1] < 0.0) & (sig[1:] >= 0.0))[0]
    t_cross = (up + sig[up] / (sig[up] - sig[up + 1])) * dt
    f = 1.0 / np.mean(np.diff(t_cross))
    assert abs(f / f_ref - 1.0) < RESONANCE_REL_TOL
    assert np.max(np.abs(sig[-per_period:])) > 0.5 * np.max(np.abs(sig[:per_period]))   # not over-damped


# --------------------------------------------------------------------------------------------- (4)
def test_closed_network_conserves_mass_to_roundoff():
    """Two plenums at different states, an orifice, an open end, a dead-end pipe, wall friction and heat
    transfer on: the total mass may only change by round-off."""
    net = Network()
    a = net.add_plenum(1.0e-3, 3.0e5, 500.0, GAMMA, R_AIR, Y=1.0)
    b = net.add_plenum(2.0e-3, 1.0e5, 300.0, GAMMA, R_AIR)
    Ap = Network.area(0.03)
    net.add_pipe(0.4, 0.03, GAMMA, R_AIR, left=(a, 0.3 * Ap), right=(b, Ap), friction_multiplier=1.0,
                 heat_multiplier=1.0, wall_T=400.0)
    net.add_pipe(0.4, 0.03, GAMMA, R_AIR, left=(b, Ap), right=WALL, friction_multiplier=1.0, heat_multiplier=1.0)
    net.finalize()
    m0 = net.total_mass()
    drift = 0.0
    for _ in range(20):
        net.advance(5.0e-3)
        drift = max(drift, abs(net.total_mass() / m0 - 1.0))
    pa, pb = net.plenum_state(a)[0], net.plenum_state(b)[0]
    assert net.nstep > 1000 and net.nbad == 0
    assert drift < MASS_ROUNDOFF
    assert abs(pa / pb - 1.0) < 0.02          # the two plenums have equalised (pipe still ringing slightly)


# --------------------------------------------------------------------------------------------- engine
@pytest.fixture(scope="module")
def gspec():
    return gasdyn_spec()


@pytest.fixture(scope="module")
def model():
    return EngineModel()


@pytest.fixture(scope="module")
def wot(gspec):
    """Coarse full-load sweep (same model as `python -m sim.wot`, fewer speeds)."""
    return run_sweep(sweep_speeds(*gspec["wot"]["coarse_sweep_rpm"]), jobs=JOBS)


def test_engine_tables_follow_adr0001_and_kin(model):
    spec, geom = load(), geometry()
    # volume table is kin.slider_crank (minimum at local firing TDC), valve timing is kin.cam at the check lift
    V = model.tabs[0]
    assert np.argmin(V) == 0 and abs(V[0] - geom.clearance_volume) < 1e-15
    assert abs(V[len(V) // 4] - float(cylinder_volume(math.pi, geom))) < 1e-15
    th = model.table_theta_deg
    lift_i, lift_e = model.lift_tables
    check = model.vspec["timing"]["check_lift"]
    val = spec["valves"]
    for lift, opens, closes in ((lift_i, val["ivo_deg"], val["ivc_deg"]), (lift_e, val["evo_deg"], val["evc_deg"])):
        assert abs(np.interp(opens, th, lift) - check) < 1e-7
        assert abs(np.interp(closes, th, lift) - check) < 1e-7
    # seat-to-seat events are longer than the check-lift events
    assert model.ivo_deg < val["ivo_deg"] and model.ivc_deg > val["ivc_deg"]
    assert model.evo_deg < val["evo_deg"] and model.evc_deg > val["evc_deg"]
    # firing order 1-3-4-2, 180 deg apart
    order = [int(c) + 1 for c in np.argsort(model.phase_deg)]
    assert order == spec["engine"]["firing_order"]
    assert sorted(model.phase_deg) == [0.0, 180.0, 360.0, 540.0]


def test_throttle_area_and_friction_correlation(gspec):
    th = gspec["intake"]["throttle"]
    areas = [throttle_area(a, th["bore"], th["closed_angle_deg"]) for a in (th["closed_angle_deg"], 20.0, 45.0, 90.0)]
    assert areas[0] == 0.0 and all(b > a for a, b in zip(areas[:-1], areas[1:]))
    assert abs(areas[-1] - 0.25 * math.pi * th["bore"] ** 2) < 1e-15
    # Chen & Flynn (SAE 650733): FMEP [bar] = 0.137 + 0.005 p_max [bar] + 0.162 Sp [m/s]
    assert abs(fmep_chen_flynn(60.0e5, 10.0, gspec["friction"]) / 1e5 - (0.137 + 0.005 * 60.0 + 0.162 * 10.0)) < 1e-12


def test_converged_point_is_periodic_and_conserves_mass(model):
    op = model.run(4500.0)
    m = op.metrics
    assert op.converged and m["floor_hits"] == 0
    assert op.history[-1]["delta"] < IMEP_CYCLE_TOL and op.history[-2]["delta"] < IMEP_CYCLE_TOL
    for r in op.per_cylinder:
        # per cycle: mass in through the intake valves = mass out through the exhaust valves
        assert abs(r["intake_net"] / r["exhaust_net"] - 1.0) < ENGINE_MASS_BALANCE
        # closed part: the EVO mass handed back equals the trapped mass of the previous cycle
        assert abs(r["m_evo"] / r["m_trapped"] - 1.0) < ENGINE_MASS_BALANCE
    assert m["imep_cyl_spread"] < ENGINE_MASS_BALANCE          # even firing: the four cylinders agree
    assert 0.33 < m["eta_ind_net"] < 0.43                       # same plausibility band as M2 (Heywood ch. 5, 9)
    # part throttle lowers the charge
    part = model.run(4500.0, throttle_deg=25.0)
    model.set_throttle(gasdyn_spec()["intake"]["throttle"]["wot_angle_deg"])
    assert part.converged and part.metrics["ve_trapped"] < 0.9 * m["ve_trapped"]


def test_wot_torque_and_power_peaks_within_10_percent_of_targets(wot):
    tg = load()["targets"]
    pk = peaks(wot)
    assert pk["all_converged"]
    assert all(r["floor_hits"] == 0 for r in wot)
    assert abs(pk["torque_nm"] / tg["torque_nm"] - 1.0) < PEAK_REL_TOL
    assert abs(pk["power_kw"] / tg["power_kw"] - 1.0) < PEAK_REL_TOL
    assert abs(pk["torque_rpm"] / tg["torque_rpm"] - 1.0) < PEAK_RPM_REL_TOL
    assert abs(pk["power_rpm"] / tg["power_rpm"] - 1.0) < PEAK_RPM_REL_TOL
    assert VE_BAND[0] < pk["ve_min"] and pk["ve_max"] < VE_BAND[1]
    # shape: torque falls away on both sides of its peak, power keeps rising past the torque peak
    T = {r["rpm"]: r["torque_nm"] for r in wot}
    assert T[min(T)] < 0.8 * pk["torque_nm"] and T[max(T)] < 0.8 * pk["torque_nm"]
    assert pk["power_rpm"] > pk["torque_rpm"]


def test_runner_length_moves_ve_peak_as_tuning_formulas_predict(model, gspec):
    """Ram tuning: a longer intake runner lowers the speed of the volumetric-efficiency peak. The peak speed
    must lie between the Helmholtz (Engelman) estimate, which scales as L^-1/2, and the quarter-wave estimate
    N = 60 c / (8 L) (pipe fundamental c / 4L at twice the crank frequency), which scales as L^-1, within one
    sweep step, and the ratio of the two peak speeds must lie between the two scalings."""
    run = gspec["intake"]["runner"]
    geom = geometry()
    area = Network.area(run["diameter"])
    c = model.c_intake
    peak = {}
    for L in RUNNER_LENGTHS:
        res = run_sweep(sweep_speeds(*TUNING_SWEEP), {"intake.runner.length": L}, jobs=JOBS)
        assert all(r["converged"] for r in res)
        peak[L] = peaks(res)["ve_peak_rpm"]
        helm = runner_tuned_rpm(c, area, L, geom.swept_volume, geom.compression_ratio)
        quarter = 60.0 * quarter_wave_frequency(c, L) / 2.0
        assert helm - TUNING_STEP_RPM <= peak[L] <= quarter + TUNING_STEP_RPM
    short, long_ = RUNNER_LENGTHS
    assert peak[short] > peak[long_]                                    # direction
    ratio = peak[short] / peak[long_]
    assert (long_ / short) ** 0.5 <= ratio <= (long_ / short) ** 1.0


@pytest.mark.slow
def test_full_sweep_peaks_at_spec_speeds(gspec):
    """The full 250 rpm sweep of `python -m sim.wot`: peaks within 10 % of the targets, at the target speeds
    to within one sweep step."""
    tg = load()["targets"]
    sw = gspec["wot"]["sweep_rpm"]
    pk = peaks(run_sweep(sweep_speeds(*sw), jobs=JOBS))
    assert pk["all_converged"]
    assert abs(pk["torque_nm"] / tg["torque_nm"] - 1.0) < PEAK_REL_TOL
    assert abs(pk["power_kw"] / tg["power_kw"] - 1.0) < PEAK_REL_TOL
    assert abs(pk["torque_rpm"] - tg["torque_rpm"]) <= sw[2]
    assert abs(pk["power_rpm"] - tg["power_rpm"]) <= sw[2]
