"""M8 exit criteria: MVEM analytic anchors, ECU decoder / scheduler through the ctypes bridge, and
closed-loop HIL behaviour (lambda, idle, tip-in, knock, rev limiter, limp). All results MODELLED.

Verify:  make -C ecu test && pytest -q tests/test_m8_ecu.py && python -m hil.run --target host --cycle tipin

Stated tolerances (never loosened silently):
  manifold filling time constant vs V / (VE Vd N / 120) ........ 1 %   (implicit Euler, dt = 1e-4 s)
  steady throttle flow vs independent orifice formula ........... 1e-12 relative (same equation),
                                                                   steady manifold balance 1e-4
  event scheduler, 7000 rpm, tooth-interpolated ................. 0.5 deg
  arm-ahead (ADR 0008) extrapolation at 12000 rpm/s, analytic .... <= 0.1 deg (a fifth of the 0.5 deg spec)
  steady-state lambda (every combustion event in the window) ..... +/- 3 %
  idle 800 rpm after a 10 Nm accessory step, after settling ...... +/- 25 rpm
  tip-in lambda excursion (2000 rpm, pedal 0.10 -> 0.40 in 0.1 s)  |lambda - 1| <= 0.12 per event
  rev limiter, free-revving WOT ................................. max speed < redline + 250 rpm
"""
from __future__ import annotations

import math
import subprocess
import sys

import numpy as np
import pytest

from mvem import maps as mmaps
from mvem import params
from mvem.plant import Engine
from redline.spec import ROOT

PY_ECU_DIR = ROOT / "ecu"


# --------------------------------------------------------------------------------------------
# fixtures and helpers
# --------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ecu():
    from hil.bridge import HostEcu
    return HostEcu()


@pytest.fixture(scope="module")
def ms():
    return params.mvem_spec()


def independent_orifice(ms, alpha_deg, p_down):
    """Throttle mass flow written out from spec/mvem.yaml without touching mvem.plant."""
    amb, thr = ms["ambient"], ms["throttle"]
    p0, T0, R, g = float(amb["p"]), float(amb["T"]), float(amb["R_air"]), float(amb["gamma_air"])
    psi0 = math.radians(float(thr["closed_plate_angle_deg"]))
    area = (math.pi * float(thr["bore_diameter"]) ** 2 / 4.0
            * (1.0 - math.cos(psi0 + math.radians(alpha_deg)) / math.cos(psi0)) + float(thr["leak_area"]))
    pr = p_down / p0
    pr_crit = (2.0 / (g + 1.0)) ** (g / (g - 1.0))
    if pr <= pr_crit:      # choked: rho* a* A
        flow = math.sqrt(g) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    else:
        flow = math.sqrt(2.0 * g / (g - 1.0) * (pr ** (2.0 / g) - pr ** ((g + 1.0) / g)))
    return float(thr["discharge_coefficient"]) * area * p0 / math.sqrt(R * T0) * flow


# --------------------------------------------------------------------------------------------
# provenance: maps fixture and generated calibration header
# --------------------------------------------------------------------------------------------
def test_generated_calibration_header_is_reproducible():
    r = subprocess.run([sys.executable, str(ROOT / "ecu" / "tools" / "gen_calib.py"), "--check"], cwd=ROOT,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_maps_fixture_matches_a_fresh_thermo_run():
    d = mmaps.load_maps()
    assert d["all_converged"]
    i, j = list(d["rpm"]).index(3000.0), list(d["map"]).index(0.50e5)
    socs = d["soc_deg"][[0, len(d["soc_deg"]) // 2]]
    pts = mmaps.run_point((3000.0, 0.50e5, socs.tolist()))
    for k, p in zip([0, len(d["soc_deg"]) // 2], pts):
        assert p["eta_gross"] == pytest.approx(d["eta_gross"][i, j, k], rel=1e-9)
    ref = int(np.argmin(np.abs(d["soc_deg"] - d["soc_ref_deg"])))
    if ref == len(d["soc_deg"]) // 2:
        assert pts[1]["ve"] == pytest.approx(d["ve"][i, j], rel=1e-9)


# --------------------------------------------------------------------------------------------
# MVEM analytic anchors
# --------------------------------------------------------------------------------------------
def test_throttle_flow_matches_orifice_formula(ms):
    eng = Engine(2000.0, 0.5e5)
    for a in (0.0, 1.0, 5.0, 20.0, 60.0, 84.0):
        for p in (0.2e5, 0.5e5, 0.8e5, 0.98e5):
            assert eng.throttle_flow(a, p) == pytest.approx(independent_orifice(ms, a, p), rel=1e-12)


def test_steady_manifold_balance_is_the_orifice_flow(ms):
    rpm, thr = 2500.0, 6.0
    eng = Engine(rpm, 0.5e5, throttle0_deg=thr)
    for _ in range(3000):
        eng.advance(1.0e-3, throttle_deg=thr, dyno_rpm=rpm)
    p = eng.p_man
    port = eng.port_flow(rpm, p)
    assert port == pytest.approx(independent_orifice(ms, thr, p), rel=1e-4)


def test_manifold_filling_time_constant_closed_form():
    """Choked throttle + constant VE: p relaxes exponentially with tau = V / (VE Vd N / 120)."""
    rpm, ve, dt = 3000.0, 0.80, 1.0e-4
    eng = Engine(rpm, 0.20e5, throttle0_deg=3.0, ve_override=ve)
    thr = 6.0
    p0 = eng.p_man
    n = 20000
    t = np.arange(1, n + 1) * dt
    p = np.empty(n)
    for k in range(n):
        eng.advance(dt, throttle_deg=thr, dyno_rpm=rpm)
        p[k] = eng.p_man
    tau_cf = eng.V_man / (ve * eng.vd * rpm / (params.SEC_PER_MIN * params.REV_PER_CYCLE))
    mdot = eng.throttle_flow(thr, p[-1])
    assert p[-1] < 0.5 * eng.p_amb                                    # stayed choked
    p_inf = mdot * eng.R * eng.T_man * tau_cf / eng.V_man             # closed-form final pressure
    assert p[-1] == pytest.approx(p_inf, rel=1e-6)
    frac = (p - p0) / (p_inf - p0)
    tau_sim = float(np.interp(1.0 - math.exp(-1.0), frac, t))
    assert tau_sim == pytest.approx(tau_cf, rel=0.01)


def test_fuel_film_x_tau_and_injector():
    eng = Engine(1500.0, 0.4e5)
    pw = 4.0e-3
    eng.inject(0, pw)
    m = eng.q_inj * (pw - eng.deadtime())
    assert eng.film[0] == pytest.approx(eng.film_x * m, rel=1e-12)
    assert eng.port[0] == pytest.approx((1.0 - eng.film_x) * m, rel=1e-12)
    tau = eng.film_tau(1500.0)
    eng.advance_flows(tau)
    assert eng.film[0] == pytest.approx(eng.film_x * m * math.exp(-1.0), rel=1e-12)
    assert eng.film[0] + eng.port[0] == pytest.approx(m, rel=1e-12)


def test_crank_inertia_from_spec_and_torque_response():
    from dyn.flywheel import rotating_inertia
    es = params.engine_spec()
    eng = Engine(2000.0, 0.5e5)
    assert eng.J == pytest.approx(rotating_inertia(), rel=1e-12)
    assert eng.J > es["masses"]["crank_inertia"] + es["masses"]["flywheel_inertia"]
    dT, dt = 20.0, 1.0e-3
    eng.T_ind = eng.friction_torque(2000.0) + dT
    w0 = eng.omega
    eng.advance(dt)
    assert eng.omega - w0 == pytest.approx(dT / eng.J * dt, rel=1e-9)


def test_etc_static_balance_against_springs(ms):
    et = {k: float(v) for k, v in ms["etc"].items()}
    eng = Engine(800.0, 0.3e5, throttle0_deg=et["limp_opening_deg"])
    volts = 2.0
    for _ in range(5000):
        eng.etc_step(volts, 1.0e-3)
    F = et["torque_constant"] * volts / et["motor_resistance"]
    lo = et["limp_opening_deg"] + math.degrees((F - et["spring_preload"] - et["static_friction"]) / et["spring_rate"])
    hi = et["limp_opening_deg"] + math.degrees((F - et["spring_preload"] + et["static_friction"]) / et["spring_rate"])
    assert lo - 1e-6 <= eng.thr <= hi + 1e-6        # inside the stiction band of the static balance


def test_lambda_sensor_delay_and_lag(ms):
    ls = ms["lambda_sensor"]
    rpm = 2000.0
    eng = Engine(rpm, 0.5e5)
    eng.lam_queue.append((0.0 + float(ls["transport_delay_s"])
                          + float(ls["transport_delay_cycles"]) * 120.0 / rpm, 1.2))
    delay = float(ls["transport_delay_s"]) + float(ls["transport_delay_cycles"]) * 120.0 / rpm
    dt = 1.0e-4
    ts, ys = [], []
    for _ in range(int((delay + 5 * float(ls["time_constant"])) / dt)):
        eng.t += dt
        eng.lambda_step(dt)
        ts.append(eng.t)
        ys.append(eng.lam_sensor)
    ts, ys = np.array(ts), np.array(ys)
    assert np.all(ys[ts < delay - dt] == 1.0)
    t63 = float(np.interp(1.0 + 0.2 * (1 - math.exp(-1)), ys, ts))
    assert t63 - delay == pytest.approx(float(ls["time_constant"]), abs=2 * dt)


def test_trigger_crossings_per_cycle(ms):
    eng = Engine(1000.0, 0.4e5)
    c = eng.crossings(0.0, 720.0)
    kinds = [k for _, k, _ in c]
    trg = ms["trigger"]
    n_real = int(trg["teeth_total"]) - int(trg["teeth_missing"])
    assert kinds.count("tooth") == 2 * n_real
    assert kinds.count("cam") == 1 and kinds.count("ivc") == 4 and kinds.count("burn") == 4
    teeth = np.array([th for th, k, _ in c if k == "tooth"])
    gaps = np.diff(teeth)
    pitch = 360.0 / int(trg["teeth_total"])
    assert np.sum(np.isclose(gaps, (int(trg["teeth_missing"]) + 1) * pitch)) == 2   # one gap per revolution
    assert np.allclose(gaps[~np.isclose(gaps, 3 * pitch)], pitch)


def test_knock_sensor_hook():
    eng = Engine(3000.0, 0.6e5)
    eng.charge[0].air, eng.charge[0].fuel = 2.0e-4, 2.0e-4 / eng.afr
    eng.charge[0].p_map, eng.charge[0].rpm = 0.6e5, 3000.0
    eng.spark_adv[0] = 20.0
    eng.knock_hook = lambda t, cyl: True
    _, counts = eng.burn(0, 0.0)
    assert counts > eng.kn["injected_counts"] * 0.5


# --------------------------------------------------------------------------------------------
# ECU through the bridge: decoder on synthetic tooth streams, scheduler accuracy
# --------------------------------------------------------------------------------------------
def c_unit_tests_pass():
    r = subprocess.run(["make", "-s", "-C", str(PY_ECU_DIR), "test"], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def test_c_unit_tests():
    code, out = c_unit_tests_pass()
    assert code == 0, out
    assert "0 failures" in out


def tooth_stream(ecu, rpm, th_start, th_end, extra_after=None, noise_every=0, cam=True, tick0=0):
    """Feed exact edges for a constant speed; returns [(theta, status)] after each real tooth."""
    k = ecu.consts
    pitch = 360.0 / k.teeth_total
    tooth0 = k.tooth0_u / k.angle_u_per_deg
    cam_deg = k.cam_edge_u / k.angle_u_per_deg
    om = rpm * 6.0
    edges = []
    kk = math.floor((th_start - tooth0) / pitch) + 1
    while True:
        th = tooth0 + kk * pitch
        if th > th_end:
            break
        if kk % k.teeth_total < k.teeth_total - k.teeth_missing:
            edges.append((th, "t"))
            if extra_after is not None and kk == extra_after:
                edges.append((th + 0.5 * pitch, "x"))
            if noise_every and kk % noise_every == 0:
                edges.append((th + 0.2 * pitch, "x"))
        kk += 1
    if cam:
        j = math.floor((th_start - cam_deg) / 720.0) + 1
        while cam_deg + 720.0 * j <= th_end:
            edges.append((cam_deg + 720.0 * j, "c"))
            j += 1
    edges.sort()
    out = []
    for th, kind in edges:
        tick = tick0 + int(round((th - th_start) / om * k.timer_hz))
        if kind == "c":
            ecu.cam(tick)
        else:
            ecu.crank(tick)
            if kind == "t":
                out.append((th, ecu.status()))
        ecu.drain()
    return out


def decoded_ok(ecu, th, st):
    k = ecu.consts
    raw = k.tooth0_u + st["tooth"] * k.angle_u_per_deg * 360 // k.teeth_total
    ok = raw % (360 * k.angle_u_per_deg) == round((th % 360.0) * k.angle_u_per_deg)
    if st["sync"] == 2:   # same convention as ecu/src/decoder.c: (tooth0 + n pitch + half 360) mod 720
        cyc = (raw + st["half"] * 360 * k.angle_u_per_deg) % (720 * k.angle_u_per_deg)
        ok &= cyc == round((th % 720.0) * k.angle_u_per_deg)
    return ok


def test_decoder_syncs_within_one_rev_after_gap(ecu):
    k = ecu.consts
    pitch = 360.0 / k.teeth_total
    tooth0 = k.tooth0_u / k.angle_u_per_deg
    for start in (3.0, 150.0, 400.0, 610.0):
        ecu.reset()
        trace = tooth_stream(ecu, 2500.0, start, start + 1440.0)
        first = next(th for th, st in trace if st["sync"] >= 1)
        gap_end = tooth0
        while gap_end <= start + 5 * pitch:        # the decoder needs 2 good teeth before a gap
            gap_end += 360.0
        assert first <= gap_end + 1e-9
        assert first - start <= 360.0 + 6 * pitch
        synced = [(th, st) for th, st in trace if st["sync"] >= 1]
        assert all(decoded_ok(ecu, th, st) for th, st in synced)
        assert trace[-1][1]["sync"] == 2 and trace[-1][1]["n_loss"] == 0


def test_decoder_noise_rejection_and_loss_of_sync_recovery(ecu):
    ecu.reset()
    trace = tooth_stream(ecu, 4000.0, 0.0, 2880.0, noise_every=9)
    st = trace[-1][1]
    assert st["n_noise"] > 20 and st["n_loss"] == 0
    assert all(decoded_ok(ecu, th, s) for th, s in trace if s["sync"] >= 1)
    ecu.reset()
    trace = tooth_stream(ecu, 3000.0, 0.0, 2880.0, extra_after=150)
    st = trace[-1][1]
    assert st["n_loss"] == 1 and st["n_sync"] == 2 and st["sync"] == 2
    # the corruption is caught at the next gap; the decoder re-syncs at the first clean gap after it
    # (at most two revolutions) and every tooth decodes correctly from there on
    k = ecu.consts
    corrupt = k.tooth0_u / k.angle_u_per_deg + 150 * 360.0 / k.teeth_total
    resync = next(th for th, s in trace if s["n_sync"] == 2)
    assert resync - corrupt <= 720.0 + 18.0
    assert all(decoded_ok(ecu, th, s) for th, s in trace if th > resync)


def test_scheduler_spark_within_half_degree_at_7000rpm(ecu):
    ecu.reset()
    k = ecu.consts
    adv = 32.0
    for c in range(k.ncyl):
        ecu.sched_set(c, int(adv * k.angle_u_per_deg), 0)
    rpm = 7000.0
    om = rpm * 6.0
    pitch = 360.0 / k.teeth_total
    tooth0 = k.tooth0_u / k.angle_u_per_deg
    cam_deg = k.cam_edge_u / k.angle_u_per_deg
    edges = [(tooth0 + i * pitch, "t") for i in range(0, 60 * 10) if i % 60 < 58]
    edges += [(cam_deg + 720.0 * j, "c") for j in range(3)]
    edges.sort()
    tdc = params.firing_tdc_deg()
    errs = []
    for th, kind in edges:
        tick = int(round(th / om * k.timer_hz))
        (ecu.cam if kind == "c" else ecu.crank)(tick)
        for tk, typ, cyl, _ in ecu.drain():
            if typ == k.ev_spark:
                ang = tk / k.timer_hz * om
                e = (ang - (tdc[cyl] - adv)) % 720.0
                errs.append(min(e, 720.0 - e))
    assert len(errs) >= 3 * k.ncyl
    assert max(errs) < 0.5


ARM_ACCEL_RPM_S = 12000.0        # M8 acceleration case (ecu/tests/test_ecu.c: WOT free-rev to redline)
ARM_EXTRAP_MAX_DEG = 0.1


def test_arm_ahead_window_generated_and_extrapolation_bounded():
    """ADR 0008: the arm-ahead guard is a generated calibration constant (a time, converted to angle
    with the tooth period and capped below one pitch), and the extra extrapolation it adds stays far
    inside the spec. An event d deg past its arming tooth, timed with the period averaged over hist
    deg, is off by alpha d (hist + d) / (2 omega^2) under a constant acceleration alpha (time error of
    a backward-looking mean speed, to second order). The largest d is across the missing-tooth gap: 3
    pitches + the guard, with 4 normal pitches of history; the slowest speed of the M8 case
    (5800 rpm) gives the largest error."""
    import re
    es = params.ecu_spec()
    upd = int(es["hardware"]["angle_units_per_deg"])
    hz = int(es["hardware"]["timer_hz"])
    hdr = (ROOT / "ecu" / "include" / "ecu" / "ecu_calib_gen.h").read_text()
    ticks = int(re.search(r"#define ECU_ARM_AHEAD_TICKS (\d+)u", hdr).group(1))
    max_u = int(re.search(r"#define ECU_ARM_AHEAD_MAX_U (\d+)u", hdr).group(1))
    assert ticks == int(round(float(es["decoder"]["arm_ahead_s"]) * hz)) and ticks > 0
    trg = params.mvem_spec()["trigger"]
    pitch = 360.0 / int(trg["teeth_total"])
    assert max_u == round(pitch * upd) - 1
    alpha = ARM_ACCEL_RPM_S * 6.0                                  # deg/s^2
    hist = int(es["decoder"]["period_average_teeth"]) * pitch
    eps = lambda rpm, d: alpha * d * (hist + d) / (2.0 * (rpm * 6.0) ** 2)
    for rpm in (5800.0, 7000.0):
        g = min(ticks / hz * rpm * 6.0, max_u / upd)
        d_old = (int(trg["teeth_missing"]) + 1) * pitch
        assert eps(rpm, d_old + g) <= ARM_EXTRAP_MAX_DEG
        assert eps(rpm, d_old + g) - eps(rpm, d_old) < 0.01        # what the window adds


# --------------------------------------------------------------------------------------------
# closed loop (C ECU core + MVEM)
# --------------------------------------------------------------------------------------------
from hil.runner import Scenario, run, summary  # noqa: E402

STEADY_POINTS = [(1500.0, 0.10), (2500.0, 0.25), (3500.0, 0.50), (4500.0, 1.00)]


@pytest.mark.parametrize("rpm,pedal", STEADY_POINTS)
def test_closed_loop_steady_lambda_within_3pct(ecu, rpm, pedal):
    res = run(Scenario(duration=4.0, pedal=[[0.0, pedal]], dyno_rpm=[[0.0, rpm]]), ecu)
    lam = res.event_lambdas(2.5, 4.0)
    assert lam.size > 20
    assert np.all(np.abs(lam - 1.0) <= 0.03), (lam.min(), lam.max())
    assert res.log["closed_loop"][-1] == 1.0
    # trims absorb the plant's +4 % injector flow error (spec/mvem.yaml)
    assert -0.08 < res.log["stft"][-1] + res.log["ltft"][-1] < 0.0


def test_idle_holds_800_under_accessory_step(ecu):
    res = run(Scenario.from_cycle("idle"), ecu)
    s = summary(res)
    idle = float(params.engine_spec()["engine"]["idle_rpm"])
    w = res.scenario.windows
    r = res.window("rpm", *w["after_settle"])
    assert np.all(np.abs(r - idle) <= 25.0), (r.min(), r.max())
    assert abs(s["idle_rpm_mean_before"] - idle) <= 25.0
    assert s["idle_rpm_min_after_step"] > 0.5 * idle          # no stall
    assert res.window("idle_active", *w["after_settle"]).min() == 1.0


def test_tipin_lambda_excursion_bounded(ecu):
    res = run(Scenario.from_cycle("tipin"), ecu)
    s = summary(res)
    assert s["steady_lambda_maxdev"] <= 0.03
    assert s["transient_lambda_maxdev"] <= 0.12
    res = run(Scenario.from_cycle("tipout"), ecu)
    assert summary(res)["transient_lambda_maxdev"] <= 0.12


def test_knock_control_retards_and_recovers(ecu):
    res = run(Scenario(duration=8.0, pedal=[[0.0, 0.30]], dyno_rpm=[[0.0, 3000.0]], knock_windows=[[2.0, 2.5]]), ecu)
    kr = res.log["knock_retard_max"]
    t = res.log["t"]
    es = params.ecu_spec()["knock"]
    assert kr[t < 2.0].max() == 0.0                            # base spark is knock free
    assert kr[(t > 2.2) & (t < 2.5)].max() == pytest.approx(float(es["max_retard_deg"]))
    assert kr[-1] == 0.0                                       # fully recovered
    # a lowered knock limit (hotter, worse fuel): the controller finds a retard that stops knock
    res = run(Scenario(duration=6.0, pedal=[[0.0, 0.30]], dyno_rpm=[[0.0, 3000.0]], knock_offset_deg=9.0), ecu)
    kr = res.log["knock_retard_max"][res.log["t"] > 4.0]
    assert 1.0 < kr.mean() < float(es["max_retard_deg"])


def test_rev_limiter_holds_below_redline_plus_margin(ecu):
    es = params.engine_spec()
    res = run(Scenario(duration=3.0, pedal=[[0.0, 0.0], [0.5, 0.0], [0.6, 1.0]], init_rpm=3000.0), ecu)
    redline = float(es["engine"]["redline_rpm"])
    assert res.log["rpm"].max() < redline + 250.0
    assert res.log["rpm"][-1000:].min() > redline - 500.0       # held near the limit, not stalled
    assert res.log["rev_cut"].max() == 1.0


def test_tps_fault_enters_limp_mode(ecu):
    lim = float(params.ecu_spec()["faults"]["limp_rpm"])
    res = run(Scenario(duration=4.0, pedal=[[0.0, 0.1]], init_rpm=1500.0, tps2_fault_at=1.0), ecu)
    t = res.log["t"]
    assert res.log["limp"][t < 1.0].max() == 0.0
    assert res.log["limp"][-1] == 1.0
    assert np.all(res.log["etc_v"][t > 1.1] == 0.0)
    assert res.log["rpm"][t > 1.5].max() < lim + 250.0          # limp rev limit (throttle springs to limp)
    assert res.log["rev_cut"][t > 1.5].max() == 1.0
