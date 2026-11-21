"""M10 tests: vehicle model, WLTP, acceleration and the three.js replay data.

Analytic anchors first (constant-speed road-load power, WLTC trace distance/duration, coast-down),
then the plausibility bands of spec/vehicle.yaml (sources there) and the replay data checks.
All results are MODELLED. Tolerances are stated next to each assertion.
"""
import json
import shutil
import subprocess

import numpy as np
import pytest

from kin.slider_crank import cylinder_volume, piston_position
from redline.spec import ROOT, geometry, load
from vehicle import params as vp
from vehicle.engine_map import engine_map
from vehicle.model import Powertrain, RoadLoad, coastdown_numeric, coastdown_time_analytic, simulate_cycle
from vehicle.wltc import distance, load_trace, phase_distances

VS = vp.vehicle_spec()
KMH = vp.KMH_PER_MS


# ------------------------------------------------------------------ analytic anchors

def test_road_load_power_closed_form():
    """Constant-speed road-load power equals m g c_r v + 1/2 rho c_d A v^3 (flat road, c_r2 = 0)
    to 1e-12 relative, and the backward simulation of a constant-speed trace books exactly that power
    at the wheel and P / eta (+ accessory) at the engine."""
    v = VS["vehicle"]
    m = vp.wltp_test_mass()
    road = RoadLoad.from_spec(m)
    for kmh in (30.0, 50.0, 100.0, 130.0):
        u = kmh / KMH
        closed = m * v["gravity"] * v["rolling_resistance"] * u + 0.5 * v["air_density"] * v["drag_coefficient"] \
            * v["frontal_area"] * u**3
        assert float(road.power(u)) == pytest.approx(closed, rel=1e-12)
    u = 100.0 / KMH
    t = np.arange(0.0, 11.0)
    res = simulate_cycle(t, np.full_like(t, u), Powertrain(), m)
    assert np.allclose(res.p_wheel, float(road.power(u)), rtol=1e-12)
    pt = Powertrain()
    p_eng = res.torque * res.rpm * vp.RAD_S_PER_RPM
    assert np.allclose(p_eng, float(road.power(u)) / pt.eta + pt.p_acc, rtol=1e-12)


def test_wltc_trace_distance_and_duration():
    """GTR 15 Class 3b: 1800 s at 1 Hz, 23.27 km (+-10 m), phase distances to 1 m, v_max 131.3 km/h."""
    sim = VS["simulation"]
    ref = sim["wltc_reference"]
    t, v = load_trace()
    assert len(t) == int(ref["duration_s"]) + 1 and t[-1] - t[0] == ref["duration_s"]
    assert np.all(np.diff(t) == 1.0)
    assert distance(t, v) / 1e3 == pytest.approx(ref["distance_km"], abs=ref["distance_tol_km"])
    assert v.max() * KMH == pytest.approx(ref["v_max_kmh"], abs=1e-9)
    for name, d in phase_distances(t, v).items():
        assert d / 1e3 == pytest.approx(sim["wltc_phases"][name]["distance_km"], abs=sim["wltc_phase_distance_tol_km"])


def test_coastdown_matches_analytic():
    """Neutral coast-down 120 -> 20 km/h under m dv/dt = -(A + C v^2): RK4 (dt 0.01 s) matches the
    closed form t = m / sqrt(A C) [atan(v0 sqrt(C/A)) - atan(v1 sqrt(C/A))] to 1e-6 relative."""
    m = vp.wltp_test_mass()
    road = RoadLoad.from_spec(m)
    A, C = road.coast_coefficients
    m_eq = vp.equivalent_mass(m, None)
    v0, v1 = VS["simulation"]["coastdown_from_kmh"] / KMH, VS["simulation"]["coastdown_to_kmh"] / KMH
    t_an = coastdown_time_analytic(A, C, m_eq, v0, v1)
    assert coastdown_numeric(road, m_eq, v0, v1, 0.01) == pytest.approx(t_an, rel=1e-6)
    # initial deceleration equals the road-load force over the equivalent mass
    assert (A + C * v0**2) / m_eq == pytest.approx(float(road.force(v0)) / m_eq, rel=1e-12)
    # sanity: a compact sedan coasts 120 -> 20 km/h in roughly a minute or more (no engine drag)
    assert 40.0 < t_an < 200.0


def test_equivalent_mass_reflects_engine_inertia():
    """m_eq(gear) - m_eq(neutral) = J_engine (i_g i_f)^2 / r^2 exactly."""
    m = vp.wltp_test_mass()
    r = VS["vehicle"]["wheel_radius"]
    for g in range(1, len(VS["driveline"]["gear_ratios"]) + 1):
        d = vp.equivalent_mass(m, g) - vp.equivalent_mass(m, None)
        assert d == pytest.approx(vp.engine_inertia() * vp.overall_ratio(g) ** 2 / r**2, rel=1e-12)


def test_willans_fit_reproduces_m7_maps():
    """The affine fuel-MEP (Willans) fit reproduces every throttled M7 map point within 3 % (rms < 1.5 %)
    and the knock-limited WOT points it was not fitted to within 6 %; the motoring torque is negative."""
    em = engine_map()
    rep = em.fit_report()
    assert rep["throttled_max_abs_rel_err"] < 0.03 and rep["throttled_rms_rel_err"] < 0.015
    assert rep["wot_max_abs_rel_err"] < 0.06
    rpm = np.linspace(em.idle_rpm, em.redline_rpm, 25)
    assert np.all(em.motoring_torque(rpm) < 0.0)
    assert np.all(em.fuel_rate(rpm, em.motoring_torque(rpm)) == pytest.approx(0.0, abs=1e-15))


def test_wot_envelope_matches_fixture_targets():
    """The WOT fixture reproduces the M4 headline (200 Nm at 4500 rpm, ~109 kW at 6500 rpm, within 10 %)."""
    em = engine_map()
    t = load()["targets"]
    assert float(em.max_torque(t["torque_rpm"])) == pytest.approx(t["torque_nm"], rel=0.10)
    assert float(em.max_power(t["power_rpm"])) / 1e3 == pytest.approx(t["power_kw"], rel=0.10)
    assert float(em.max_torque(em.redline_rpm + 1.0)) == 0.0


# ------------------------------------------------------------------ plausibility (MODELLED results)

@pytest.fixture(scope="module")
def wltp():
    from vehicle.wltp import run
    return run()


def test_wltp_fuel_consumption_plausible(wltp):
    """WLTP combined in the band for 2.0 L NA compact petrol manuals (spec/vehicle.yaml plausibility,
    6.0-8.5 L/100 km), trace followed with no torque shortfall, fuel cut books zero fuel."""
    s, res = wltp
    lo, hi = VS["plausibility"]["wltp_l_per_100km"]
    assert lo <= s["l_per_100km"] <= hi
    assert s["seconds_with_shortfall"] == 0
    assert s["distance_km"] == pytest.approx(VS["simulation"]["wltc_reference"]["distance_km"], abs=0.01)
    cut = np.array([m == "fuel_cut" for m in res.mode])
    assert cut.any() and np.all(res.fuel_rate[cut] == 0.0)
    # CO2 follows from carbon balance: g/km = fuel g/km x 44.01 x 8 / 114.23 (iso-octane), to 1e-9
    f = VS["fuel"]
    ratio = f["carbon_atoms"] * (f["atomic_mass_C"] + 2 * f["atomic_mass_O"]) / (
        f["carbon_atoms"] * f["atomic_mass_C"] + f["hydrogen_atoms"] * f["atomic_mass_H"])
    assert s["g_co2_per_km"] == pytest.approx(s["fuel_kg"] * 1e3 * ratio / s["distance_km"], rel=1e-9)
    # the engine never exceeds redline on the cycle and idles at the spec idle speed
    assert s["max_engine_rpm"] < load()["engine"]["redline_rpm"]
    idle = np.array([m == "idle" for m in res.mode])
    assert np.all(res.rpm[idle] == load()["engine"]["idle_rpm"])


def test_zero_to_100_plausible():
    """0-100 km/h of a ~110 kW / ~1350 kg FWD compact: 7.5-10.5 s (spec/vehicle.yaml plausibility);
    80-120 km/h in 5th is slower than 0-100 but finite (< 30 s)."""
    from vehicle.accel import run
    s, standing, elastic = run()
    lo, hi = VS["plausibility"]["zero_to_100_s"]
    assert lo <= s["zero_to_100_s"] <= hi
    assert s["zero_to_100_s"] < s["in_gear_80_120_s"] < 30.0
    assert np.all(np.diff(standing.v) > -1.0)          # speed never collapses (shift interruptions only)


# ------------------------------------------------------------------ replay data

@pytest.fixture(scope="module")
def replay():
    from viz.build import config
    cfg = config()
    path = ROOT / cfg["output"]["dir"] / cfg["output"]["data"]
    if not path.exists():
        from viz.build import build
        build(cfg)
    return cfg, json.loads(path.read_text())


def test_replay_data_consistent_with_kinematics(replay):
    """Replay json vs kin.slider_crank and the CAD placement: wrist-pin height, crankpin position,
    conrod angle and cylinder volume to 1e-6 (json rounding), rod length closed at every frame,
    crank rotation = -theta (ADR 0002), and the pressure / torque arrays cover the cycle."""
    cfg, d = replay
    geom = geometry(load())
    r, l = d["geometry"]["crank_radius"], d["geometry"]["conrod_length"]
    assert r == pytest.approx(geom.crank_radius, rel=1e-12) and l == pytest.approx(geom.conrod_length, rel=1e-12)
    th = np.array(d["theta_deg"])
    assert th[0] == 0.0 and np.allclose(np.diff(th), cfg["frames"]["step_deg"]) and th[-1] < 720.0
    assert np.allclose(d["crank_rotation_x"], -np.radians(th), atol=1e-6)
    for c, ph in enumerate(d["geometry"]["phase_deg"]):
        a = np.radians(th - ph)
        assert np.allclose(d["piston_z"][c], piston_position(a, r, l), atol=1e-6)
        assert np.allclose(d["pin_y"][c], r * np.sin(a), atol=1e-6)
        assert np.allclose(d["pin_z"][c], r * np.cos(a), atol=1e-6)
        assert np.allclose(d["rod_angle"][c], np.arcsin(r * np.sin(a) / l), atol=1e-6)
        rod = np.hypot(np.array(d["pin_y"][c]), np.array(d["piston_z"][c]) - np.array(d["pin_z"][c]))
        assert np.allclose(rod, l, atol=2e-6)
        assert np.allclose(np.array(d["volume_cc"][c]) * 1e-6, cylinder_volume(a, geom), rtol=1e-6)
        p = np.array(d["p_bar"][c])
        assert p.shape == th.shape and p.min() > 0.0
        # firing peak of cylinder c within 40 deg after its own firing TDC
        assert 0.0 <= (th[np.argmax(p)] - ph) % 720.0 < 40.0
    assert np.mean(d["torque_total_nm"]) == pytest.approx(d["torque_mean_nm"], rel=1e-6)
    assert d["modelled"] is True


def test_replay_page_is_self_contained_and_pinned(replay):
    """index.html embeds the same data and every CAD part, pins three.js to the exact configured version
    for both the module and the addons, and stays under the size budget."""
    cfg, d = replay
    out = ROOT / cfg["output"]["dir"]
    html = (out / cfg["output"]["html"]).read_text()
    ver = cfg["three"]["version"]
    assert f"three@{ver}/build/three.module.js" in html and f"three@{ver}/examples/jsm/" in html
    assert "three@latest" not in html and "__DATA_JSON__" not in html
    start = html.index('id="replay-data">') + len('id="replay-data">')
    embedded = json.loads(html[start:html.index("</script>", start)])
    assert embedded == d
    gstart = html.index('id="replay-glb">') + len('id="replay-glb">')
    glb = json.loads(html[gstart:html.index("</script>", gstart)])
    assert set(glb) == set(cfg["cad"]["parts"])
    size_mb = sum((out / cfg["output"][k]).stat().st_size for k in ("html", "data")) / 1e6
    assert size_mb < cfg["output"]["max_total_mb"]


def test_replay_loads_headless(replay):
    """Headless Chromium (Playwright via node): no console or page errors, WebGL up, the crank angle and
    frame counter advance between two samples; screenshot at out/viz/replay.png."""
    if shutil.which("node") is None:
        pytest.skip("node not available")
    if subprocess.run(["node", "-e", "require('playwright')"], capture_output=True).returncode != 0:
        pytest.skip("playwright (node) not available")
    from viz.build import check, fetch_vendor
    cfg, _ = replay
    try:
        fetch_vendor(cfg)
    except subprocess.CalledProcessError:
        pytest.skip("three.js files not reachable for the vendor cache")
    rep = check(cfg)
    assert rep["ok"], rep
    assert rep["webgl"] and rep["frames1"] > rep["frames0"] and rep["theta1"] != rep["theta0"]
