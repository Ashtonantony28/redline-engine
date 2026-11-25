"""M7: knock, quasi-dimensional flame model, emissions, spark optimisation and speed-load maps.

Verify: python -m pytest -q tests/test_m7_maps.py
        (the fixture regeneration test is marked slow, ~30 s; deselect with -m "not slow")
        python -m sim.maps
All engine numbers are MODELLED (no physical engine exists, CLAUDE.md rule 5).

Analytic / published anchors come first, then the M7 exit criteria (PLAN.md sec. 5).
"""
import copy
import math

import numpy as np
import pytest

from gasdyn.engine import EngineModel
from gasdyn.spec import gasdyn_spec
from redline.spec import geometry, load
from thermo import twozone as tz
from thermo.cycle import run_cycle
from thermo.emissions import Emissions, EquilibriumTable, integrate_no_fixed, zeldovich_params
from thermo.flame import laminar_speed, mk_params, sphere_cyl_area, sphere_cyl_volume
from thermo.knock import klsa_search, livengood_wu, tau_douaud_eyzat
from thermo.m7spec import combustion_spec
from thermo.props import build_tables, thermo_spec
from thermo.spark import SparkModel, Trapped

# ------------------------------------------------------------------ stated tolerances (never loosened silently)
LW_REL_TOL = 1e-12            # Livengood-Wu integral vs t / tau at constant state (trapezoid is exact)
SL_REL_TOL = 1e-12            # Metghalchi-Keck reference point (closed form)
GEOM_REL_TOL = 1e-12          # flame geometry closed forms vs elementary shapes
GEOM_QUAD_REL_TOL = 1e-6      # clipped flame volume vs numerical quadrature
NO_EQ_TABLE_REL_TOL = 1e-2    # tabulated (ln X bilinear) equilibrium NO vs direct Cantera, off-grid point
NO_KIN_REL_TOL = 5e-3         # Zeldovich kinetics at fixed T, p after many relaxation times vs [NO]_e
NO_RATE_REL_TOL = 1e-3        # first backward-Euler step (1 ns) vs the initial rate 2 R1 (Heywood eq. 11.9)
# Flame calibration (spec flame.entrainment_coeff): 3000 rpm WOT burn vs the engine.yaml Wiebe burn
CAL_CA50_TOL_DEG = 1.5
CAL_1090_TOL_DEG = 3.0
# Burn duration (spark -> CA90, crank degrees) vs speed: grows, but less than linearly. Log-log slope band.
DUR_SLOPE = (0.0, 0.5)
# MBT CA50: Heywood 1988 sec. 9.2.2 (MBT timing puts the 50 % mass-burned point roughly 6-10 deg ATDC,
# peak pressure ~15 deg ATDC). Stated model tolerance: 3 deg below the cited band. The calibrated flame
# model's burn-up tail (tau_b = lambda_T / S_L) is short, so the burn is front-loaded towards CA50 and its
# heat-release centroid lies later than CA50; at high speed and part load (lower heat loss per cycle)
# the modelled MBT CA50 therefore sits at 3.4-6 deg ATDC instead of 6-10. Not tuned away; reported.
MBT_CA50_BAND = (3.0, 12.0)
# Minimum BSFC of a modern NA port-injected 2.0 L SI engine: ~235-260 g/kWh (Heywood 1988 Fig. 15-9 shows
# best BSFC near 250 g/kWh for older SI engines; the PLAN target is 245). Plausibility band:
BSFC_MIN_BAND = (230.0, 270.0)
CSV_ABS_TOL = 1e-3            # fixture values are written with 6 significant digits

N_DUR_SPEEDS = (1000.0, 3000.0, 7000.0)


@pytest.fixture(scope="module")
def cspec():
    return combustion_spec()


@pytest.fixture(scope="module")
def tables():
    return build_tables()


@pytest.fixture(scope="module")
def eqtab(cspec, tables):
    return EquilibriumTable(thermo_spec(), cspec, tables)


@pytest.fixture(scope="module")
def engine(tables):
    return EngineModel(tables=tables)


@pytest.fixture(scope="module")
def sm(engine, cspec):
    return SparkModel(engine.spec, engine.tspec, cspec, engine.gspec, engine.tables, engine.geom)


@pytest.fixture(scope="module")
def wot3000(engine):
    """M4 3000 rpm WOT point with the spec Wiebe burn and its cylinder-mean trapped state."""
    op = engine.run(3000.0)
    rec = op.per_cylinder
    mean = lambda k: float(np.mean([r[k] for r in rec]))
    st = Trapped(mean("m_trapped"), mean("xr"), mean("T_ivc"), engine.ivc_deg, engine.evo_deg,
                 float(engine.gspec["numerics"]["closed_step_deg"]))
    return op, st


# =========================================================================== analytic anchors
def test_livengood_wu_constant_state_is_t_over_tau(cspec):
    k = cspec["knock"]
    p, T = 40.0e5, 850.0
    t = np.linspace(0.0, 3.0e-3, 301)
    tau = tau_douaud_eyzat(np.full_like(t, p), np.full_like(t, T), k)
    I = livengood_wu(t, tau)
    assert np.allclose(I, t / tau[0], rtol=LW_REL_TOL, atol=0.0)
    # Douaud-Eyzat closed form at ON = ON_ref, 1 atm: tau = A exp(T_a / T)
    k2 = copy.deepcopy(k)
    k2["octane_number"] = k["douaud_eyzat"]["on_ref"]
    de = k["douaud_eyzat"]
    t0 = float(tau_douaud_eyzat(de["p_unit"], 1000.0, k2))
    assert t0 == pytest.approx(de["A"] * math.exp(de["T_activation"] / 1000.0), rel=LW_REL_TOL)


def test_klsa_search_finds_first_crossing():
    # synthetic non-monotonic margin: rises through 1 at 12.3 deg, falls again beyond 30 deg
    f = lambda a: 1.0 + 0.05 * (a - 12.3) if a < 30.0 else 0.5
    klsa, status = klsa_search(f, -10.0, 50.0, 1e-3, 1.0, 2.5)
    assert status == "limited" and klsa == pytest.approx(12.3, abs=1e-3)
    assert klsa_search(lambda a: 0.1, -10.0, 50.0, 1e-3, 1.0, 2.5) == (50.0, "unlimited")
    assert klsa_search(lambda a: 2.0, -10.0, 50.0, 1e-3, 1.0, 2.5) == (-10.0, "unavoidable")


def test_laminar_speed_reproduces_metghalchi_keck_reference(cspec):
    lam = cspec["flame"]["laminar"]
    mk = mk_params(cspec, 1.0)
    s_ref = lam["B_M"] + lam["B_phi"] * (1.0 - lam["phi_M"]) ** 2       # 24.89 cm/s for iso-octane
    assert laminar_speed(lam["T0"], lam["p0"], 0.0, mk) == pytest.approx(s_ref, rel=SL_REL_TOL)
    assert s_ref == pytest.approx(0.2489, abs=1e-4)
    # power laws and dilution factor at phi = 1
    assert laminar_speed(2 * lam["T0"], lam["p0"], 0.0, mk) / s_ref == pytest.approx(2 ** lam["alpha_0"], rel=SL_REL_TOL)
    assert laminar_speed(lam["T0"], 10 * lam["p0"], 0.0, mk) / s_ref == pytest.approx(10 ** lam["beta_0"], rel=SL_REL_TOL)
    x = 0.1
    assert laminar_speed(lam["T0"], lam["p0"], x, mk) / s_ref == pytest.approx(
        1.0 - lam["dilution_coeff"] * x ** lam["dilution_exp"], rel=SL_REL_TOL)


def test_flame_geometry_closed_forms():
    R, h = 0.043, 0.02
    r = 0.01      # free hemisphere
    assert sphere_cyl_volume(r, h, R) == pytest.approx(2.0 / 3.0 * math.pi * r**3, rel=GEOM_REL_TOL)
    assert sphere_cyl_area(r, h, R) == pytest.approx(2.0 * math.pi * r**2, rel=GEOM_REL_TOL)
    big = 1.0     # sphere engulfs the chamber
    assert sphere_cyl_volume(big, h, R) == pytest.approx(math.pi * R**2 * h, rel=GEOM_REL_TOL)
    assert sphere_cyl_area(big, h, R) == 0.0
    # clipped by both walls: quadrature of the disc slices
    r = 0.047
    z = np.linspace(0.0, min(r, h), 200001)
    a2 = np.minimum(r**2 - z**2, R**2)
    vq = np.trapezoid(math.pi * a2, z)
    assert sphere_cyl_volume(r, h, R) == pytest.approx(vq, rel=GEOM_QUAD_REL_TOL)


def test_zeldovich_equilibrium_and_kinetics(cspec, tables, eqtab):
    T, p = 2537.0, 43.0e5          # off the table grid in both T and p
    x_tab = eqtab.X("NO", T, p)
    x_ct = eqtab.direct_X("NO", T, p)
    assert x_tab == pytest.approx(x_ct, rel=NO_EQ_TABLE_REL_TOL)
    em = Emissions(thermo_spec(), cspec, tables, eqtab)
    kp = zeldovich_params(cspec)
    # many NO relaxation times (tau_NO ~ 1e-4 s here; Heywood eq. 11.10)
    no, NOe, R1 = integrate_no_fixed(T, p, 0.2, 20000, *eqtab.args(), *em.idx, kp, em.tol, em.maxit)
    assert no[-1] == pytest.approx(NOe, rel=NO_KIN_REL_TOL)
    assert no[-1] <= NOe * (1.0 + 1e-12)
    assert np.all(np.diff(no) >= 0.0)
    # [NO]_e consistent with the direct Cantera mole fraction
    assert NOe == pytest.approx(x_ct * p / (8.314462618 * T), rel=NO_EQ_TABLE_REL_TOL)
    # initial rate = 2 R1 (Heywood eq. 11.9); first backward-Euler step of a tiny dt
    dt = 1.0e-9
    no1, _, _ = integrate_no_fixed(T, p, dt, 1, *eqtab.args(), *em.idx, kp, em.tol, em.maxit)
    assert no1[-1] / dt == pytest.approx(2.0 * R1, rel=NO_RATE_REL_TOL)


# =========================================================================== flame model
def test_flame_calibrated_to_wiebe_at_3000_wot(sm, wot3000):
    op, st = wot3000
    comb = load()["combustion"]
    a, m_, dur, soc = comb["wiebe_a"], comb["wiebe_m"], comb["burn_duration_deg"], comb["heat_release_start_deg"]
    f = lambda x: soc + dur * (-math.log(1.0 - x) / a) ** (1.0 / (m_ + 1.0))
    ev = sm.evaluate(st, 3000.0, comb["spark_advance_deg"])
    assert ev["ca50"] == pytest.approx(f(0.5), abs=CAL_CA50_TOL_DEG)
    assert ev["ca90"] - ev["ca10"] == pytest.approx(f(0.9) - f(0.1), abs=CAL_1090_TOL_DEG)
    assert ev["cl"].col(tz.O_XB)[-1] > 0.999


def test_burn_duration_grows_less_than_linearly_with_speed(sm, wot3000):
    _, st = wot3000
    adv = load()["combustion"]["spark_advance_deg"]
    dur = [sm.evaluate(st, n, adv)["ca90"] + adv for n in N_DUR_SPEEDS]
    assert all(np.diff(dur) > 0.0)
    slope = np.polyfit(np.log(N_DUR_SPEEDS), np.log(dur), 1)[0]
    assert DUR_SLOPE[0] < slope < DUR_SLOPE[1]


# =========================================================================== exit criteria
def _klsa(sm_, st, rpm, cspec):
    s = cspec["spark"]
    return klsa_search(lambda a: sm_.evaluate(st, rpm, a)["knock"]["margin"], s["advance_min_deg"],
                       s["advance_max_deg"], s["klsa_tol_deg"], cspec["knock"]["integral_threshold"],
                       s["klsa_scan_step_deg"])


def test_knock_limit_moves_with_cr_and_intake_temperature(cspec, tables):
    """Quasi-steady WOT cycle at low speed: KLSA falls with higher CR and with hotter intake air."""
    rpm = 1500.0
    tspec0 = thermo_spec()
    qs = cspec["maps"]["quasi_steady"]
    base = load()

    def klsa(cr_delta, dT):
        spec = copy.deepcopy(base)
        spec["engine"]["compression_ratio"] += cr_delta
        ts = copy.deepcopy(tspec0)
        ts["loads"]["m7"] = dict(ts["loads"]["full"], T_intake=float(qs["T_intake"]) + dT)
        r = run_cycle(rpm, "m7", float(qs["step_deg"]), tables=tables, spec=spec, tspec=ts)
        val = spec["valves"]
        st = Trapped(r.metrics["m_trapped_kg"], r.metrics["residual_fraction"], r.metrics["T_ivc_K"],
                     float(val["ivc_deg"]), float(val["evo_deg"]), float(qs["step_deg"]))
        s = SparkModel(spec, ts, cspec, gasdyn_spec(), tables, geometry(spec))
        return _klsa(s, st, rpm, cspec)

    k0, s0 = klsa(0.0, 0.0)
    k_cr, s_cr = klsa(1.0, 0.0)
    k_t, s_t = klsa(0.0, 30.0)
    assert s0 == "limited"
    assert k_cr < k0
    assert k_t < k0


def test_mbt_ca50_and_nox_rise_with_advance(sm, wot3000, cspec, tables):
    _, st = wot3000
    opt = sm.optimise(st, 3000.0)
    assert not opt["knock_limited"]
    assert MBT_CA50_BAND[0] <= opt["ca50_mbt"] <= MBT_CA50_BAND[1]
    em = Emissions(thermo_spec(), cspec, tables)
    nox = [em.closed(sm.evaluate(st, 3000.0, a)["cl"], 3000.0)["m_NOx_kg"] for a in (15.0, 25.0, 35.0)]
    assert nox[0] < nox[1] < nox[2]


# ---------------------------------------------------------------------------- committed maps fixture
@pytest.fixture(scope="module")
def fixture_rows():
    from sim.maps import COLUMNS, FIXTURE, read_csv
    rows = read_csv(FIXTURE)
    assert list(rows[0].keys()) == COLUMNS
    return rows


def test_maps_fixture_covers_grid_and_converged(fixture_rows, cspec):
    from sim.maps import grid, level_name
    keys = {(float(r["rpm"]), r["level"]) for r in fixture_rows}
    want = {(r, level_name(l)) for r, l in grid(cspec, False)}
    assert keys <= want
    # only map levels above the WOT plenum pressure may be missing, never WOT itself
    assert all((r, "wot") in keys for r, _ in want)
    assert all(r["converged"] == "True" for r in fixture_rows)


def test_maps_fixture_exit_criteria(fixture_rows, cspec):
    f = lambda r, k: float(r[k])
    margin = cspec["knock"]["klsa_margin_deg"]
    adv_min = cspec["spark"]["advance_min_deg"]
    bsfc = [f(r, "bsfc_g_per_kwh") for r in fixture_rows if np.isfinite(f(r, "bsfc_g_per_kwh"))]
    assert BSFC_MIN_BAND[0] <= min(bsfc) <= BSFC_MIN_BAND[1]
    for r in fixture_rows:
        if r["knock_limited"] == "False":
            assert MBT_CA50_BAND[0] <= f(r, "ca50_deg") <= MBT_CA50_BAND[1], (r["rpm"], r["level"])
        else:
            assert f(r, "spark_deg") < f(r, "mbt_deg")
            assert f(r, "spark_deg") == pytest.approx(max(f(r, "klsa_deg") - margin, adv_min), abs=CSV_ABS_TOL)
    # knock-limited region sits at low speed and high load (the fixture must show one)
    kl = [r for r in fixture_rows if r["knock_limited"] == "True"]
    assert kl and all(r["level"] in ("wot", "0.90bar", "0.80bar") for r in kl)


@pytest.mark.slow
def test_maps_fixture_regenerates_identically(fixture_rows, cspec):
    from sim.maps import COLUMNS, csv_rows, level_name, run_points
    pts = [(float(p["rpm"]), p["level"] if p["level"] == "wot" else float(p["level"]))
           for p in cspec["maps"]["fixture_points"]]
    rows, skipped = run_points(pts, False, jobs=1)
    assert not skipped
    by_key = {(float(r["rpm"]), r["level"]): [r[c] for c in COLUMNS] for r in fixture_rows}
    for (rpm, lv), new in zip(pts, csv_rows(rows)):
        assert new == by_key[(rpm, level_name(lv))]
