"""M2 exit criteria: two-zone Wiebe cycle with Cantera properties and Woschni heat transfer.

Verify: /home/user/redline-venv/bin/python -m pytest -q tests/test_m2_thermo.py
All engine numbers are MODELLED (no physical engine exists).
"""
import numpy as np
import pytest

from redline.spec import geometry, load
from thermo import gasexchange as gx
from thermo.cycle import energy_balance, otto_anchor, run_cycle
from thermo.ideal import otto_efficiency
from thermo.props import Direct, build_tables, burned, frozen, thermo_spec

# Stated tolerances (CLAUDE.md: never loosened silently)
TABLE_REL_TOL = 1e-3          # table interpolation vs direct Cantera
OTTO_REL_TOL = 1e-4           # analytic Otto anchor
ENERGY_REL_TOL = 5e-3         # closed-cycle energy balance (0.5 %)
MASS_REL_TOL = 1e-6           # mass conservation through the cycle
STEP_REL_TOL = 2e-3           # IMEP 0.1 vs 0.05 deg (0.2 %)
LHV_REL_TOL = 5e-3            # Cantera vapour LHV - h_fg vs spec (liquid) LHV
AFR_REL_TOL = 1e-2            # stoichiometric AFR from species data vs spec

# Plausible ranges for a modern NA SI engine at WOT, mid speed, CR ~10. Source: Heywood, "Internal
# Combustion Engine Fundamentals" (McGraw-Hill, 1988), ch. 9 (SI combustion: near-MBT timing places the
# 50 % mass-burned point roughly 5-10 deg ATDC; WOT peak pressures of several tens of bar) and ch. 12/13
# with sec. 5.7 (real-engine gross indicated fuel conversion efficiency well below the fuel-air-cycle
# value, typically ~0.33-0.42; full-load IMEP of NA SI engines ~10-13 bar). These are plausibility bands
# stated in the M2 plan, not calibration targets; CA50 band widened to 12 deg for the spec SOC.
P_MAX_BAR = (50.0, 80.0)
ETA_GROSS = (0.33, 0.42)
IMEP_BAR = (10.0, 13.0)
CA50_DEG = (6.0, 12.0)


@pytest.fixture(scope="module")
def tables():
    return build_tables()


@pytest.fixture(scope="module")
def wot(tables):
    return run_cycle(3000.0, "full", tables=tables)


# --------------------------------------------------------------------------------------------- tables
def test_frozen_tables_match_cantera(tables):
    d = Direct(tables=tables)
    rng = np.random.default_rng(1)
    for _ in range(40):
        T = rng.uniform(260.0, 3400.0)
        xr = rng.uniform(0.0, 1.0)
        h, cp, R = frozen(T, xr, *tables.args()[:5])
        h0, cp0, R0 = d.unburned(T, 1.0e5, xr)
        # enthalpy is absolute (can cross zero), so its error is scaled by the sensible scale cp*T
        assert abs(h - h0) / (cp0 * T) < TABLE_REL_TOL
        assert abs(cp - cp0) / cp0 < TABLE_REL_TOL
        assert abs(R - R0) / R0 < TABLE_REL_TOL


def test_burned_equilibrium_tables_match_cantera(tables):
    d = Direct(tables=tables)
    rng = np.random.default_rng(2)
    for _ in range(40):
        T = rng.uniform(800.0, 3400.0)
        p = np.exp(rng.uniform(np.log(0.5e5), np.log(1.5e7)))
        h, R, hT, hp, RT, Rp = burned(T, p, *tables.args()[5:])
        h0, R0 = d.burned(T, p)
        assert abs(h - h0) / (hT * T) < TABLE_REL_TOL
        assert abs(R - R0) / R0 < TABLE_REL_TOL


def test_fuel_data_consistent_with_spec(tables):
    spec, ts = load(), thermo_spec()
    lhv_liquid = tables.lhv_vapour - float(ts["mixture"]["fuel_h_fg"])
    assert abs(lhv_liquid / float(spec["combustion"]["lhv"]) - 1.0) < LHV_REL_TOL
    assert abs(tables.afr_stoich / float(spec["combustion"]["afr_stoich"]) - 1.0) < AFR_REL_TOL


# --------------------------------------------------------------------------------------------- (1)
def test_otto_anchor_instantaneous_and_wiebe():
    r = geometry().compression_ratio
    gamma = 1.4
    eta_inst, _ = otto_anchor(gamma=gamma, burn="instant")
    assert abs(eta_inst / otto_efficiency(r, gamma) - 1.0) < OTTO_REL_TOL
    eta_wiebe, _ = otto_anchor(gamma=gamma, burn="wiebe")
    assert eta_wiebe < eta_inst


# --------------------------------------------------------------------------------------------- (2)
def test_energy_balance_closes(wot, tables):
    eb = energy_balance(wot.closed, tables)
    assert abs(eb["rel_residual"]) < ENERGY_REL_TOL
    assert eb["Q_wall"] > 0.0 and eb["E_incomplete"] >= 0.0


# --------------------------------------------------------------------------------------------- (3)
def test_wot_3000_plausible_ranges(wot):
    m = wot.metrics
    assert m["converged"]
    assert P_MAX_BAR[0] <= m["p_max_bar"] <= P_MAX_BAR[1]
    assert ETA_GROSS[0] <= m["eta_ind_gross"] <= ETA_GROSS[1]
    assert IMEP_BAR[0] <= m["imep_gross_bar"] <= IMEP_BAR[1]
    assert IMEP_BAR[0] <= m["imep_net_bar"] <= IMEP_BAR[1]
    assert CA50_DEG[0] <= m["ca50_deg"] <= CA50_DEG[1]
    assert 0.0 < m["residual_fraction"] < 0.15
    assert 0.7 < m["volumetric_efficiency"] < 1.05


# --------------------------------------------------------------------------------------------- (4)
def test_mass_conservation(wot):
    m = wot.metrics
    oo = wot.open_out
    # gas exchange bookkeeping: cylinder mass change equals the sum of valve flows
    flows = (oo[-1, gx.GO_MF_IN] + oo[-1, gx.GO_MS_IN] + oo[-1, gx.GO_ME_IN]
             - oo[-1, gx.GO_MI_OUT] - oo[-1, gx.GO_ME_OUT])
    m_evo = oo[0, gx.GO_M]
    assert abs((oo[-1, gx.GO_M] - m_evo) - flows) / m_evo < MASS_REL_TOL
    # converged cycle: everything inducted from the manifold leaves through the exhaust
    assert abs(m["fresh_in_kg"] - m["exhaust_out_net_kg"]) / m["fresh_in_kg"] < MASS_REL_TOL
    # trapped mass is periodic (IVC of the final cycle vs the previous one)
    h = wot.history
    assert abs(h[-1]["m_ivc"] / h[-2]["m_ivc"] - 1.0) < MASS_REL_TOL
    # closed part: EVO mass equals IVC mass
    assert abs(m_evo / m["m_trapped_kg"] - 1.0) < MASS_REL_TOL


# --------------------------------------------------------------------------------------------- (5)
def test_step_size_convergence(wot, tables):
    fine = run_cycle(3000.0, "full", step_deg=0.05, tables=tables)
    for k in ("imep_gross_bar", "imep_net_bar"):
        assert abs(fine.metrics[k] / wot.metrics[k] - 1.0) < STEP_REL_TOL
