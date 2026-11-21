"""Full 720-degree single-cylinder cycle: closed two-zone part (IVC -> EVO) + filling/emptying gas exchange
(EVO -> IVC), iterated until the trapped state (mass, residual fraction, temperature) is periodic.

Crank angle follows ADR 0001 (0 = firing TDC, 360 = gas-exchange TDC). Internally the cycle runs from
IVC - 720 (negative, before firing TDC) to IVC. All results are MODELLED (there is no physical engine).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from kin.slider_crank import cylinder_volume, dV_dtheta
from redline.spec import geometry, load
from thermo import gasexchange as gx
from thermo import twozone as tz
from thermo.props import Direct, Tables, build_tables, frozen, burned, thermo_spec, T_from_u_frozen, ideal_tables
from thermo.woschni import heat_params

CYCLE_DEG = 720.0


# ---------------------------------------------------------------------------------------------
# prescribed laws
# ---------------------------------------------------------------------------------------------
def wiebe(theta_deg, soc_deg, dur_deg, a, m):
    """Wiebe mass fraction burned x_b and dx_b/dtheta [1/deg] (Heywood eq. 9.23)."""
    f = np.clip((np.asarray(theta_deg, dtype=float) - soc_deg) / dur_deg, 0.0, None)
    x = 1.0 - np.exp(-a * f ** (m + 1))
    dx = a * (m + 1) * f**m * np.exp(-a * f ** (m + 1)) / dur_deg
    return x, dx


def valve_area(theta_deg, open_deg, close_deg, n_valves, diameter, max_lift, cd, ld_cap):
    """Effective flow area Cd * n * pi * d * min(L, cap*d) with a sine-squared lift between open and close."""
    dur = (close_deg - open_deg) % CYCLE_DEG
    phi = ((np.asarray(theta_deg, dtype=float) - open_deg) % CYCLE_DEG) / dur
    lift = np.where(phi < 1.0, max_lift * np.sin(np.pi * np.clip(phi, 0.0, 1.0)) ** 2, 0.0)
    return cd * n_valves * np.pi * diameter * np.minimum(lift, ld_cap * diameter)


def _grid(theta_s, theta_e, step_deg):
    n = max(1, int(round((theta_e - theta_s) / step_deg)))
    hdeg = (theta_e - theta_s) / n
    half = theta_s + 0.5 * hdeg * np.arange(2 * n + 1)
    return n, hdeg, half


# ---------------------------------------------------------------------------------------------
# closed part
# ---------------------------------------------------------------------------------------------
@dataclass
class ClosedResult:
    theta: np.ndarray
    V: np.ndarray
    out: np.ndarray
    m: float
    xr: float
    T0: float

    def col(self, c):
        return self.out[:, c]


def run_closed(tables: Tables, geom, m, xr, T0, theta_s, theta_e, step_deg, hp, burn="wiebe",
               soc_deg=None, dur_deg=None, wiebe_a=None, wiebe_m=None, inst_deg=0.0, num=None):
    num = num or thermo_spec()["numerics"]
    n, hdeg, half = _grid(theta_s, theta_e, step_deg)
    th = np.radians(half)
    Vh = cylinder_volume(th, geom)
    dVh = dV_dtheta(th, geom)
    xb = np.zeros_like(half)
    dxb = np.zeros_like(half)
    k_act = k_inst = -1
    eps = 1e-9
    if burn == "wiebe":
        k_act = int(math.ceil((soc_deg - theta_s) / hdeg + 1.0 - eps))
        x, dx = wiebe(half, soc_deg, dur_deg, wiebe_a, wiebe_m)
        xb[2 * k_act:] = x[2 * k_act:]
        dxb[2 * k_act:] = np.degrees(dx[2 * k_act:])   # per rad
    elif burn == "instant":
        k_inst = int(round((inst_deg - theta_s) / hdeg))
        xb[2 * k_inst:] = 1.0
    elif burn != "none":
        raise ValueError(burn)
    _, _, Ru = frozen(T0, xr, tables.Fz, tables.T0f, tables.dTf, tables.R_F, tables.R_P)
    p0 = m * Ru * T0 / Vh[0]
    y0 = np.array([p0, T0, T0, 0.0, 0.0, 0.0])
    out = tz.integrate(y0, n, math.radians(hdeg), Vh, dVh, xb, dxb, m, xr, k_act, k_inst, *tables.args(), hp,
                       float(num["newton_tol"]), int(num["newton_max_iter"]), float(num["T_flame_guess"]))
    return ClosedResult(half[::2], Vh[::2], out, m, xr, T0)


# ---------------------------------------------------------------------------------------------
# full cycle
# ---------------------------------------------------------------------------------------------
@dataclass
class CycleResult:
    rpm: float
    load: str
    step_deg: float
    theta: np.ndarray          # deg, IVC-720 .. IVC
    V: np.ndarray
    p: np.ndarray
    T_u: np.ndarray            # unburned zone (closed) / bulk (open)
    T_b: np.ndarray            # burned zone (closed), NaN elsewhere
    T_bulk: np.ndarray
    x_b: np.ndarray
    hrr_chem: np.ndarray       # J/deg, fuel chemical energy release rate (Wiebe)
    hrr_apparent: np.ndarray   # J/deg, single-zone apparent net heat release from p, V
    closed: ClosedResult
    open_out: np.ndarray
    metrics: dict = field(default_factory=dict)
    history: list = field(default_factory=list)


def _apparent_hrr(theta_deg, p, V, gamma):
    dp = np.gradient(p, theta_deg)
    dV = np.gradient(V, theta_deg)
    return gamma / (gamma - 1.0) * p * dV + 1.0 / (gamma - 1.0) * V * dp


def run_cycle(rpm: float = 3000.0, load_name: str = "full", step_deg: float | None = None,
              heat_transfer: bool = True, tables: Tables | None = None, spec: dict | None = None,
              tspec: dict | None = None) -> CycleResult:
    spec = spec or load()
    tspec = tspec or thermo_spec()
    tables = tables or build_tables(tspec)
    geom = geometry(spec)
    num = tspec["numerics"]
    step_deg = float(step_deg or num["step_deg"])
    comb = spec["combustion"]
    val = spec["valves"]
    gex = tspec["gas_exchange"]
    bcs = tspec["loads"][load_name]
    tol, maxit = float(num["newton_tol"]), int(num["newton_max_iter"])
    targs = tables.args()

    ivc, evo = float(val["ivc_deg"]), float(val["evo_deg"])
    th_s, th_e = ivc - CYCLE_DEG, evo

    # gas-exchange grid and valve areas
    n_o, hdeg_o, half_o = _grid(evo, ivc, step_deg)
    tho = np.radians(half_o)
    Vo, dVo = cylinder_volume(tho, geom), dV_dtheta(tho, geom)
    nv = val["per_cylinder"]
    Ai = valve_area(half_o, val["ivo_deg"], ivc, nv["intake"], val["intake_head_diameter"],
                    val["intake_max_lift"], gex["cd_intake"], gex["lift_to_diameter_cap"])
    Ae = valve_area(half_o, evo, val["evc_deg"], nv["exhaust"], val["exhaust_head_diameter"],
                    val["exhaust_max_lift"], gex["cd_exhaust"], gex["lift_to_diameter_cap"])
    bc = np.array([bcs["p_intake"], bcs["T_intake"], bcs["p_exhaust"], bcs["T_exhaust"],
                   gex["linear_flow_pr"], tol, maxit, 0.0], dtype=float)
    _, _, R_F = frozen(bcs["T_intake"], 0.0, *targs[:5])
    bc[gx.BC_MEPS] = gex["store_blend_mass_fraction"] * bcs["p_intake"] * geom.swept_volume / (R_F * bcs["T_intake"])
    hp_open = heat_params(tspec, geom, rpm, "gas_exchange", on=heat_transfer)

    # initial guess at IVC
    xr = float(num["xr_guess"])
    T = float(num["T_ivc_guess"])
    V_ivc = cylinder_volume(math.radians(ivc), geom)
    _, _, Ru = frozen(T, xr, *targs[:5])
    m = bcs["p_intake"] * V_ivc / (Ru * T)
    store = np.zeros(3)
    history = []
    converged = False
    for it in range(int(num["max_cycles"])):
        hu, cpu, Ru = frozen(T, xr, *targs[:5])
        p_ivc = m * Ru * T / V_ivc
        hp_c = heat_params(tspec, geom, rpm, "closed", on=heat_transfer,
                           ref=(p_ivc, V_ivc, T, cpu / (cpu - Ru)))
        cl = run_closed(tables, geom, m, xr, T, th_s, th_e, step_deg, hp_c, "wiebe",
                        comb["heat_release_start_deg"], comb["burn_duration_deg"], comb["wiebe_a"],
                        comb["wiebe_m"], num=num)
        # EVO: merge zones at fixed total internal energy, burned gas relaxes to frozen products
        e = cl.out[-1]
        xb = e[tz.O_XB]
        mb, mu = xb * m, (1.0 - xb) * m
        hu, _, Ru = frozen(e[tz.O_TU], xr, *targs[:5])
        U = mu * (hu - Ru * e[tz.O_TU])
        if mb > 0.0:
            hb, Rb, *_ = burned(e[tz.O_TB], e[tz.O_P], *targs[5:])
            U += mb * (hb - Rb * e[tz.O_TB])
        x_open = (mb + mu * xr) / m
        T_open = T_from_u_frozen(U / m, x_open, e[tz.O_TB], *targs[:5], tol, maxit)
        y0 = np.zeros(gx.NG)
        y0[gx.G_M], y0[gx.G_MX], y0[gx.G_U] = m, m * x_open, U
        y0[gx.G_MS], y0[gx.G_MXS], y0[gx.G_HS] = store
        oo, yend = gx.integrate(y0, T_open, n_o, math.radians(hdeg_o), Vo, dVo, Ai, Ae, bc, *targs[:5], hp_open)
        m_new = yend[gx.G_M]
        xr_new = yend[gx.G_MX] / m_new
        T_new = oo[-1, gx.GO_T]
        store_new = np.array([yend[gx.G_MS], yend[gx.G_MXS], yend[gx.G_HS]])
        d = max(abs(m_new / m - 1.0), abs(xr_new - xr), abs(T_new / T - 1.0),
                abs(store_new[0] - store[0]) / m_new)
        store = store_new
        history.append({"iter": it, "m_ivc": m_new, "xr": xr_new, "T_ivc": T_new, "delta": d})
        last = (cl, oo, m, xr, T, yend)
        m, xr, T = m_new, xr_new, T_new
        if d < float(num["cycle_tol"]):
            converged = True
            break
    cl, oo, m_c, xr_c, T_c, yend = last

    # assemble 720-degree traces (drop duplicated EVO point)
    theta = np.concatenate([cl.theta, half_o[::2][1:]])
    V = np.concatenate([cl.V, Vo[::2][1:]])
    p = np.concatenate([cl.col(tz.O_P), oo[1:, gx.GO_P]])
    Tu = np.concatenate([cl.col(tz.O_TU), oo[1:, gx.GO_T]])
    Tb = np.concatenate([cl.col(tz.O_TB), np.full(n_o, np.nan)])
    xbc = cl.col(tz.O_XB)
    xbt = np.concatenate([xbc, np.full(n_o, np.nan)])
    # bulk temperature in the closed part from p V = m R_mix T
    Ru_c = frozen(T_c, xr_c, *targs[:5])[2]   # frozen mixture: R independent of T
    Rb_c = np.array([burned(tb, pp, *targs[5:])[1] if np.isfinite(tb) else 0.0
                     for tb, pp in zip(cl.col(tz.O_TB), cl.col(tz.O_P))])
    Tbulk_c = cl.col(tz.O_P) * cl.V / (m_c * ((1 - xbc) * Ru_c + xbc * Rb_c))
    Tbulk = np.concatenate([Tbulk_c, oo[1:, gx.GO_T]])

    m_fuel = m_c * (1.0 - xr_c) * tables.fuel_mass_fraction
    _, dx = wiebe(theta, comb["heat_release_start_deg"], comb["burn_duration_deg"], comb["wiebe_a"], comb["wiebe_m"])
    hrr_chem = np.where(theta <= evo, m_fuel * tables.lhv_vapour * dx, 0.0)
    hu, cpu, Ru = frozen(T_c, xr_c, *targs[:5])
    gamma_ref = cpu / (cpu - Ru)
    hrr_app = _apparent_hrr(theta, p, V, gamma_ref)

    res = CycleResult(rpm, load_name, step_deg, theta, V, p, Tu, Tb, Tbulk, xbt, hrr_chem, hrr_app, cl, oo,
                      history=history)
    res.metrics = _metrics(res, spec, tspec, tables, geom, m_c, xr_c, T_c, m_fuel, bcs, yend, converged,
                           hdeg_o, evo)
    return res


def _cum_at(theta0, hdeg, cum, angle):
    s = (angle - theta0) / hdeg
    i = min(max(int(math.floor(s)), 0), len(cum) - 2)
    f = s - i
    return cum[i] + f * (cum[i + 1] - cum[i])


def _metrics(res, spec, tspec, tables, geom, m, xr, T_ivc, m_fuel, bcs, yend, converged, hdeg_o, evo):
    cl, oo = res.closed, res.open_out
    Vd = geom.swept_volume
    lhv = float(spec["combustion"]["lhv"])
    W_closed = cl.col(tz.O_W)[-1]
    Wo = oo[:, gx.GO_W]
    bdc_c = -180.0 + CYCLE_DEG    # BDC before compression, in the gas-exchange window
    bdc_e = 180.0                 # BDC after expansion
    W_gross = (W_closed + _cum_at(evo, hdeg_o, Wo, bdc_e) - Wo[0]
               + Wo[-1] - _cum_at(evo, hdeg_o, Wo, bdc_c))
    W_net = W_closed + Wo[-1]
    pc = cl.col(tz.O_P)
    ipk = int(np.argmax(pc))
    xb = cl.col(tz.O_XB)
    ca50 = float(np.interp(0.5, xb[xb < 0.999], cl.theta[xb < 0.999]))
    ca10 = float(np.interp(0.1, xb[xb < 0.999], cl.theta[xb < 0.999]))
    ca90 = float(np.interp(0.9, xb[xb < 0.999], cl.theta[xb < 0.999]))
    rho_air = bcs["p_intake"] / (tables.R_air * bcs["T_intake"])
    fresh_in = yend[gx.G_MF_IN]
    m_air_in = fresh_in * (1.0 - tables.fuel_mass_fraction)
    ex_net = yend[gx.G_ME_OUT] - yend[gx.G_ME_IN]
    Q_closed = cl.col(tz.O_QU)[-1] + cl.col(tz.O_QB)[-1]
    Q_open = oo[-1, gx.GO_Q]
    return {
        "converged": converged,
        "cycles": len(res.history),
        "imep_gross_bar": W_gross / Vd / 1e5,
        "imep_net_bar": W_net / Vd / 1e5,
        "pmep_bar": (W_gross - W_net) / Vd / 1e5,
        "W_gross_J": W_gross,
        "W_net_J": W_net,
        "eta_ind_gross": W_gross / (m_fuel * lhv),
        "eta_ind_net": W_net / (m_fuel * lhv),
        "p_max_bar": pc[ipk] / 1e5,
        "theta_p_max_deg": float(cl.theta[ipk]),
        "ca10_deg": ca10,
        "ca50_deg": ca50,
        "ca90_deg": ca90,
        "T_b_max_K": float(np.nanmax(cl.col(tz.O_TB))),
        "T_u_max_K": float(np.max(cl.col(tz.O_TU))),
        "residual_fraction": xr,
        "volumetric_efficiency": m_air_in / (rho_air * Vd),
        "m_trapped_kg": m,
        "m_fuel_kg": m_fuel,
        "T_ivc_K": T_ivc,
        "p_ivc_bar": pc[0] / 1e5,
        "Q_wall_closed_J": Q_closed,
        "Q_wall_open_J": Q_open,
        "fresh_in_kg": fresh_in,
        "exhaust_out_net_kg": ex_net,
        "port_store_kg": yend[gx.G_MS],
        "lhv_spec": lhv,
        "lhv_vapour_cantera": tables.lhv_vapour,
    }


# ---------------------------------------------------------------------------------------------
# energy balance of the closed part
# ---------------------------------------------------------------------------------------------
def energy_balance(cl: ClosedResult, tables: Tables, tspec: dict | None = None, direct: Direct | None = None):
    """Fuel chemical energy = W + Q_wall + dU_sensible + incomplete combustion, IVC -> EVO.

    End-state energies are evaluated with direct Cantera calls (independent of the tables). Fuel chemical
    energy uses the constant-volume heating value at T_ref (closed system, Heywood sec. 3.5).
    """
    tspec = tspec or thermo_spec()
    direct = direct or Direct(tspec, tables)
    Tref = float(tspec["mixture"]["T_ref"])
    pref = float(tspec["mixture"]["p_ref"])
    m, xr = cl.m, cl.xr
    e = cl.out[-1]
    xb = e[tz.O_XB]
    mb, mu = xb * m, (1.0 - xb) * m

    def u_unb(T):
        h, cp, R = direct.unburned(T, pref, xr)
        return h - R * T

    def u_prod(T):
        h, cp, R = direct.unburned(T, pref, 1.0)
        return h - R * T

    hb, Rb = direct.burned(e[tz.O_TB], e[tz.O_P])
    u_beq = hb - Rb * e[tz.O_TB]
    m_fuel = m * (1.0 - xr) * tables.fuel_mass_fraction
    E_fuel = m_fuel * tables.lhv_cv_vapour
    W = e[tz.O_W]
    Q = e[tz.O_QU] + e[tz.O_QB]
    U_s0 = m * (u_unb(cl.T0) - u_unb(Tref))
    U_s1 = mu * (u_unb(e[tz.O_TU]) - u_unb(Tref)) + mb * (u_prod(e[tz.O_TB]) - u_prod(Tref))
    E_inc = mb * (u_beq - u_prod(e[tz.O_TB])) + mu * (u_unb(Tref) - u_prod(Tref))
    resid = E_fuel - (W + Q + (U_s1 - U_s0) + E_inc)
    return {"E_fuel": E_fuel, "W": W, "Q_wall": Q, "dU_sensible": U_s1 - U_s0, "E_incomplete": E_inc,
            "residual": resid, "rel_residual": resid / E_fuel}


# ---------------------------------------------------------------------------------------------
# analytic anchor
# ---------------------------------------------------------------------------------------------
def otto_anchor(gamma=1.4, R=287.0, p1=1.0e5, T1=300.0, q_mix=2.0e6, step_deg=0.1, burn="instant",
                spec=None, tspec=None):
    """Two-zone integrator with constant-gamma, identical-property tables (collapses to single zone),
    heat transfer off, BDC -> BDC. Returns (efficiency, ClosedResult)."""
    spec = spec or load()
    tspec = tspec or thermo_spec()
    geom = geometry(spec)
    tab = ideal_tables(gamma, R, q_mix)
    comb = spec["combustion"]
    hp = heat_params(tspec, geom, 1.0, "closed", on=False)
    V1 = cylinder_volume(-math.pi, geom)
    m = p1 * V1 / (R * T1)
    cl = run_closed(tab, geom, m, 0.0, T1, -180.0, 180.0, step_deg, hp, burn,
                    comb["heat_release_start_deg"], comb["burn_duration_deg"], comb["wiebe_a"], comb["wiebe_m"],
                    inst_deg=0.0, num=tspec["numerics"])
    W = cl.col(tz.O_W)[-1]
    xb_end = cl.col(tz.O_XB)[-1]
    return W / (m * q_mix * xb_end), cl
