"""Burned-zone emissions (milestone M7): extended Zeldovich NO kinetics frozen at EVO, equilibrium CO.

Equilibrium burned-gas mole fractions come from a Cantera table on the spec/thermo.yaml burned-gas grid
(T, ln p), interpolated bilinearly in ln X. NO follows Heywood 1988 eq. 11.8,

    d[NO]/dt = 2 R1 (1 - alpha^2) / (1 + alpha R1 / (R2 + R3)),   alpha = [NO] / [NO]_e
    R1 = k1f [O]_e [N2]_e,  R2 = k2r [NO]_e [O]_e,  R3 = k3r [NO]_e [H]_e,

with the Table 11.1 rate constants (spec/combustion_ext.yaml). The burned zone is split into equal-mass
elements burned in sequence; each starts at the constant-pressure adiabatic flame temperature of the
unburned gas it came from and is then compressed / expanded isentropically on the equilibrium burned-gas
tables, with temperatures scaled so the element mass average equals the two-zone burned-zone temperature
(heat loss is shared in proportion). NO is integrated per element with backward Euler (stiff near
3000 K) and frozen at EVO. CO is reported at equilibrium at EVO and, as an upper bound, at the peak
burned-zone temperature. MODELLED results only.
"""
from __future__ import annotations

import hashlib
import json
import math

import cantera as ct
import numpy as np
from numba import njit

from thermo import twozone as tz
from thermo.props import CACHE_DIR, T_adiabatic_p, burned, frozen, make_gas

_TABLE_VERSION = "m7-eq-1"
R_MOL = ct.gas_constant / 1.0e3      # J/(mol K)
KG_PER_G = 1.0e-3
J_PER_KWH = 3.6e6
SEC_PER_MIN = 60.0
REV_PER_CYCLE = 2.0     # four-stroke


class EquilibriumTable:
    """ln(mole fraction) of the burned-gas species and mean molar mass [kg/mol] on (T, ln p)."""

    def __init__(self, tspec: dict, cspec: dict, tables, use_cache: bool = True):
        tab = tspec["tables"]
        em = cspec["emissions"]
        self.species = list(em["table_species"])
        self.floor = float(em["log_floor"])
        self.T0, self.dT = float(tab["Tb_min"]), float(tab["dT_burned"])
        T = np.arange(self.T0, float(tab["Tb_max"]) + 0.5 * self.dT, self.dT)
        lp = np.linspace(np.log(float(tab["p_min"])), np.log(float(tab["p_max"])), int(tab["n_p"]))
        self.lp0, self.dlp = float(lp[0]), float(lp[1] - lp[0])
        key = hashlib.sha256(json.dumps({"tab": tab, "mix": tspec["mixture"], "sp": self.species,
                                         "v": _TABLE_VERSION, "ct": ct.__version__}, sort_keys=True).encode()
                             ).hexdigest()[:16]
        cache = CACHE_DIR / f"eq_{key}.npz"
        self.gas = make_gas(tspec)
        self.Y = tables.Y_P
        if use_cache and cache.exists():
            self.L = np.load(cache)["L"]
        else:
            L = np.empty((len(self.species) + 1, T.size, lp.size))
            idx = [self.gas.species_index(s) for s in self.species]
            for j, l in enumerate(lp):
                self.gas.TPY = T[0], math.exp(l), self.Y
                for i, t in enumerate(T):
                    self.gas.TP = t, math.exp(l)
                    self.gas.equilibrate("TP")
                    X = self.gas.X
                    for k, n in enumerate(idx):
                        L[k, i, j] = math.log(max(X[n], self.floor))
                    L[-1, i, j] = self.gas.mean_molecular_weight / 1.0e3
            self.L = L
            if use_cache:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                np.savez(cache, L=L)
        self.ix = {s: k for k, s in enumerate(self.species)}

    def args(self):
        return self.L, self.T0, self.dT, self.lp0, self.dlp

    def X(self, species: str, T: float, p: float) -> float:
        out = np.empty(self.L.shape[0])
        _interp(T, p, *self.args(), out)
        return float(math.exp(out[self.ix[species]]))

    def direct_X(self, species: str, T: float, p: float) -> float:
        """Direct Cantera equilibrium (reference for the tests)."""
        self.gas.TPY = T, p, self.Y
        self.gas.equilibrate("TP")
        return float(self.gas.X[self.gas.species_index(species)])


@njit(cache=True)
def _interp(T, p, L, T0, dT, lp0, dlp, out):
    nT, nP = L.shape[1], L.shape[2]
    s = (T - T0) / dT
    i = min(max(int(math.floor(s)), 0), nT - 2)
    f = s - i
    q = (math.log(p) - lp0) / dlp
    j = min(max(int(math.floor(q)), 0), nP - 2)
    g = q - j
    for k in range(L.shape[0] - 1):
        out[k] = ((1 - f) * (1 - g) * L[k, i, j] + f * (1 - g) * L[k, i + 1, j]
                  + (1 - f) * g * L[k, i, j + 1] + f * g * L[k, i + 1, j + 1])
    k = L.shape[0] - 1
    out[k] = ((1 - f) * (1 - g) * L[k, i, j] + f * (1 - g) * L[k, i + 1, j]
              + (1 - f) * g * L[k, i, j + 1] + f * g * L[k, i + 1, j + 1])


# rate-constant vector: (A, n, E) for k1f, k1r, k2f, k2r, k3f, k3r, already in m^3/(mol s)
def zeldovich_params(cspec: dict) -> np.ndarray:
    z = cspec["emissions"]["zeldovich"]
    c = float(z["cm3_to_m3"])
    rows = []
    for name in ("k1f", "k1r", "k2f", "k2r", "k3f", "k3r"):
        rows.append([float(z[name]["A"]) * c, float(z[name]["n"]), float(z[name]["E"])])
    return np.array(rows)


@njit(cache=True)
def _k(kp, r, T):
    return kp[r, 0] * T ** kp[r, 1] * math.exp(-kp[r, 2] / T)


@njit(cache=True)
def zeldovich_terms(T, p, lnX, iN2, iO, iOH, iH, iNO, kp):
    """(R1, R2 + R3, [NO]_e) [mol/(m^3 s), mol/m^3] at equilibrium (Heywood eq. 11.8)."""
    c = p / (R_MOL * T)
    N2 = c * math.exp(lnX[iN2])
    O = c * math.exp(lnX[iO])
    OH = c * math.exp(lnX[iOH])
    H = c * math.exp(lnX[iH])
    NOe = c * math.exp(lnX[iNO])
    R1 = _k(kp, 0, T) * O * N2
    R23 = _k(kp, 3, T) * NOe * O + _k(kp, 5, T) * NOe * H
    return R1, R23, NOe


@njit(cache=True)
def no_rate(NO, R1, R23, NOe):
    a = NO / NOe
    return 2.0 * R1 * (1.0 - a * a) / (1.0 + a * R1 / R23)


@njit(cache=True)
def _be_step(NO_old, dt, R1, R23, NOe, tol, maxit):
    """Backward Euler for d[NO]/dt at fixed (T, p) over dt: solve x - NO_old - dt f(x) = 0 (Newton)."""
    x = NO_old
    for _ in range(maxit):
        a = x / NOe
        den = 1.0 + a * R1 / R23
        f = 2.0 * R1 * (1.0 - a * a) / den
        df = (2.0 * R1 * (-2.0 * a) / den - 2.0 * R1 * (1.0 - a * a) * (R1 / R23) / (den * den)) / NOe
        r = x - NO_old - dt * f
        dx = r / (1.0 - dt * df)
        x -= dx
        if x < 0.0:
            x = 0.0
        if abs(dx) <= tol * max(NOe, 1e-300):
            break
    return x


@njit(cache=True)
def integrate_no_fixed(T, p, t_end, n, L, T0, dT, lp0, dlp, iN2, iO, iOH, iH, iNO, kp, tol, maxit):
    """[NO](t) at constant T, p from zero (anchor test): n backward-Euler steps up to t_end."""
    lnX = np.empty(L.shape[0])
    _interp(T, p, L, T0, dT, lp0, dlp, lnX)
    R1, R23, NOe = zeldovich_terms(T, p, lnX, iN2, iO, iOH, iH, iNO, kp)
    out = np.zeros(n + 1)
    dt = t_end / n
    for i in range(n):
        out[i + 1] = _be_step(out[i], dt, R1, R23, NOe, tol, maxit)
    return out, NOe, R1


@njit(cache=True)
def _elements(theta, p, Tu, Tb, xb, xr, omega, n_el, Fz, T0f, dTf, R_F, R_P, B, T0b, dTb, lp0, dlp,
              L, LT0, LdT, Llp0, Ldlp, iN2, iO, iOH, iH, iNO, kp, tol, maxit, Tguess, MW_NO):
    """Per-element NO mass fraction at the end of the closed part (EVO). Returns (y_NO[n_el], formed[n_el],
    element peak temperature[n_el])."""
    n = theta.size
    xend = xb[n - 1]
    y = np.zeros(n_el)
    Te = np.zeros(n_el)
    Tpk = np.zeros(n_el)
    formed = np.zeros(n_el, dtype=np.bool_)
    lnX = np.empty(L.shape[0])
    nf = 0
    for i in range(n):
        if not np.isfinite(Tb[i]):
            continue
        # isentropic update of formed elements over the step i-1 -> i (midpoint rule)
        if i > 0 and nf > 0:
            dp = p[i] - p[i - 1]
            for j in range(nf):
                T = Te[j]
                h, R, hT, hp, RT, Rp = burned(T, p[i - 1], B, T0b, dTb, lp0, dlp)
                Tm = T + 0.5 * (R * T / p[i - 1] - hp) / hT * dp
                pm = 0.5 * (p[i - 1] + p[i])
                h, R, hT, hp, RT, Rp = burned(Tm, pm, B, T0b, dTb, lp0, dlp)
                Te[j] = T + (R * Tm / pm - hp) / hT * dp
        # new elements whose formation (mid-slice) burned fraction has been reached
        while nf < n_el and xb[i] >= (nf + 0.5) / n_el * xend:
            hu, cpu, Ru = frozen(Tu[i], xr, Fz, T0f, dTf, R_F, R_P)
            Te[nf] = T_adiabatic_p(hu, p[i], Tguess, B, T0b, dTb, lp0, dlp, tol, maxit)
            formed[nf] = True
            nf += 1
        if nf == 0 or i == 0:
            continue
        # scale element temperatures so their mean is the two-zone burned-zone temperature
        s = 0.0
        for j in range(nf):
            s += Te[j]
        scale = Tb[i] / (s / nf)
        dt = (theta[i] - theta[i - 1]) * math.pi / 180.0 / omega
        for j in range(nf):
            T = Te[j] * scale
            if T > Tpk[j]:
                Tpk[j] = T
            _interp(T, p[i], L, LT0, LdT, Llp0, Ldlp, lnX)
            R1, R23, NOe = zeldovich_terms(T, p[i], lnX, iN2, iO, iOH, iH, iNO, kp)
            h, R, hT, hp, RT, Rp = burned(T, p[i], B, T0b, dTb, lp0, dlp)
            rho = p[i] / (R * T)
            c_old = y[j] * rho / MW_NO
            c_new = _be_step(c_old, dt, R1, R23, NOe, tol, maxit)
            y[j] = c_new * MW_NO / rho
    return y, formed, Tpk


class Emissions:
    MW = {"NO": 30.006e-3, "NO2": 46.006e-3, "CO": 28.010e-3}   # kg/mol (IUPAC atomic weights)

    def __init__(self, tspec: dict, cspec: dict, tables, eq: EquilibriumTable | None = None):
        self.tspec, self.cspec, self.tables = tspec, cspec, tables
        self.eq = eq or EquilibriumTable(tspec, cspec, tables)
        self.kp = zeldovich_params(cspec)
        em = cspec["emissions"]
        self.n_el = int(em["n_elements"])
        self.tol, self.maxit = float(em["newton_tol"]), int(em["newton_max_iter"])
        self.nox_as = em["nox_as"]
        ix = self.eq.ix
        self.idx = (ix["N2"], ix["O"], ix["OH"], ix["H"], ix["NO"])

    def closed(self, cl, rpm: float) -> dict:
        """NO (frozen at EVO) and equilibrium CO of a closed two-zone run, per cylinder and cycle."""
        num = self.tspec["numerics"]
        omega = rpm * 2.0 * math.pi / SEC_PER_MIN
        targs = self.tables.args()
        p, Tu, Tb, xb = cl.col(tz.O_P), cl.col(tz.O_TU), cl.col(tz.O_TB), cl.col(tz.O_XB)
        y, formed, Tpk = _elements(cl.theta, p, Tu, Tb, xb, cl.xr, omega, self.n_el, *targs, *self.eq.args(),
                                   *self.idx, self.kp, self.tol, self.maxit, float(num["T_flame_guess"]),
                                   self.MW["NO"])
        m_el = cl.m * xb[-1] / self.n_el
        m_no = float(np.sum(y[formed]) * m_el)
        mb, mu = cl.m * xb[-1], cl.m * (1.0 - xb[-1])
        out = np.empty(self.eq.L.shape[0])
        _interp(Tb[-1], p[-1], *self.eq.args(), out)
        mw_b = out[-1]
        n_b = mb / mw_b
        n_u = mu * frozen(Tu[-1], cl.xr, *targs[:5])[2] / R_MOL
        n_tot = n_b + n_u
        X_co_evo = math.exp(out[self.eq.ix["CO"]])
        ipk = int(np.nanargmax(Tb))
        _interp(Tb[ipk], p[ipk], *self.eq.args(), out)
        X_co_pk = math.exp(out[self.eq.ix["CO"]])
        n_no = m_no / self.MW["NO"]
        return {"m_NO_kg": m_no, "m_NOx_kg": n_no * self.MW[self.nox_as],
                "NO_ppm": n_no / n_tot * 1.0e6,
                "m_CO_kg": X_co_evo * n_b * self.MW["CO"], "CO_pct_evo": X_co_evo * n_b / n_tot * 100.0,
                "CO_pct_peak": X_co_pk * n_b / n_tot * 100.0, "m_CO_peak_kg": X_co_pk * n_b * self.MW["CO"],
                "T_element_max_K": float(np.max(Tpk)), "T_b_max_K": float(Tb[ipk])}


def specific(mass_per_cycle_kg: float, cylinders: int, rpm: float, brake_power_w: float) -> float:
    """Brake-specific emission [g/kWh] from a per-cylinder, per-cycle mass (four-stroke)."""
    if brake_power_w <= 0.0:
        return float("nan")
    rate = mass_per_cycle_kg * cylinders * rpm / SEC_PER_MIN / REV_PER_CYCLE
    return rate / KG_PER_G / (brake_power_w / J_PER_KWH)
