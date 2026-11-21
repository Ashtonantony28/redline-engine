"""Thermodynamic property tables from Cantera for the M2 two-zone model (PLAN.md sec. 3, risk "Cantera cost").

Three gas mixtures, all built from one Cantera ideal-gas phase:

* ``F``  fresh charge: iso-octane vapour + air at the spec lambda (frozen composition).
* ``P``  complete-combustion products of ``F`` (CO2, H2O, N2, excess O2), frozen. Used for the residual
         gas in the unburned zone and for burned gas during gas exchange.
* ``B``  burned gas in chemical equilibrium at (T, p) with the elemental composition of ``F``.

Species thermo: combustion products from Cantera's bundled ``gri30.yaml`` (GRI-Mech 3.0 NASA-7 fits);
iso-octane (2,2,4-trimethylpentane) vapour from Cantera's bundled ``nasa_gas.yaml`` entry
``C8H18,isooctane`` (NASA Glenn thermo database, McBride, Zehe & Gordon, NASA/TP-2002-211556, "X 4/85").

The unburned zone is a frozen mixture of F and P, so its properties are exactly mass-fraction-weighted
blends of the 1-D F and P tables. The burned table is 2-D in (T, ln p) and stores h, R and their partial
derivatives (central differences on the grid) for the zone-temperature ODEs. All interpolation is
(bi)linear and runs inside numba kernels.

Enthalpies are absolute (formation included), so chemical energy release is implicit in h_u - h_b.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import cantera as ct
import numpy as np
from numba import njit

from redline.spec import ROOT, load

THERMO_SPEC = ROOT / "spec" / "thermo.yaml"
CACHE_DIR = ROOT / "out" / "thermo" / "cache"
_TABLE_VERSION = "m2-props-2"   # bump to invalidate cached tables when this module's maths changes

# burned-table layers
B_H, B_R, B_HT, B_HLP, B_RT, B_RLP = range(6)
N_B_LAYERS = 6
# frozen-table layers
F_H, F_CP, P_H, P_CP = range(4)
_KMOL_PER_MOL = 1.0e3   # Cantera molar quantities are per kmol


def thermo_spec(path: str | None = None) -> dict:
    return load(str(path or THERMO_SPEC))


# ---------------------------------------------------------------------------------------------
# Cantera phase and mixture compositions
# ---------------------------------------------------------------------------------------------
def make_gas(tspec: dict) -> ct.Solution:
    mix = tspec["mixture"]
    gri = {s.name: s for s in ct.Species.list_from_file("gri30.yaml")}
    fuel = [s for s in ct.Species.list_from_file("nasa_gas.yaml") if s.name == mix["fuel_species"]]
    if len(fuel) != 1:
        raise KeyError(f"fuel species {mix['fuel_species']!r} not found in nasa_gas.yaml")
    species = [gri[n] for n in tspec["tables"]["product_species"]] + fuel
    return ct.Solution(thermo="ideal-gas", species=species)


@dataclass(frozen=True)
class Mixture:
    fuel: str
    X_fresh: dict
    X_products: dict
    lam: float
    o2_stoich: float          # mol O2 per mol fuel at stoichiometry
    afr_stoich: float         # kg air / kg fuel from the species data
    fuel_mass_fraction: float  # in the fresh charge


def mixture(tspec: dict, gas: ct.Solution | None = None) -> Mixture:
    gas = gas or make_gas(tspec)
    mix = tspec["mixture"]
    lam = float(mix["lambda"])
    if lam < 1.0:
        raise NotImplementedError("M2 handles lean/stoichiometric mixtures only (complete products undefined)")
    fuel = mix["fuel_species"]
    nC = gas.n_atoms(fuel, "C")
    nH = gas.n_atoms(fuel, "H")
    nO = gas.n_atoms(fuel, "O")
    a_st = nC + nH / 4.0 - nO / 2.0
    rN = float(mix["o2_n2_ratio"])
    a = lam * a_st
    X_fresh = {fuel: 1.0, "O2": a, "N2": a * rN}
    X_prod = {"CO2": nC, "H2O": nH / 2.0, "N2": a * rN}
    if a > a_st:
        X_prod["O2"] = a - a_st
    W = {n: gas.molecular_weights[gas.species_index(n)] for n in (fuel, "O2", "N2")}
    afr_st = a_st * (W["O2"] + rN * W["N2"]) / W[fuel]
    yf = W[fuel] / (W[fuel] + a * (W["O2"] + rN * W["N2"]))
    return Mixture(fuel, X_fresh, X_prod, lam, a_st, afr_st, yf)


# ---------------------------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------------------------
@dataclass
class Tables:
    T0f: float
    dTf: float
    Fz: np.ndarray        # (4, n): hF, cpF, hP, cpP  [J/kg, J/kg/K]
    R_F: float
    R_P: float
    T0b: float
    dTb: float
    lp0: float
    dlp: float
    B: np.ndarray         # (6, nT, nP): h, R, dh/dT, dh/dlnp, dR/dT, dR/dlnp
    Y_F: np.ndarray       # mass fractions (for direct Cantera checks)
    Y_P: np.ndarray
    lhv_vapour: float     # J/kg fuel, from the species data at T_ref (water as vapour)
    lhv_cv_vapour: float  # constant-volume heating value at T_ref
    fuel_mass_fraction: float
    afr_stoich: float
    R_air: float
    meta: dict

    @property
    def Tf_grid(self):
        return self.T0f + self.dTf * np.arange(self.Fz.shape[1])

    @property
    def Tb_grid(self):
        return self.T0b + self.dTb * np.arange(self.B.shape[1])

    @property
    def lp_grid(self):
        return self.lp0 + self.dlp * np.arange(self.B.shape[2])

    def args(self):
        """Positional tuple consumed by the numba kernels."""
        return (self.Fz, self.T0f, self.dTf, self.R_F, self.R_P, self.B, self.T0b, self.dTb, self.lp0, self.dlp)


def _key(tspec: dict) -> str:
    blob = json.dumps({"mix": tspec["mixture"], "tab": tspec["tables"], "v": _TABLE_VERSION,
                       "ct": ct.__version__}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def _frozen_props(gas, Y, T):
    gas.TPY = T, gas.P, Y
    return gas.enthalpy_mass, gas.cp_mass


def build_tables(tspec: dict | None = None, use_cache: bool = True) -> Tables:
    tspec = tspec or thermo_spec()
    key = _key(tspec)
    cache = CACHE_DIR / f"props_{key}.npz"
    if use_cache and cache.exists():
        d = np.load(cache, allow_pickle=False)
        meta = json.loads(str(d["meta"]))
        return Tables(B=d["B"], Fz=d["Fz"], Y_F=d["Y_F"], Y_P=d["Y_P"], meta=meta, **meta["scalars"])

    gas = make_gas(tspec)
    mx = mixture(tspec, gas)
    tab = tspec["tables"]
    Tref = float(tspec["mixture"]["T_ref"])
    pref = float(tspec["mixture"]["p_ref"])

    gas.TPX = Tref, pref, mx.X_fresh
    Y_F = gas.Y.copy()
    R_F = ct.gas_constant / gas.mean_molecular_weight
    hF_ref, uF_ref = gas.enthalpy_mass, gas.int_energy_mass
    gas.TPX = Tref, pref, mx.X_products
    Y_P = gas.Y.copy()
    R_P = ct.gas_constant / gas.mean_molecular_weight
    hP_ref, uP_ref = gas.enthalpy_mass, gas.int_energy_mass
    yf = Y_F[gas.species_index(mx.fuel)]
    lhv = (hF_ref - hP_ref) / yf
    lhv_cv = (uF_ref - uP_ref) / yf

    T0f, T1f, dTf = float(tab["T_min"]), float(tab["T_max"]), float(tab["dT_frozen"])
    Tf = np.arange(T0f, T1f + 0.5 * dTf, dTf)
    Fz = np.empty((4, Tf.size))
    for i, T in enumerate(Tf):
        Fz[F_H, i], Fz[F_CP, i] = _frozen_props(gas, Y_F, T)
        Fz[P_H, i], Fz[P_CP, i] = _frozen_props(gas, Y_P, T)

    T0b, T1b, dTb = float(tab["Tb_min"]), float(tab["Tb_max"]), float(tab["dT_burned"])
    Tb = np.arange(T0b, T1b + 0.5 * dTb, dTb)
    lp = np.linspace(np.log(float(tab["p_min"])), np.log(float(tab["p_max"])), int(tab["n_p"]))
    dlp = lp[1] - lp[0]
    H = np.empty((Tb.size, lp.size))
    R = np.empty_like(H)
    for j, l in enumerate(lp):
        gas.TPY = Tb[0], np.exp(l), Y_P
        for i, T in enumerate(Tb):   # warm start: previous equilibrium composition as the initial guess
            gas.TP = T, np.exp(l)
            gas.equilibrate("TP")
            H[i, j] = gas.enthalpy_mass
            R[i, j] = ct.gas_constant / gas.mean_molecular_weight
    B = np.empty((N_B_LAYERS, Tb.size, lp.size))
    B[B_H], B[B_R] = H, R
    B[B_HT] = np.gradient(H, dTb, axis=0)
    B[B_HLP] = np.gradient(H, dlp, axis=1)
    B[B_RT] = np.gradient(R, dTb, axis=0)
    B[B_RLP] = np.gradient(R, dlp, axis=1)

    W_air = float(tspec["mixture"]["air_molar_mass"])
    scalars = dict(T0f=T0f, dTf=dTf, R_F=R_F, R_P=R_P, T0b=T0b, dTb=dTb, lp0=float(lp[0]), dlp=float(dlp),
                   lhv_vapour=lhv, lhv_cv_vapour=lhv_cv, fuel_mass_fraction=float(yf),
                   afr_stoich=mx.afr_stoich, R_air=ct.gas_constant / (W_air * _KMOL_PER_MOL))
    meta = {"scalars": scalars, "key": key, "fuel": mx.fuel}
    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(cache, B=B, Fz=Fz, Y_F=Y_F, Y_P=Y_P, meta=json.dumps(meta))
    return Tables(B=B, Fz=Fz, Y_F=Y_F, Y_P=Y_P, meta=meta, **scalars)


def ideal_tables(gamma: float, R: float, q_mix: float, T0: float = 200.0, T1: float = 6000.0,
                 dT: float = 10.0, p0: float = 1.0e4, p1: float = 1.0e9) -> Tables:
    """Constant-cp ideal-gas tables (same cp and R for F, P and B) for the analytic Otto anchor.

    ``q_mix`` is the chemical energy released per kg of mixture (h_F - h_P at any T). With identical cp
    and R in both zones the two-zone equations collapse exactly to a single-zone model.
    """
    cp = gamma * R / (gamma - 1.0)
    T = np.arange(T0, T1 + 0.5 * dT, dT)
    Fz = np.empty((4, T.size))
    Fz[F_H], Fz[F_CP] = cp * T + q_mix, cp
    Fz[P_H], Fz[P_CP] = cp * T, cp
    lp = np.linspace(np.log(p0), np.log(p1), 2)
    B = np.zeros((N_B_LAYERS, T.size, lp.size))
    B[B_H] = (cp * T)[:, None]
    B[B_R] = R
    B[B_HT] = cp
    nan = float("nan")
    return Tables(T0f=T0, dTf=dT, Fz=Fz, R_F=R, R_P=R, T0b=T0, dTb=dT, lp0=float(lp[0]), dlp=float(lp[1] - lp[0]),
                  B=B, Y_F=np.empty(0), Y_P=np.empty(0), lhv_vapour=nan, lhv_cv_vapour=nan,
                  fuel_mass_fraction=nan, afr_stoich=nan, R_air=R, meta={"ideal": True})


# ---------------------------------------------------------------------------------------------
# Direct Cantera evaluations (reference for the table tests)
# ---------------------------------------------------------------------------------------------
class Direct:
    """Direct Cantera property calls, used only to check the tables."""

    def __init__(self, tspec: dict | None = None, tables: Tables | None = None):
        self.tspec = tspec or thermo_spec()
        self.gas = make_gas(self.tspec)
        self.t = tables or build_tables(self.tspec)

    def unburned(self, T: float, p: float, xr: float):
        Y = (1.0 - xr) * self.t.Y_F + xr * self.t.Y_P
        self.gas.TPY = T, p, Y
        return self.gas.enthalpy_mass, self.gas.cp_mass, ct.gas_constant / self.gas.mean_molecular_weight

    def burned(self, T: float, p: float):
        self.gas.TPY = T, p, self.t.Y_P
        self.gas.equilibrate("TP")
        return self.gas.enthalpy_mass, ct.gas_constant / self.gas.mean_molecular_weight


# ---------------------------------------------------------------------------------------------
# numba interpolation kernels
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def _cell(x, x0, dx, n):
    s = (x - x0) / dx
    i = int(np.floor(s))
    if i < 0:
        i = 0
    elif i > n - 2:
        i = n - 2
    return i, s - i


@njit(cache=True)
def frozen(T, xr, Fz, T0f, dTf, R_F, R_P):
    """Unburned (fresh + residual products, frozen) h, cp, R at temperature T."""
    i, f = _cell(T, T0f, dTf, Fz.shape[1])
    hF = Fz[0, i] + f * (Fz[0, i + 1] - Fz[0, i])
    cF = Fz[1, i] + f * (Fz[1, i + 1] - Fz[1, i])
    hP = Fz[2, i] + f * (Fz[2, i + 1] - Fz[2, i])
    cP = Fz[3, i] + f * (Fz[3, i + 1] - Fz[3, i])
    a = 1.0 - xr
    return a * hF + xr * hP, a * cF + xr * cP, a * R_F + xr * R_P


@njit(cache=True)
def burned(T, p, B, T0b, dTb, lp0, dlp):
    """Equilibrium burned gas: h, R, dh/dT, dh/dp, dR/dT, dR/dp at (T, p)."""
    i, f = _cell(T, T0b, dTb, B.shape[1])
    j, g = _cell(np.log(p), lp0, dlp, B.shape[2])
    w00 = (1 - f) * (1 - g)
    w10 = f * (1 - g)
    w01 = (1 - f) * g
    w11 = f * g
    out = np.empty(6)
    for k in range(6):
        out[k] = w00 * B[k, i, j] + w10 * B[k, i + 1, j] + w01 * B[k, i, j + 1] + w11 * B[k, i + 1, j + 1]
    return out[0], out[1], out[2], out[3] / p, out[4], out[5] / p


@njit(cache=True)
def T_from_u_frozen(u, xr, Tguess, Fz, T0f, dTf, R_F, R_P, tol, maxit):
    """Invert u(T) = h - R T for the frozen mixture (Newton)."""
    T = Tguess
    for _ in range(maxit):
        h, cp, R = frozen(T, xr, Fz, T0f, dTf, R_F, R_P)
        r = h - R * T - u
        dT = r / (cp - R)
        T -= dT
        if abs(dT) < tol * T:
            break
    return T


@njit(cache=True)
def T_adiabatic_p(h_target, p, Tguess, B, T0b, dTb, lp0, dlp, tol, maxit):
    """Burned-gas temperature with h_b(T, p) = h_target (constant-pressure adiabatic flame)."""
    T = Tguess
    for _ in range(maxit):
        h, R, hT, hp, RT, Rp = burned(T, p, B, T0b, dTb, lp0, dlp)
        dT = (h - h_target) / hT
        T -= dT
        if abs(dT) < tol * T:
            break
    return T


@njit(cache=True)
def burn_constant_volume(u_target, m, V, Tguess, B, T0b, dTb, lp0, dlp, tol, maxit):
    """Adiabatic constant-volume combustion: find (T, p) with u_b(T, p) = u_target, p = m R_b T / V."""
    T = Tguess
    p = np.exp(lp0)   # any in-table pressure: R_b depends only weakly on p, refined below
    for _ in range(maxit):
        h, R, hT, hp, RT, Rp = burned(T, p, B, T0b, dTb, lp0, dlp)
        p = m * R * T / V
        h, R, hT, hp, RT, Rp = burned(T, p, B, T0b, dTb, lp0, dlp)
        r = h - R * T - u_target
        dT = r / (hT - R - T * RT)
        T -= dT
        if abs(dT) < tol * T:
            break
    h, R, hT, hp, RT, Rp = burned(T, p, B, T0b, dTb, lp0, dlp)
    p = m * R * T / V
    return T, p
