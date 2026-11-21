"""Piston cooling jets (M6, ADR 0006): jet state from the nozzle hydraulics, the Stevens & Webb (1991)
impingement heat-transfer coefficient on the crown underside, and the oil-pump power the jets cost.

Hydraulics (spec/thermo.yaml `piston_cooling_jets`):
    Q = Cd A sqrt(2 dp / rho)      per orifice,  u = Q / A  (mean nozzle-exit velocity, below the
                                   vena-contracta velocity, so conservative for heat transfer)
    Re = u d / nu,  Pr = mu cp / k,  mu from spec/dynamics.yaml bearings.oil_viscosity (single source)

Impingement (free, single-phase, round liquid jet on a flat plate):
    J. Stevens and B. W. Webb, "Local heat transfer coefficients under an axisymmetric, single-phase liquid
    jet", J. Heat Transfer 113(1), 71-78 (1991), doi:10.1115/1.2910554:
        Nu_0 = 2.67 Re^0.567 Pr^0.4 (z/d)^-0.0336 (u/d)^-0.237     (u/d in 1/s; Re 4,000 - 52,000)
        Nu / Nu_0 = [1 + f(r/d)^-9]^(-1/9),  f = a exp(-b r/d)
    as restated for piston oil-jet cooling by M. B. Varghese and A. K. Agarwal, "Numerical investigations of
    piston cooling using oil jet", SAE 2004-28-0061 (eq. 3-4, Table 1: a = 1.13, b = 0.23 for d = 2.2 mm)
    and in the same form by Varghese, Goyal and Agarwal, SAE 2005-01-1382. Nu = h d / k_oil. The jet
    velocity is the cycle-averaged relative velocity (Varghese & Agarwal): the piston velocity averages to
    zero over the cycle. The restatement prints the radial function without the minus sign; only the
    decaying form gives Nu <= Nu_0 away from the stagnation point (ASSUMED reading, see spec/fea.yaml).

Pump cost: the jets are extra delivery of a variable-displacement pump at the regulated supply pressure,
    P_shaft = n_jets Q dp / eta_pump,   FMEP_jets = P_shaft / (V_d,total N / 120)   (four-stroke)
when the jets flow (check valve open above `on_above_rpm`), zero below. All values are MODELLED.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from redline.spec import ROOT, geometry, load

THERMO_SPEC = ROOT / "spec" / "thermo.yaml"
DYN_SPEC = ROOT / "spec" / "dynamics.yaml"
SECONDS_PER_MIN = 60.0
REV_PER_CYCLE = 2.0                 # four-stroke
RPM_REL_TOL = 1.0e-9                # round-off allowance when the speed comes back from a mean piston speed


def jet_spec() -> dict:
    return load(str(THERMO_SPEC))["piston_cooling_jets"]


@dataclass(frozen=True)
class JetState:
    d: float                        # m, orifice diameter
    u: float                        # m/s, mean exit velocity
    q_orifice: float                # m^3/s per orifice
    n_jets: int                     # orifices in the engine
    rho: float
    mu: float
    cp: float
    k: float
    re: float
    pr: float
    supply_pressure: float          # Pa gauge

    @property
    def q_total(self) -> float:
        return self.q_orifice * self.n_jets


@lru_cache(maxsize=None)
def jet_state() -> JetState:
    js = jet_spec()
    oil = js["oil"]
    rho, cp, k = float(oil["rho"]), float(oil["cp"]), float(oil["k"])
    mu = float(load(str(DYN_SPEC))["bearings"]["oil_viscosity"])
    d = float(js["orifice_diameter"])
    dp = float(js["supply_pressure"])
    u = float(js["discharge_coefficient"]) * math.sqrt(2.0 * dp / rho)
    q = u * math.pi * d**2 / 4.0
    n = int(js["orifices_per_cylinder"]) * geometry().cylinders
    return JetState(d, u, q, n, rho, mu, cp, k, rho * u * d / mu, mu * cp / k, dp)


def stagnation_nusselt(st: JetState, z: float, corr: dict) -> float:
    """Stevens & Webb stagnation Nusselt number (coefficients from spec/fea.yaml, signs explicit)."""
    lo, hi = (float(v) for v in corr["re_range"])
    if not lo <= st.re <= hi:
        raise ValueError(f"jet Re {st.re:.0f} outside the Stevens & Webb range {lo:.0f}-{hi:.0f}")
    return (float(corr["C"]) * st.re ** float(corr["re_exp"]) * st.pr ** float(corr["pr_exp"])
            * (z / st.d) ** float(corr["zd_exp"]) * (st.u / st.d) ** float(corr["ud_exp"]))


def local_h(r, st: JetState, z: float, corr: dict):
    """Local impingement h [W/(m^2 K)] at radius r [m] (scalar or numpy array) from the jet axis."""
    import numpy as np
    h0 = stagnation_nusselt(st, z, corr) * st.k / st.d
    rad = corr["radial"]
    f = float(rad["a"]) * np.exp(float(rad["b_exp"]) * np.asarray(r, dtype=float) / st.d)
    n = float(rad["n"])
    return h0 * (1.0 + f ** (-n)) ** (-1.0 / n)


def jets_flowing(rpm: float) -> bool:
    js = jet_spec()
    return bool(js["enabled"]) and rpm >= float(js["on_above_rpm"]) * (1.0 - RPM_REL_TOL)


def pump_power(rpm: float) -> float:
    """Oil-pump shaft power [W] spent on the jets (whole engine)."""
    if not jets_flowing(rpm):
        return 0.0
    st = jet_state()
    return st.q_total * st.supply_pressure / float(jet_spec()["pump_efficiency"])


def fmep_jets(rpm: float) -> float:
    """FMEP [Pa] equivalent of the jet pump power at engine speed rpm."""
    if rpm <= 0.0:
        return 0.0
    vd = geometry().displacement
    return pump_power(rpm) / (vd * rpm / (SECONDS_PER_MIN * REV_PER_CYCLE))


def fmep_jets_from_piston_speed(mean_piston_speed: float) -> float:
    """Same, from the mean piston speed S_p = 2 s N / 60 used by the friction correlation."""
    rpm = mean_piston_speed * SECONDS_PER_MIN / (2.0 * geometry().stroke)
    return fmep_jets(rpm)
