"""SYNTHETIC compressor and turbine maps from documented non-dimensional models (M10 turbo stretch).

Openly published maps of a turbocharger sized for this engine were not available, so both maps are
generated from physically consistent non-dimensional models whose constants live in spec/turbo.yaml.
They are clearly SYNTHETIC; all numbers derived from them are MODELLED (CLAUDE.md rule 5).

Compressor (centrifugal, backswept impeller)
--------------------------------------------
    U2 = omega D2 / 2,  Mu = U2 / a01,  phi = m / (rho01 U2 D2^2)
    Euler work:      dh0 = lambda U2^2,  lambda = lambda_0 - lambda_slope phi    (Dixon & Hall, "Fluid
                     Mechanics and Thermodynamics of Turbomachinery", 7th ed., ch. 7: backsweep makes the
                     work coefficient fall linearly with flow)
    choke:           phi_ch(Mu) = phi_choke_0 / sqrt(1 + (Mu / mach_choke)^2)   (corrected choke flow rises
                     with speed and saturates once the inducer throat chokes: m_ch -> const rho01 a01 D2^2)
    efficiency:      eta = eta_peak (1 - k_M (Mu - Mu_opt)^2) (1 - ((x - x_pk) / x_w)^2) (1 - x^n),  x = phi/phi_ch,
                     x_w = surge-side or choke-side half-width (asymmetric islands)
    pressure ratio:  PR = (1 + eta dh0 / (cp T01))^(g / (g - 1))
    surge:           the peak of each speed line (dPR/dphi = 0), the static-instability (Greitzer) limit.
Compressor power is the Euler work m lambda U2^2, which by construction equals the closed form
m cp T01 (PR^((g-1)/g) - 1) / eta. Beyond choke (PR < 1, including the stopped wheel) the compressor is a
passive restriction of effective area `bypass_area` in parallel with the choked flow. Left of surge the
flow falls linearly to zero over `surge_band_pr` and the point is flagged (deep surge is not modelled).

Turbine (radial, nozzle-equivalent)
-----------------------------------
    flow:        isentropic nozzle of area A_eq from (p03, T03) to p4 (the gasdyn kernel's orifice boundary)
    efficiency:  eta_ts = eta_max (1 - ((BSR - BSR_opt) / BSR_opt)^2),  BSR = U_T / c_is,
                 c_is = sqrt(2 cp T03 (1 - (p4 / p03)^((g-1)/g)))     (Watson & Janota 1982, ch. 4-5)
    power:       P_T = m cp T03 eta_ts (1 - (p4 / p03)^((g-1)/g))
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from thermo.gasexchange import orifice

NO_LINEARISATION = 1.0      # thermo.gasexchange.orifice: pure isentropic orifice law


class CompressorMap:
    def __init__(self, cfg: dict, gamma: float, R: float):
        self.cfg = cfg
        self.gamma, self.R = float(gamma), float(R)
        self.cp = self.gamma * self.R / (self.gamma - 1.0)
        self.k = (self.gamma - 1.0) / self.gamma
        c = cfg
        self.D2 = float(c["impeller_diameter"])
        self.lam0, self.lam_s = float(c["lambda_0"]), float(c["lambda_slope"])
        self.eta_pk, self.mu_opt, self.k_mu = float(c["eta_peak"]), float(c["mach_opt"]), float(c["eta_mach_coeff"])
        self.phi_ch0, self.mu_ch = float(c["phi_choke_0"]), float(c["mach_choke"])
        self.x_pk, self.n_ch = float(c["x_peak"]), float(c["choke_exp"])
        self.x_wl, self.x_wr = float(c["x_width_surge_side"]), float(c["x_width_choke_side"])
        self.surge_band = float(c["surge_band_pr"])
        self.A_bypass = float(c["bypass_area"])
        self.min_pr = float(c["min_pr"])
        # surge position x_s = phi_s / phi_ch on a tip-Mach grid (the line shape depends on Mu via lambda)
        self._mu_tab = np.linspace(0.0, float(c["mach_table_max"]), int(c["mach_table_n"]))
        self._xs_tab = np.array([self._surge_x(mu) for mu in self._mu_tab])

    # ---- non-dimensional speed line --------------------------------------------------------------
    def phi_choke(self, mu):
        return self.phi_ch0 / np.sqrt(1.0 + (np.asarray(mu, dtype=float) / self.mu_ch) ** 2)

    def eta_max(self, mu):
        return self.eta_pk * (1.0 - self.k_mu * (np.asarray(mu, dtype=float) - self.mu_opt) ** 2)

    def eta_x(self, x, mu):
        x = np.asarray(x, dtype=float)
        w = np.where(x < self.x_pk, self.x_wl, self.x_wr)
        return self.eta_max(mu) * (1.0 - ((x - self.x_pk) / w) ** 2) * (1.0 - x ** self.n_ch)

    def lam(self, phi):
        return self.lam0 - self.lam_s * np.asarray(phi, dtype=float)

    def _work_term(self, x, mu):
        """eta * lambda along a speed line (PR is monotone in it at fixed speed)."""
        return self.eta_x(x, mu) * self.lam(x * self.phi_choke(mu))

    def _surge_x(self, mu):
        xs = np.linspace(1.0e-3, 1.0, 2001)
        f = self._work_term(xs, mu)
        i = int(np.argmax(f))
        lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
        gr = (math.sqrt(5.0) - 1.0) / 2.0
        for _ in range(60):
            a = hi - gr * (hi - lo)
            b = lo + gr * (hi - lo)
            if self._work_term(a, mu) > self._work_term(b, mu):
                hi = b
            else:
                lo = a
        return 0.5 * (lo + hi)

    def surge_x(self, mu):
        return np.interp(mu, self._mu_tab, self._xs_tab)

    # ---- dimensional helpers ---------------------------------------------------------------------
    def tip(self, omega, T01):
        U2 = omega * self.D2 / 2.0
        return U2, U2 / math.sqrt(self.gamma * self.R * T01)

    def flow_scale(self, omega, p01, T01):
        """rho01 U2 D2^2 [kg/s]: mass flow per unit flow coefficient."""
        U2, _ = self.tip(omega, T01)
        return p01 / (self.R * T01) * U2 * self.D2 ** 2

    def line(self, omega, p01, T01, m):
        """Speed-line point at mass flow m inside [surge, choke]: dict(pr, eta, lam, x, phi)."""
        U2, mu = self.tip(omega, T01)
        phi = m / self.flow_scale(omega, p01, T01)
        x = phi / float(self.phi_choke(mu))
        eta = float(self.eta_x(x, mu))
        lam = float(self.lam(phi))
        pr = (1.0 + eta * lam * U2 ** 2 / (self.cp * T01)) ** (1.0 / self.k)
        return {"pr": pr, "eta": eta, "lam": lam, "x": x, "phi": phi, "U2": U2, "mu": mu}

    def limits(self, omega, p01, T01):
        """(m_surge, pr_surge, m_choke) of the speed line at omega."""
        U2, mu = self.tip(omega, T01)
        fs = self.flow_scale(omega, p01, T01)
        phc = float(self.phi_choke(mu))
        xs = float(self.surge_x(mu))
        m_s = xs * phc * fs
        pr_s = self.line(omega, p01, T01, m_s)["pr"] if m_s > 0.0 else 1.0
        return m_s, pr_s, phc * fs

    def pr_of_flow(self, omega, p01, T01, m):
        """Monotone (for m >= m_surge) total-pressure ratio of the extended characteristic."""
        m_s, pr_s, m_ch = self.limits(omega, p01, T01)
        if m <= m_s:
            return pr_s
        if m < m_ch:
            return self.line(omega, p01, T01, m)["pr"]
        rho1 = p01 / (self.R * T01)
        return 1.0 - (m - m_ch) ** 2 / (2.0 * rho1 * self.A_bypass ** 2 * p01)

    def power(self, omega, p01, T01, m):
        """Shaft power [W] absorbed at flow m (Euler work of the bladed part of the flow)."""
        m_s, pr_s, m_ch = self.limits(omega, p01, T01)
        U2, _ = self.tip(omega, T01)
        mb = min(max(m, 0.0), m_ch)
        if mb <= 0.0:
            return 0.0
        phi = mb / self.flow_scale(omega, p01, T01)
        return mb * float(self.lam(phi)) * U2 ** 2

    def operating_point(self, omega, p01, T01, demand) -> dict:
        """Mass flow where the characteristic meets the downstream demand PR(m) (increasing in m)."""
        m_s, pr_s, m_ch = self.limits(omega, p01, T01)
        rho1 = p01 / (self.R * T01)
        m_max = m_ch + self.A_bypass * math.sqrt(2.0 * rho1 * p01 * (1.0 - self.min_pr))
        g = lambda m: self.pr_of_flow(omega, p01, T01, m) - demand(m)
        surge = False
        if g(m_s) <= 0.0:
            # demand above the surge pressure ratio: left of surge (flagged)
            excess = demand(m_s) - pr_s
            m = m_s * max(0.0, 1.0 - excess / self.surge_band) if m_s > 0.0 else 0.0
            surge = m_s > 0.0
        elif g(m_max) >= 0.0:
            m = m_max
        else:
            m = brentq(g, m_s, m_max, xtol=1.0e-9 * max(m_max, 1.0e-6), rtol=1.0e-10)
        P = self.power(omega, p01, T01, m)
        pr = self.pr_of_flow(omega, p01, T01, max(m, m_s))
        T2 = T01 + (P / (m * self.cp) if m > 0.0 else 0.0)
        eta = (self.cp * T01 * (pr ** self.k - 1.0) * m / P) if (P > 0.0 and pr > 1.0) else 0.0
        return {"mdot": m, "pr": pr, "power": P, "T2": T2, "eta": eta, "m_surge": m_s, "pr_surge": pr_s,
                "m_choke": m_ch, "surge": surge,
                "surge_margin": (m / m_s - 1.0) if m_s > 0.0 else float("inf")}

    # ---- corrected-flow map for plots and tests ----------------------------------------------------
    def speed_line(self, omega, p01, T01, n: int = 60):
        m_s, _, m_ch = self.limits(omega, p01, T01)
        ms = np.linspace(m_s, m_ch * (1.0 - 1.0e-6), n)
        pts = [self.line(omega, p01, T01, m) for m in ms]
        return ms, np.array([p["pr"] for p in pts]), np.array([p["eta"] for p in pts])


class TurbineMap:
    def __init__(self, cfg: dict, gamma: float, R: float):
        self.gamma, self.R = float(gamma), float(R)
        self.cp = self.gamma * self.R / (self.gamma - 1.0)
        self.k = (self.gamma - 1.0) / self.gamma
        self.A = float(cfg["equivalent_area"])
        self.D = float(cfg["rotor_diameter"])
        self.eta_max, self.bsr_opt = float(cfg["eta_max"]), float(cfg["bsr_opt"])
        self.eta_floor = float(cfg["eta_floor"])

    def flow(self, p03, T03, p4):
        """Nozzle-equivalent mass flow [kg/s] (same law as the kernel's orifice boundary)."""
        return orifice(self.A, p03, T03, self.R, self.gamma, p4, NO_LINEARISATION)

    def efficiency(self, bsr):
        e = self.eta_max * (1.0 - ((np.asarray(bsr, dtype=float) - self.bsr_opt) / self.bsr_opt) ** 2)
        return np.maximum(e, self.eta_floor)

    def power(self, m, p03, T03, p4, omega) -> dict:
        if m <= 0.0 or p4 >= p03:
            return {"power": 0.0, "eta": 0.0, "bsr": float("nan"), "pr": 1.0}
        dh_is = self.cp * T03 * (1.0 - (p4 / p03) ** self.k)
        c_is = math.sqrt(2.0 * dh_is)
        bsr = omega * self.D / 2.0 / c_is
        eta = float(self.efficiency(bsr))
        return {"power": m * eta * dh_is, "eta": eta, "bsr": bsr, "pr": p03 / p4}
