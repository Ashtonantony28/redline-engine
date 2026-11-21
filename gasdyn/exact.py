"""Closed-form references for the M4 analytic tests (CLAUDE.md rule 2: analytic before numeric)."""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------------------------------------
# exact Riemann solution of the 1D Euler equations (Toro 2009, ch. 4)
# ---------------------------------------------------------------------------------------------
def _fk(p, rk, pk, g):
    ck = math.sqrt(g * pk / rk)
    if p > pk:   # shock
        a, b = 2.0 / ((g + 1.0) * rk), (g - 1.0) / (g + 1.0) * pk
        return (p - pk) * math.sqrt(a / (p + b))
    return 2.0 * ck / (g - 1.0) * ((p / pk) ** ((g - 1.0) / (2.0 * g)) - 1.0)


def riemann_exact(left, right, gamma, xi):
    """Self-similar solution at xi = x / t for left / right states (rho, u, p). Returns rho, u, p arrays."""
    rl, ul, pl = left
    rr, ur, pr = right
    g = gamma
    cl, cr = math.sqrt(g * pl / rl), math.sqrt(g * pr / rr)
    f = lambda p: _fk(p, rl, pl, g) + _fk(p, rr, pr, g) + ur - ul
    ps = brentq(f, 1e-12 * min(pl, pr), 1e3 * max(pl, pr), xtol=1e-14, rtol=1e-14)
    us = 0.5 * (ul + ur) + 0.5 * (_fk(ps, rr, pr, g) - _fk(ps, rl, pl, g))
    gm, gp = g - 1.0, g + 1.0
    xi = np.asarray(xi, dtype=float)
    rho, u, p = np.empty_like(xi), np.empty_like(xi), np.empty_like(xi)
    for k, s in enumerate(xi.flat):
        if s <= us:   # left of the contact
            if ps > pl:
                sl = ul - cl * math.sqrt(gp / (2 * g) * ps / pl + gm / (2 * g))
                w = (rl, ul, pl) if s < sl else (rl * (ps / pl + gm / gp) / (gm / gp * ps / pl + 1.0), us, ps)
            else:
                csl = cl * (ps / pl) ** (gm / (2 * g))
                if s < ul - cl:
                    w = (rl, ul, pl)
                elif s > us - csl:
                    w = (rl * (ps / pl) ** (1.0 / g), us, ps)
                else:
                    c = 2.0 / gp * (cl + gm / 2.0 * (ul - s))
                    w = (rl * (c / cl) ** (2.0 / gm), 2.0 / gp * (cl + gm / 2.0 * ul + s), pl * (c / cl) ** (2 * g / gm))
        else:
            if ps > pr:
                sr = ur + cr * math.sqrt(gp / (2 * g) * ps / pr + gm / (2 * g))
                w = (rr, ur, pr) if s > sr else (rr * (ps / pr + gm / gp) / (gm / gp * ps / pr + 1.0), us, ps)
            else:
                csr = cr * (ps / pr) ** (gm / (2 * g))
                if s > ur + cr:
                    w = (rr, ur, pr)
                elif s < us + csr:
                    w = (rr * (ps / pr) ** (1.0 / g), us, ps)
                else:
                    c = 2.0 / gp * (cr - gm / 2.0 * (ur - s))
                    w = (rr * (c / cr) ** (2.0 / gm), 2.0 / gp * (-cr + gm / 2.0 * ur + s), pr * (c / cr) ** (2 * g / gm))
        rho.flat[k], u.flat[k], p.flat[k] = w
    return rho, u, p


# ---------------------------------------------------------------------------------------------
# steady isentropic orifice / nozzle flow (Heywood App. C, eq. C.8 and C.9)
# ---------------------------------------------------------------------------------------------
def critical_pressure_ratio(gamma):
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


def orifice_mass_flow(area, p0, T0, p_back, gamma, R):
    """Isentropic mass flow [kg/s] through effective area `area` from stagnation (p0, T0) to p_back."""
    pr = max(p_back / p0, critical_pressure_ratio(gamma))
    if pr >= 1.0:
        return 0.0
    return area * p0 / math.sqrt(R * T0) * math.sqrt(
        2.0 * gamma / (gamma - 1.0) * (pr ** (2.0 / gamma) - pr ** ((gamma + 1.0) / gamma)))


# ---------------------------------------------------------------------------------------------
# acoustics and intake tuning
# ---------------------------------------------------------------------------------------------
def quarter_wave_frequency(c, length):
    """Fundamental of a pipe closed at one end and open at the other [Hz] (no end correction)."""
    return c / (4.0 * length)


def helmholtz_frequency(c, area, length, volume):
    """Helmholtz resonator [Hz]: neck of area and length, cavity volume."""
    return c / (2.0 * math.pi) * math.sqrt(area / (length * volume))


def runner_tuned_rpm(c, area, length, swept_volume, compression_ratio, ratio=2.0):
    """Engine speed of the induction ram peak for one cylinder on one runner (Engelman 1973, as given in
    Heywood sec. 6.6 / Stone, "Introduction to Internal Combustion Engines", ch. 6): the Helmholtz
    frequency of the runner on the mean cylinder volume V_eff = V_d (r_c + 1) / (2 (r_c - 1)) is about
    `ratio` (~2) times the piston frequency N / 60."""
    v_eff = swept_volume * (compression_ratio + 1.0) / (2.0 * (compression_ratio - 1.0))
    return 60.0 * helmholtz_frequency(c, area, length, v_eff) / ratio
