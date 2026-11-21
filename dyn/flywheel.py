"""Flywheel sizing by the energy-fluctuation method.

Load torque is taken equal to the mean crank torque (steady speed). The fluctuation energy is
dE = max - min of E(theta) = integral (T - T_mean) dtheta, and the coefficient of cyclic speed
variation is C_s = (w_max - w_min) / w_mean ~= dE / (I w^2) for small C_s.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from dyn import params


def fluctuation_energy(theta, torque) -> float:
    th = np.asarray(theta, dtype=float)
    T = np.asarray(torque, dtype=float)
    # close the cycle so the periodic integral of (T - T_mean) returns to zero
    th_c = np.append(th, th[0] + params.FULL_CYCLE)
    T_c = np.append(T, T[0])
    t_mean = np.trapezoid(T_c, th_c) / params.FULL_CYCLE
    E = cumulative_trapezoid(T_c - t_mean, th_c, initial=0.0)
    return float(E.max() - E.min())


def rotating_inertia(flywheel_inertia: float | None = None, spec=None, dspec=None) -> float:
    """Crank + flywheel (+ mean equivalent of the conrod/piston masses if configured) [kg m^2]."""
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    ms = spec["masses"]
    I = ms["crank_inertia"] + (ms["flywheel_inertia"] if flywheel_inertia is None else flywheel_inertia)
    if dspec["flywheel"]["include_reciprocating_equivalent"]:
        g = params.geometry(spec)
        m = params.masses(spec)
        r2 = g.crank_radius**2
        I += spec["engine"]["cylinders"] * (m.rod_big_end * r2 + 0.5 * m.reciprocating * r2)
    return I


def cyclic_irregularity(dE: float, inertia: float, omega: float) -> float:
    return dE / (inertia * omega**2)


def required_flywheel_inertia(dE: float, omega: float, cs_target: float, spec=None, dspec=None) -> float:
    """Flywheel inertia that meets `cs_target` (other rotating inertia subtracted)."""
    other = rotating_inertia(0.0, spec, dspec)
    return dE / (cs_target * omega**2) - other


def simulate_speed(theta, torque, inertia: float, omega0: float, cycles: int = 3):
    """Direct check: integrate I w dw/dtheta = T(theta) - T_mean; returns (theta, w) of the last cycle."""
    th = np.asarray(theta, dtype=float)
    th_c = np.append(th, th[0] + params.FULL_CYCLE)
    T_c = np.append(np.asarray(torque, dtype=float), torque[0])
    t_mean = np.trapezoid(T_c, th_c) / params.FULL_CYCLE

    def rhs(t, y):
        Tq = np.interp(np.mod(t, params.FULL_CYCLE), th_c, T_c)
        return [(Tq - t_mean) / (inertia * y[0])]

    span = (0.0, cycles * params.FULL_CYCLE)
    t_eval = th + (cycles - 1) * params.FULL_CYCLE
    max_step = float(np.min(np.diff(th_c))) * 4
    sol = solve_ivp(rhs, span, [omega0], t_eval=t_eval, rtol=1e-9, atol=1e-9, max_step=max_step)
    return th, sol.y[0]
