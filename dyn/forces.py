"""Gas and inertia forces, crankpin forces and crank torque (single and multi-cylinder).

Conventions (ADR 0001): local crank angle t of a cylinder = 0 at its firing TDC; piston
displacement s positive downward; forces along the bore positive DOWNWARD (towards the crank);
positive torque drives the crank forward. Conrod lumped as two masses (small end reciprocating,
big end rotating); the rod's own moment-of-inertia correction is neglected (documented limit).

Piston force   F_p = (p - p_cc) A_p - m_rec * d2s/dt2
Rod angle      beta = asin(lam sin t)
Rod force      F_rod = F_p / cos(beta)            (compression positive)
Side thrust    F_n = F_p tan(beta)                (on the liner, positive towards +y, the side the
                                                   crankpin swings to in 0 < t < pi)
Crankpin       F_t = F_rod sin(t + beta)          (tangential, positive in rotation direction)
               F_r = F_rod cos(t + beta)          (radial, positive towards the crank axis)
Torque         T = F_t r = F_p ds/dtheta
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kin.slider_crank import conrod_angle, cylinder_volume, piston_acceleration, piston_velocity
from dyn import params
from dyn.pressure import PressureTrace


@dataclass
class CylinderForces:
    theta: np.ndarray        # local crank angle [rad]
    F_gas: np.ndarray        # N, down
    F_inertia: np.ndarray    # N, down (d'Alembert of reciprocating mass)
    F_piston: np.ndarray     # N, down
    beta: np.ndarray         # rad
    F_rod: np.ndarray        # N, compression +
    F_side: np.ndarray       # N
    F_tangential: np.ndarray
    F_radial: np.ndarray
    T_gas: np.ndarray        # N m
    T_inertia: np.ndarray    # N m
    T: np.ndarray            # N m


def cylinder_forces(theta, trace: PressureTrace, omega: float, spec=None, dspec=None) -> CylinderForces:
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    g = params.geometry(spec)
    m = params.masses(spec)
    r, l = g.crank_radius, g.conrod_length
    th = np.asarray(theta, dtype=float)
    p_cc = dspec["layout"]["crankcase_pressure"]

    F_gas = (trace.at(th) - p_cc) * g.piston_area
    F_in = -m.reciprocating * piston_acceleration(th, omega, r, l)
    F_p = F_gas + F_in
    beta = conrod_angle(th, r, l)
    F_rod = F_p / np.cos(beta)
    dsdth = piston_velocity(th, 1.0, r, l)
    return CylinderForces(
        theta=th, F_gas=F_gas, F_inertia=F_in, F_piston=F_p, beta=beta, F_rod=F_rod,
        F_side=F_p * np.tan(beta),
        F_tangential=F_rod * np.sin(th + beta),
        F_radial=F_rod * np.cos(th + beta),
        T_gas=F_gas * dsdth, T_inertia=F_in * dsdth, T=F_p * dsdth,
    )


def multi_cylinder_torque(theta, trace: PressureTrace, omega: float, spec=None, dspec=None):
    """Crank torque of every cylinder on the global crank-angle grid (cylinder 1 local angle).

    Returns (T_cyl [n_cyl, N], T_total [N]). Cylinder k sees local angle theta - phi_k with
    phi_k its firing offset (firing order 1-3-4-2, ADR 0001). All cylinders share `trace`.
    """
    th = np.asarray(theta, dtype=float)
    phases = params.firing_offsets(spec)
    T_cyl = np.array([cylinder_forces(th - ph, trace, omega, spec, dspec).T for ph in phases])
    return T_cyl, T_cyl.sum(axis=0)


def imep(trace: PressureTrace, spec=None) -> float:
    """Net indicated mean effective pressure over the 720-deg cycle: (closed integral p dV) / V_s."""
    g = params.geometry(spec or params.engine_spec())
    V = cylinder_volume(trace.theta, g)
    return float(np.trapezoid(trace.p, V) / g.swept_volume)


def mean_torque_identity(trace: PressureTrace, spec=None) -> float:
    """Analytic mean crank torque of the whole engine: IMEP * V_d / (4 pi)."""
    g = params.geometry(spec or params.engine_spec())
    return imep(trace, spec) * g.displacement / params.FULL_CYCLE


def cycle_grid(n: int) -> np.ndarray:
    """Uniform grid over one cycle [0, 4 pi), endpoint excluded (for exact periodic means/FFTs)."""
    return np.linspace(0.0, params.FULL_CYCLE, n, endpoint=False)
