"""Shaking forces and moments of the flat-plane I4, and Lanchester balance-shaft sizing.

Body axes (ADR 0001): x along the crank (front to rear), z up the bore axis, y completes a
right-handed set. The shaking force is the reaction of the moving masses on the structure,
F = -sum(m a); gas forces are internal and cancel. With s downward from TDC the vertical shaking
force of one cylinder is F_z = +m_rec d2s/dt2 (the piston pulls the engine up at TDC).
Moments about the engine centre: M = sum(r_k x F_k), so M_y = -sum(x_k F_z,k), M_z = sum(x_k F_y,k).

Secondary series (kin.slider_crank.fourier_coefficients): d2s/dt2 = w^2 r [cos t + 4 A2 cos 2t
+ 16 A4 cos 4t + ...], 4 A2 = lam + lam^3/4 + ..., so to leading order the I4 secondary vertical
force is 4 m_rec r w^2 lam cos(2 theta).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kin.slider_crank import fourier_coefficients, piston_acceleration
from dyn import params

BALANCE_SHAFT_SPEED_RATIO = 2      # Lanchester shafts run at twice crank speed
BALANCE_SHAFT_COUNT = 2            # counter-rotating pair


@dataclass
class Shaking:
    theta: np.ndarray
    Fz: np.ndarray             # exact vertical force [N]
    Fy: np.ndarray             # lateral force from rotating masses [N]
    Fz_rot: np.ndarray         # vertical force from rotating masses [N]
    My: np.ndarray             # pitching moment [N m]
    Mz: np.ndarray             # yaw moment [N m]
    Fz_primary: np.ndarray
    Fz_secondary: np.ndarray   # series secondary (A2 to lam^5) [N]
    My_primary: np.ndarray
    My_secondary: np.ndarray


def rotating_unbalance_mr(spec=None, dspec=None) -> float:
    """Net rotating unbalance per throw after counterweights [kg m] (throw + rod big end)."""
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    g = params.geometry(spec)
    m = params.masses(spec)
    total = dspec["crank"]["throw_unbalance_mr"] + m.rod_big_end * g.crank_radius
    return (1.0 - dspec["crank"]["counterweight_fraction"]) * total


def shaking(theta, omega: float, spec=None, dspec=None) -> Shaking:
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    g = params.geometry(spec)
    m_rec = params.masses(spec).reciprocating
    r, l, lam = g.crank_radius, g.conrod_length, g.lam
    a2, _, _ = fourier_coefficients(lam)
    th = np.asarray(theta, dtype=float)
    thr = params.throw_angles(spec)
    x = params.cylinder_positions(spec, dspec)
    u = rotating_unbalance_mr(spec, dspec) * omega**2

    t_k = th[None, :] - thr[:, None]                      # local geometric angle of each throw
    fz = m_rec * piston_acceleration(t_k, omega, r, l)
    fz_p = m_rec * omega**2 * r * np.cos(t_k)
    fz_s = m_rec * omega**2 * r * 4 * a2 * np.cos(2 * t_k)
    # throw k sits at angle t_k from +z towards +y (rotation direction)
    fz_rot = u * np.cos(t_k)
    fy_rot = u * np.sin(t_k)
    xc = x[:, None]
    return Shaking(
        theta=th,
        Fz=fz.sum(0) + fz_rot.sum(0), Fy=fy_rot.sum(0), Fz_rot=fz_rot.sum(0),
        My=-(xc * (fz + fz_rot)).sum(0), Mz=(xc * fy_rot).sum(0),
        Fz_primary=fz_p.sum(0), Fz_secondary=fz_s.sum(0),
        My_primary=-(xc * fz_p).sum(0), My_secondary=-(xc * fz_s).sum(0),
    )


def secondary_force_analytic(theta, omega: float, spec=None):
    """Leading-order I4 secondary vertical force 4 m_rec r w^2 lam cos(2 theta) [N]."""
    spec = spec or params.engine_spec()
    g = params.geometry(spec)
    n = spec["engine"]["cylinders"]
    m_rec = params.masses(spec).reciprocating
    return n * m_rec * g.crank_radius * omega**2 * g.lam * np.cos(2 * np.asarray(theta, dtype=float))


def harmonic(signal, theta, order: float) -> complex:
    """Complex amplitude C of order `order` (signal ~ Re(C e^{i order theta})) over a uniform
    periodic grid covering whole periods of that order (endpoint excluded)."""
    s = np.asarray(signal, dtype=float)
    return complex(2.0 * np.mean(s * np.exp(-1j * order * np.asarray(theta))))


@dataclass
class LanchesterPair:
    me_per_shaft: float        # kg m eccentric unbalance of each shaft
    speed_ratio: int
    phase: float               # rad, eccentric angle (from +z) when theta = 0, both shafts


def lanchester_sizing(spec=None, use_exact_harmonic: bool = True, omega_ref: float = 1.0) -> LanchesterPair:
    """Size a counter-rotating pair at 2x speed whose vertical force cancels the 2nd-order force.

    Pair force: 2 m e (2w)^2 cos(2 theta + phase), lateral components cancel. With the eccentrics
    pointing DOWN (phase = pi) at theta = 0 it opposes F_2 cos(2 theta), giving
    m e = F_2 / (2 * (2w)^2). F_2 is the exact 2nd harmonic (default) or the 4 m r w^2 4A2 series.
    """
    spec = spec or params.engine_spec()
    g = params.geometry(spec)
    m_rec = params.masses(spec).reciprocating
    n = spec["engine"]["cylinders"]
    if use_exact_harmonic:
        th = np.linspace(0.0, params.REV, 4096, endpoint=False)
        f2 = abs(harmonic(shaking(th, omega_ref, spec).Fz, th, 2)) / omega_ref**2
    else:
        a2, _, _ = fourier_coefficients(g.lam)
        f2 = n * m_rec * g.crank_radius * 4 * a2
    k = BALANCE_SHAFT_SPEED_RATIO
    return LanchesterPair(me_per_shaft=f2 / (BALANCE_SHAFT_COUNT * k**2), speed_ratio=k, phase=np.pi)


def lanchester_force(theta, omega: float, pair: LanchesterPair):
    """(F_z, F_y) of the balance-shaft pair [N]; shafts at +/- k*theta + phase."""
    th = np.asarray(theta, dtype=float)
    k = pair.speed_ratio
    c = pair.me_per_shaft * (k * omega) ** 2
    a1 = k * th + pair.phase
    a2 = -k * th - pair.phase
    return c * (np.cos(a1) + np.cos(a2)), c * (np.sin(a1) + np.sin(a2))
