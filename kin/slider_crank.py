"""Exact slider-crank kinematics (zero pin offset) and the 4th-order Fourier series.

Crank angle theta in radians, 0 = TDC. Piston position x is measured from the crank
axis along the bore axis; s = distance below TDC.
"""
from __future__ import annotations

import numpy as np


def piston_position(theta, r, l):
    """Distance from crank axis to wrist pin."""
    theta = np.asarray(theta, dtype=float)
    return r * np.cos(theta) + np.sqrt(l**2 - (r * np.sin(theta)) ** 2)


def displacement_from_tdc(theta, r, l):
    return (r + l) - piston_position(theta, r, l)


def piston_velocity(theta, omega, r, l):
    """ds/dt (positive moving away from TDC), exact."""
    theta = np.asarray(theta, dtype=float)
    s, c = np.sin(theta), np.cos(theta)
    root = np.sqrt(l**2 - (r * s) ** 2)
    return omega * r * s * (1 + r * c / root)


def piston_acceleration(theta, omega, r, l):
    """d2s/dt2 at constant omega, exact."""
    theta = np.asarray(theta, dtype=float)
    s, c = np.sin(theta), np.cos(theta)
    root2 = l**2 - (r * s) ** 2
    root = np.sqrt(root2)
    term = (r * (c**2 - s**2)) / root + (r**3 * s**2 * c**2) / root2**1.5
    return omega**2 * r * (c + term)


def fourier_coefficients(lam: float):
    """Coefficients A2, A4, A6 of x/r ~ cos(t) + sum A_n cos(n t) (+ const), series in lam.

    Standard expansion: x = r[ cos t + (1/lam)(1 - lam^2/4 - 3lam^4/64) + (lam/4 + lam^3/16) cos 2t
    - (lam^3/64) cos 4t + ... ]; acceleration uses the same harmonics.
    """
    a2 = lam / 4 + lam**3 / 16 + 15 * lam**5 / 512
    a4 = -(lam**3 / 64 + 3 * lam**5 / 256)
    a6 = lam**5 / 512
    return a2, a4, a6


def acceleration_series(theta, omega, r, l, order: int = 4):
    """Harmonic series for d2s/dt2 truncated at `order` (2, 4 or 6)."""
    lam = r / l
    a2, a4, a6 = fourier_coefficients(lam)
    theta = np.asarray(theta, dtype=float)
    acc = np.cos(theta) + 4 * a2 * np.cos(2 * theta)
    if order >= 4:
        acc = acc + 16 * a4 * np.cos(4 * theta)
    if order >= 6:
        acc = acc + 36 * a6 * np.cos(6 * theta)
    return omega**2 * r * acc


def cylinder_volume(theta, geom):
    """Instantaneous cylinder volume [m^3]."""
    return geom.clearance_volume + geom.piston_area * displacement_from_tdc(
        theta, geom.crank_radius, geom.conrod_length
    )


def dV_dtheta(theta, geom):
    """dV/dtheta [m^3/rad]."""
    return geom.piston_area * piston_velocity(theta, 1.0, geom.crank_radius, geom.conrod_length)


def conrod_angle(theta, r, l):
    return np.arcsin(r * np.sin(np.asarray(theta, dtype=float)) / l)
