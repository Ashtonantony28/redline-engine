"""Fatigue, buckling and stress-concentration helpers (closed forms, cited).

- Goodman (modified Goodman line, Shigley "Mechanical Engineering Design", 10th ed., eq. 6-41):
  1/n = sigma_a/S_e + sigma_m/S_ut for sigma_m >= 0; for compressive mean stress the mean term is
  dropped, n = S_e / sigma_a (Shigley sec. 6-12, conservative for ductile metals).
- Multiaxial: equivalent alternating / mean stresses are von Mises values of the alternating and mean
  stress tensors (Shigley sec. 6-14), the mean signed by the first stress invariant; for load histories
  the signed von Mises range at each point gives the alternating and mean values.
- Column buckling: Euler P_cr = pi^2 E I / (K L)^2 for slenderness above C_c = sqrt(2 pi^2 E / S_y);
  J.B. Johnson parabola P_cr = A S_y (1 - S_y (K L / rho)^2 / (4 pi^2 E)) below it (Shigley eq. 4-43, 4-49).
- Peterson-type Kt of a stepped round bar in bending (Pilkey & Pilkey, "Peterson's Stress Concentration
  Factors", 3rd ed., chart 3.10 curve fit), used only as a sanity check of the FEA fillet Kt.
"""
from __future__ import annotations

import math

import numpy as np

from fea.ccx import von_mises


def signed_von_mises(s: np.ndarray) -> np.ndarray:
    """Von Mises magnitude with the sign of the first stress invariant (hydrostatic part), (..., 6) -> (...).
    Uniaxial and near-uniaxial states (fillets, free surfaces) keep the sign of the dominant stress."""
    tr = s[..., 0] + s[..., 1] + s[..., 2]
    return np.where(tr >= 0.0, 1.0, -1.0) * von_mises(s)


def goodman_sf(sigma_a, sigma_m, se, sut):
    """Modified Goodman safety factor (elementwise); se, sut may be arrays (temperature derating)."""
    sigma_a = np.asarray(sigma_a, dtype=float)
    sigma_m = np.asarray(sigma_m, dtype=float)
    tiny = 1.0e-30
    inv = sigma_a / se + np.clip(sigma_m, 0.0, None) / sut
    return 1.0 / np.maximum(inv, tiny)


def goodman_from_history(svm_hist):
    """From a signed-von-Mises history (..., n_time): (sigma_a, sigma_m) per point."""
    hi, lo = svm_hist.max(axis=-1), svm_hist.min(axis=-1)
    return 0.5 * (hi - lo), 0.5 * (hi + lo)


def johnson_euler(A, I, L_eff, E, sy):
    """Critical load of a column; returns (P_cr, regime, slenderness, C_c)."""
    rho = math.sqrt(I / A)
    sr = L_eff / rho
    cc = math.sqrt(2 * math.pi**2 * E / sy)
    if sr >= cc:
        return math.pi**2 * E * I / L_eff**2, "Euler", sr, cc
    return A * sy * (1 - sy * sr**2 / (4 * math.pi**2 * E)), "Johnson", sr, cc


def peterson_kt_bending(D, d, r):
    """Kt of a shoulder fillet in a round bar under bending (chart 3.10 fit), h = (D - d) / 2."""
    h = (D - d) / 2
    q = h / r
    s = math.sqrt(q)
    if q <= 2.0:
        c1 = 0.947 + 1.206 * s - 0.131 * q
        c2 = 0.022 - 3.405 * s + 0.915 * q
        c3 = 0.869 + 1.777 * s - 0.555 * q
        c4 = -0.810 + 0.422 * s - 0.260 * q
    else:
        c1 = 1.232 + 0.832 * s - 0.008 * q
        c2 = -3.813 + 0.968 * s - 0.260 * q
        c3 = 7.423 - 4.868 * s + 0.869 * q
        c4 = -3.839 + 3.070 * s - 0.600 * q
    t = 2 * h / D
    return c1 + c2 * t + c3 * t**2 + c4 * t**3


def richardson(f, ratio):
    """Richardson extrapolation of three solutions f = (coarse, medium, fine) on a constant refinement
    ratio. Returns (observed order p or None, extrapolated value or None, GCI_fine or None)."""
    f0, f1, f2 = (float(v) for v in f)
    d1, d2 = f1 - f0, f2 - f1
    if d1 == 0.0 or d2 == 0.0 or d1 * d2 < 0.0:
        return None, None, None          # oscillatory or exact: no monotone extrapolation
    p = math.log(abs(d1 / d2)) / math.log(ratio)
    fe = f2 + d2 / (ratio**p - 1.0)
    fs = 1.25                            # Roache safety factor for three-grid studies
    gci = fs * abs(d2 / f2) / (ratio**p - 1.0)
    return p, fe, gci
