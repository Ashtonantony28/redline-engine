"""End-gas knock (milestone M7): Livengood-Wu integral of the Douaud-Eyzat autoignition delay.

    tau(p, T) = A (ON / ON_ref)^a p[atm]^b exp(T_a / T)        (Douaud & Eyzat, SAE 780080, 1978)
    I(t)      = integral_{IVC}^{t} dt' / tau(p(t'), T_u(t'))   (Livengood & Wu, 1957)

Autoignition (knock onset) is the angle where I first reaches the threshold (1), counted as knock only if
the mass fraction burned there is below `end_gas_xb_limit`. The knock-limited spark advance (KLSA) is the
first advance (from the retarded side) at which the knock margin I(theta at x_b = limit) reaches the
threshold, found by a scan plus bisection. Constants: spec/combustion_ext.yaml `knock`. MODELLED results only.
"""
from __future__ import annotations

import math

import numpy as np

from thermo import twozone as tz

SEC_PER_MIN = 60.0


def tau_douaud_eyzat(p, T, kspec: dict):
    """Autoignition delay [s] at pressure p [Pa], unburned temperature T [K]."""
    de = kspec["douaud_eyzat"]
    on = float(kspec["octane_number"])
    return (float(de["A"]) * (on / float(de["on_ref"])) ** float(de["on_exp"])
            * (np.asarray(p, dtype=float) / float(de["p_unit"])) ** float(de["p_exp"])
            * np.exp(float(de["T_activation"]) / np.asarray(T, dtype=float)))


def livengood_wu(t, tau):
    """Cumulative Livengood-Wu integral of 1/tau over the time samples t (trapezoidal rule, exact for a
    constant tau: I(t) = (t - t0) / tau)."""
    t = np.asarray(t, dtype=float)
    r = 1.0 / np.asarray(tau, dtype=float)
    I = np.zeros_like(t)
    I[1:] = np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(t))
    return I


def knock_analysis(cl, rpm: float, kspec: dict) -> dict:
    """Knock integral of a closed two-zone run (thermo.cycle.ClosedResult) at engine speed rpm.

    Returns the integral at x_b = end_gas_xb_limit (`margin`, knock if >= threshold), the onset angle
    [deg, ADR 0001] where the integral reaches the threshold (NaN if never), x_b at onset, and `knock`."""
    omega = rpm * 2.0 * math.pi / SEC_PER_MIN
    th = cl.theta
    t = np.radians(th - th[0]) / omega
    tau = tau_douaud_eyzat(cl.col(tz.O_P), cl.col(tz.O_TU), kspec)
    I = livengood_wu(t, tau)
    xb = cl.col(tz.O_XB)
    xlim = float(kspec["end_gas_xb_limit"])
    thr = float(kspec["integral_threshold"])
    if xb[-1] >= xlim:
        j = int(np.argmax(xb >= xlim))
        f = 0.0 if j == 0 else (xlim - xb[j - 1]) / (xb[j] - xb[j - 1])
        margin = I[j - 1] + f * (I[j] - I[j - 1]) if j > 0 else I[0]
    else:
        margin = I[-1]
    onset = float("nan")
    xb_on = float("nan")
    if I[-1] >= thr:
        j = int(np.argmax(I >= thr))
        f = (thr - I[j - 1]) / (I[j] - I[j - 1]) if j > 0 else 0.0
        onset = float(th[j - 1] + f * (th[j] - th[j - 1])) if j > 0 else float(th[0])
        xb_on = float(np.interp(onset, th, xb))
    return {"margin": float(margin), "onset_deg": onset, "xb_onset": xb_on,
            "knock": bool(margin >= thr), "I_end": float(I[-1])}


def klsa_search(margin_fn, adv_min: float, adv_max: float, tol: float, threshold: float, scan_step: float):
    """Knock-limited spark advance: the first advance, scanning from adv_min towards adv_max in scan_step
    increments, at which margin_fn(adv) reaches threshold, refined by bisection to tol.

    The margin is not monotonic over the whole window: with extreme advance most of the charge burns
    before the end gas is compressed, so the margin at x_b = limit falls again. The physical limit is
    the first crossing from the retarded side. Returns (klsa_deg, status) with status "limited",
    "unlimited" (no knock anywhere in the window, klsa = adv_max) or "unavoidable" (knock even at
    adv_min; klsa = adv_min)."""
    if margin_fn(adv_min) >= threshold:
        return adv_min, "unavoidable"
    n = max(1, int(math.ceil((adv_max - adv_min) / scan_step - 1e-9)))
    lo = adv_min
    for i in range(1, n + 1):
        hi = min(adv_min + i * scan_step, adv_max)
        if margin_fn(hi) >= threshold:
            break
        lo = hi
    else:
        return adv_max, "unlimited"
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if margin_fn(mid) < threshold:
            lo = mid
        else:
            hi = mid
    return lo, "limited"
