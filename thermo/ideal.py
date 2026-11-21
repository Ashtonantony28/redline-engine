"""Air-standard cycles used as analytic anchors (M0)."""
from __future__ import annotations

import numpy as np

from kin.slider_crank import cylinder_volume


def otto_efficiency(r: float, gamma: float = 1.4) -> float:
    return 1.0 - r ** (1.0 - gamma)


def otto_cycle_numeric(geom, gamma=1.4, p1=1.0e5, T1=300.0, q_in_per_kg=1.8e6, n=20001, R=287.0):
    """Integrate an air-standard Otto cycle over the real slider-crank volume trace.

    Isentropic compression BDC->TDC, constant-volume heat addition at TDC, isentropic
    expansion back to BDC. Returns (efficiency, imep, (V, p) traces).
    """
    th_c = np.linspace(-np.pi, 0.0, n)
    th_e = np.linspace(0.0, np.pi, n)
    Vc = cylinder_volume(th_c, geom)
    Ve = cylinder_volume(th_e, geom)
    V1 = Vc[0]
    m = p1 * V1 / (R * T1)
    pc = p1 * (V1 / Vc) ** gamma
    cv = R / (gamma - 1)
    T2 = T1 * (V1 / Vc[-1]) ** (gamma - 1)
    T3 = T2 + q_in_per_kg / cv
    p3 = m * R * T3 / Ve[0]
    pe = p3 * (Ve[0] / Ve) ** gamma
    w = np.trapezoid(pe, Ve) - np.trapezoid(pc[::-1], Vc[::-1])
    q = m * q_in_per_kg
    swept = Ve[-1] - Ve[0]
    return w / q, w / swept, (np.concatenate([Vc, Ve]), np.concatenate([pc, pe]))
