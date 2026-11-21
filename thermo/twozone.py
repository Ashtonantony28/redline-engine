"""Closed-cycle (IVC -> EVO) two-zone model with a prescribed mass-fraction-burned law.

Zones share the pressure p. Unknowns p, T_u, T_b; with m_b = x_b m and ' = d/dtheta [rad]:

  unburned (frozen, loses mass at its own enthalpy):     m_u cp_u T_u' = V_u p' - Q_u'
  burned (equilibrium, gains mass at h_u):  m_b (h_T T_b' + (h_p - v_b) p') = -Q_b' + (h_u - h_b) m_b'
  volume:   V = m_u R_u T_u / p + m_b R_b(T_b, p) T_b / p

Eliminating T_u' and T_b' from the differentiated volume constraint gives p' = B / A (see `_rhs`).
In the single-zone limit (m_b = 0) this reduces to p' = -gamma p V'/V - (gamma - 1) Q'/V.
Heat loss is split between zones by volume fraction (keeps Q_b'/m_b finite as m_b -> 0).
Integration: classical fixed-step RK4 on a uniform crank-angle grid; V, V', x_b, x_b' are
sampled by the caller on the half-step grid (2N+1 points).
"""
from __future__ import annotations

import numpy as np
from numba import njit

from thermo.props import burned, frozen, T_adiabatic_p, burn_constant_volume
from thermo.woschni import coefficient, wall_terms, HP_ON, HP_OMEGA

# state indices
Y_P, Y_TU, Y_TB, Y_W, Y_QU, Y_QB = range(6)
NY = 6
# output columns
O_P, O_TU, O_TB, O_XB, O_W, O_QU, O_QB, O_HC = range(8)
NO = 8


@njit(cache=True)
def _rhs(y, V, dV, xb, dxb, active, m, xr, Fz, T0f, dTf, R_F, R_P, B, T0b, dTb, lp0, dlp, hp, dy):
    p, Tu, Tb = y[Y_P], y[Y_TU], y[Y_TB]
    if not active:   # burned zone not yet created: stages that straddle its creation see m_b = 0
        xb = 0.0
        dxb = 0.0
    mb = xb * m
    mu = m - mb
    dmb = dxb * m
    hu, cpu, Ru = frozen(Tu, xr, Fz, T0f, dTf, R_F, R_P)
    vu = Ru * Tu / p
    Vu = mu * vu
    if active:
        hb, Rb, hT, hpb, RT, Rp = burned(Tb, p, B, T0b, dTb, lp0, dlp)
        vb = Rb * Tb / p
        vT = (Rb + Tb * RT) / p
        vp = vb * (Rp / Rb - 1.0 / p)
        rb = dxb / xb
    else:
        hb = Rb = hT = hpb = RT = Rp = vb = vT = vp = rb = 0.0
        hT = 1.0
        Rb = 1.0
    Vb = mb * vb
    # heat loss per unit zone mass [J/kg/rad]
    qu = 0.0
    qb = 0.0
    hc = 0.0
    if hp[HP_ON] > 0.0:
        Tm = p * V / (mu * Ru + mb * Rb * (1.0 if active else 0.0))
        hc = coefficient(p, Tm, V, active, hp)
        Aw, ATw = wall_terms(V, hp)
        k = hc / (hp[HP_OMEGA] * V)
        qu = k * vu * (Aw * Tu - ATw)
        if active:
            qb = k * vb * (Aw * Tb - ATw)
    g = Ru / (cpu * p)
    A = Vu * g - Vu / p
    Bn = dV + g * mu * qu
    if active:
        A += mb * vp - vT * (hpb - vb) * mb / hT
        Bn += -(vb - vu) * dmb - vT * (-mb * qb + (hu - hb) * dmb) / hT
    dp = Bn / A
    dy[Y_P] = dp
    dy[Y_TU] = (vu * dp - qu) / cpu
    if active:
        dy[Y_TB] = (-qb + (hu - hb) * rb - (hpb - vb) * dp) / hT
    else:
        dy[Y_TB] = 0.0
    dy[Y_W] = p * dV
    dy[Y_QU] = mu * qu
    dy[Y_QB] = mb * qb
    return hc


@njit(cache=True)
def integrate(y0, N, h, Vh, dVh, xbh, dxbh, m, xr, k_act, k_inst,
              Fz, T0f, dTf, R_F, R_P, B, T0b, dTb, lp0, dlp, hp, tol, maxit, T_flame_guess):
    """RK4 over N steps of size h [rad]. k_act: grid index where the Wiebe burned zone is created
    (T_b = constant-pressure adiabatic flame temperature); k_inst: grid index of an instantaneous
    constant-volume burn. Pass -1 to disable either. Returns (N+1, NO) outputs."""
    out = np.zeros((N + 1, NO))
    y = y0.copy()
    k1 = np.empty(NY)
    k2 = np.empty(NY)
    k3 = np.empty(NY)
    k4 = np.empty(NY)
    yt = np.empty(NY)
    active = False
    for k in range(N + 1):
        i = 2 * k
        if k == k_act:
            hu, cpu, Ru = frozen(y[Y_TU], xr, Fz, T0f, dTf, R_F, R_P)
            y[Y_TB] = T_adiabatic_p(hu, y[Y_P], T_flame_guess, B, T0b, dTb, lp0, dlp, tol, maxit)
            active = True
        if k == k_inst:
            hu, cpu, Ru = frozen(y[Y_TU], xr, Fz, T0f, dTf, R_F, R_P)
            Tb, pb = burn_constant_volume(hu - Ru * y[Y_TU], m, Vh[i], T_flame_guess,
                                          B, T0b, dTb, lp0, dlp, tol, maxit)
            y[Y_TB] = Tb
            y[Y_P] = pb
            active = True
        hc = _rhs(y, Vh[i], dVh[i], xbh[i], dxbh[i], active, m, xr, Fz, T0f, dTf, R_F, R_P,
                  B, T0b, dTb, lp0, dlp, hp, k1)
        out[k, O_P] = y[Y_P]
        out[k, O_TU] = y[Y_TU]
        out[k, O_TB] = y[Y_TB] if active else np.nan
        out[k, O_XB] = xbh[i]
        out[k, O_W] = y[Y_W]
        out[k, O_QU] = y[Y_QU]
        out[k, O_QB] = y[Y_QB]
        out[k, O_HC] = hc
        if k == N:
            break
        for j in range(NY):
            yt[j] = y[j] + 0.5 * h * k1[j]
        _rhs(yt, Vh[i + 1], dVh[i + 1], xbh[i + 1], dxbh[i + 1], active, m, xr, Fz, T0f, dTf, R_F, R_P,
             B, T0b, dTb, lp0, dlp, hp, k2)
        for j in range(NY):
            yt[j] = y[j] + 0.5 * h * k2[j]
        _rhs(yt, Vh[i + 1], dVh[i + 1], xbh[i + 1], dxbh[i + 1], active, m, xr, Fz, T0f, dTf, R_F, R_P,
             B, T0b, dTb, lp0, dlp, hp, k3)
        for j in range(NY):
            yt[j] = y[j] + h * k3[j]
        _rhs(yt, Vh[i + 2], dVh[i + 2], xbh[i + 2], dxbh[i + 2], active, m, xr, Fz, T0f, dTf, R_F, R_P,
             B, T0b, dTb, lp0, dlp, hp, k4)
        for j in range(NY):
            y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
    return out
