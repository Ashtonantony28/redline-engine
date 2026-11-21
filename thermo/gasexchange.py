"""Filling-and-emptying gas exchange (EVO -> IVC), single zone, quasi-steady valve flow.

The cylinder is one well-mixed open system of fresh charge (F) and frozen products (P) with burned
mass fraction x = m_P / m. Conservative states (m, m x, U) are integrated with fixed-step RK4, so
mass and energy bookkeeping are exact to round-off:

    m'  = sum(in) - out
    (m x)' = sum(x_in in) - x out
    U'  = -p V' - Q' + sum(h_in in) - h out

Manifolds are constant-pressure, constant-temperature reservoirs (spec/thermo.yaml `loads`). Gas
pushed back into the intake port is held in a port store (mass, products, enthalpy) and is
re-inducted first (blended smoothly into fresh charge over the last `store_blend_mass_fraction` of it); exhaust backflow enters as products at the exhaust manifold temperature.
Valve flow: compressible orifice (Heywood App. C) with Cd on curtain area; above the pressure ratio
`linear_flow_pr` the orifice law is linearised to remove the infinite d(mdot)/d(dp) at dp = 0.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from thermo.props import frozen, T_from_u_frozen
from thermo.woschni import coefficient, wall_terms, HP_ON, HP_OMEGA

(G_M, G_MX, G_U, G_MF_IN, G_MS_IN, G_MI_OUT, G_ME_OUT, G_ME_IN, G_W, G_Q, G_MS, G_MXS, G_HS) = range(13)
NG = 13
(BC_PIM, BC_TIM, BC_PEM, BC_TEM, BC_PRLIN, BC_TOL, BC_MAXIT, BC_MEPS) = range(8)
NBC = 8
(GO_P, GO_T, GO_X, GO_M, GO_MF_IN, GO_MS_IN, GO_MI_OUT, GO_ME_OUT, GO_ME_IN, GO_W, GO_Q,
 GO_MS, GO_DMI, GO_DME, GO_HC) = range(15)
NGO = 15


@njit(cache=True)
def _psi(pr, g):
    return np.sqrt(2.0 * g / (g - 1.0) * (pr ** (2.0 / g) - pr ** ((g + 1.0) / g)))


@njit(cache=True)
def orifice(cdA, p0, T0, R, g, pd, pr_lin):
    """Mass flow [kg/s] from upstream stagnation (p0, T0) to downstream static pd."""
    if cdA <= 0.0 or pd >= p0:
        return 0.0
    pr = pd / p0
    crit = (2.0 / (g + 1.0)) ** (g / (g - 1.0))
    if pr < crit:
        pr = crit
    if pr > pr_lin:
        psi = _psi(pr_lin, g) * (1.0 - pr) / (1.0 - pr_lin)
    else:
        psi = _psi(pr, g)
    return cdA * p0 / np.sqrt(R * T0) * psi


@njit(cache=True)
def _rhs(y, V, dV, Ai, Ae, Tg, bc, Fz, T0f, dTf, R_F, R_P, hp, dy, aux):
    m = y[G_M]
    x = y[G_MX] / m
    T = T_from_u_frozen(y[G_U] / m, x, Tg[0], Fz, T0f, dTf, R_F, R_P, bc[BC_TOL], int(bc[BC_MAXIT]))
    Tg[0] = T
    h, cp, R = frozen(T, x, Fz, T0f, dTf, R_F, R_P)
    p = m * R * T / V
    g = cp / (cp - R)
    om = hp[HP_OMEGA]
    pim, Tim, pem, Tem, prl = bc[BC_PIM], bc[BC_TIM], bc[BC_PEM], bc[BC_TEM], bc[BC_PRLIN]

    mf_in = 0.0   # fresh from manifold
    ms_in = 0.0   # from port store
    mi_out = 0.0  # cylinder -> intake port store
    me_out = 0.0
    me_in = 0.0
    h_s = 0.0
    x_s = 0.0
    if Ai > 0.0:
        if pim >= p:
            hF, cpF, RF = frozen(Tim, 0.0, Fz, T0f, dTf, R_F, R_P)
            mdot = orifice(Ai, pim, Tim, RF, cpF / (cpF - RF), p, prl) / om
            # draw from the port store first; blend smoothly to fresh charge as the store empties
            # (store fraction phi = clip(m_s / m_eps, 0, 1) keeps the RHS continuous for RK4)
            phi = min(max(y[G_MS] / bc[BC_MEPS], 0.0), 1.0)
            if phi > 0.0:
                h_s = y[G_HS] / y[G_MS]
                x_s = min(max(y[G_MXS] / y[G_MS], 0.0), 1.0)
            ms_in = phi * mdot
            mf_in = (1.0 - phi) * mdot
        else:
            mi_out = orifice(Ai, p, T, R, g, pim, prl) / om
    hFin, cpFin, RFin = frozen(Tim, 0.0, Fz, T0f, dTf, R_F, R_P)
    hPex = 0.0
    if Ae > 0.0:
        if p >= pem:
            me_out = orifice(Ae, p, T, R, g, pem, prl) / om
        else:
            hPex, cpP, RP = frozen(Tem, 1.0, Fz, T0f, dTf, R_F, R_P)
            me_in = orifice(Ae, pem, Tem, RP, cpP / (cpP - RP), p, prl) / om
    qdot = 0.0
    hc = 0.0
    if hp[HP_ON] > 0.0:
        hc = coefficient(p, T, V, False, hp)
        Aw, ATw = wall_terms(V, hp)
        qdot = hc * (Aw * T - ATw) / om
    out = mi_out + me_out
    dy[G_M] = mf_in + ms_in + me_in - out
    dy[G_MX] = x_s * ms_in + me_in - x * out
    dy[G_U] = -p * dV - qdot + hFin * mf_in + h_s * ms_in + hPex * me_in - h * out
    dy[G_MF_IN] = mf_in
    dy[G_MS_IN] = ms_in
    dy[G_MI_OUT] = mi_out
    dy[G_ME_OUT] = me_out
    dy[G_ME_IN] = me_in
    dy[G_W] = p * dV
    dy[G_Q] = qdot
    dy[G_MS] = mi_out - ms_in
    dy[G_MXS] = x * mi_out - x_s * ms_in
    dy[G_HS] = h * mi_out - h_s * ms_in
    aux[0] = p
    aux[1] = T
    aux[2] = x
    aux[3] = mf_in + ms_in - mi_out
    aux[4] = me_out - me_in
    aux[5] = hc


@njit(cache=True)
def integrate(y0, T0, N, h, Vh, dVh, Aih, Aeh, bc, Fz, T0f, dTf, R_F, R_P, hp):
    out = np.zeros((N + 1, NGO))
    y = y0.copy()
    k1 = np.empty(NG)
    k2 = np.empty(NG)
    k3 = np.empty(NG)
    k4 = np.empty(NG)
    yt = np.empty(NG)
    aux = np.empty(6)
    Tg = np.array([T0])
    for k in range(N + 1):
        i = 2 * k
        _rhs(y, Vh[i], dVh[i], Aih[i], Aeh[i], Tg, bc, Fz, T0f, dTf, R_F, R_P, hp, k1, aux)
        out[k, GO_P] = aux[0]
        out[k, GO_T] = aux[1]
        out[k, GO_X] = aux[2]
        out[k, GO_M] = y[G_M]
        out[k, GO_MF_IN] = y[G_MF_IN]
        out[k, GO_MS_IN] = y[G_MS_IN]
        out[k, GO_MI_OUT] = y[G_MI_OUT]
        out[k, GO_ME_OUT] = y[G_ME_OUT]
        out[k, GO_ME_IN] = y[G_ME_IN]
        out[k, GO_W] = y[G_W]
        out[k, GO_Q] = y[G_Q]
        out[k, GO_MS] = y[G_MS]
        out[k, GO_DMI] = aux[3]
        out[k, GO_DME] = aux[4]
        out[k, GO_HC] = aux[5]
        if k == N:
            break
        for j in range(NG):
            yt[j] = y[j] + 0.5 * h * k1[j]
        _rhs(yt, Vh[i + 1], dVh[i + 1], Aih[i + 1], Aeh[i + 1], Tg, bc, Fz, T0f, dTf, R_F, R_P, hp, k2, aux)
        for j in range(NG):
            yt[j] = y[j] + 0.5 * h * k2[j]
        _rhs(yt, Vh[i + 1], dVh[i + 1], Aih[i + 1], Aeh[i + 1], Tg, bc, Fz, T0f, dTf, R_F, R_P, hp, k3, aux)
        for j in range(NG):
            yt[j] = y[j] + h * k3[j]
        _rhs(yt, Vh[i + 2], dVh[i + 2], Aih[i + 2], Aeh[i + 2], Tg, bc, Fz, T0f, dTf, R_F, R_P, hp, k4, aux)
        for j in range(NG):
            y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
    return out, y
