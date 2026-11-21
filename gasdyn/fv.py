"""1D finite-volume compressible-flow kernels for a network of pipes and 0D volumes (milestone M4).

Scheme
------
* Euler equations for an ideal gas with constant gamma per pipe, plus one passive scalar Y (burned-gas
  mass fraction). Conserved cell state: (rho, rho u, E, rho Y), E = p / (gamma - 1) + rho u^2 / 2.
* HLLC approximate Riemann solver (Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics",
  3rd ed., Springer 2009, sec. 10.4-10.6) with pressure-based wave-speed estimates (Toro sec. 10.5.2).
* MUSCL reconstruction of the primitives (rho, u, p, Y) with a slope limiter (van Leer or minmod),
  second-order SSP Runge-Kutta (Heun) in time (Shu & Osher 1988), CFL-controlled global time step.
* Positivity: a cell whose reconstructed face density or pressure would fall below a fraction of the cell
  value drops to first order; after each stage density and pressure are floored (floor hits are counted
  and returned, and are expected to be zero in normal runs).
* Wall friction (Darcy factor: laminar 64/Re or Blasius 0.3164 Re^-0.25, whichever is larger) and wall
  heat transfer by the Reynolds / Chilton-Colburn analogy St = (f / 8) Pr^(-2/3) as source terms.

Pipe ends
---------
Every pipe end is either a wall or connects to a 0D volume through an orifice of effective area A_eff
(A_eff = A_pipe is a plain open end). The boundary state is found by intersecting

* the Riemann invariant arriving from the pipe interior, v_b = v_i + 2 c_i / (g - 1) ((p_b / p_i)^beta - 1),
  beta = (g - 1) / (2 g), v = velocity directed into the pipe, with
* the quasi-steady isentropic orifice law (Heywood App. C; thermo.gasexchange.orifice) between the volume
  and the boundary pressure p_b (the "constant-pressure" valve boundary of Benson, "The Thermodynamics and
  Gas Dynamics of Internal Combustion Engines", vol. 1, OUP 1982, ch. 7, here with one scalar root solve).

The scalar equation in p_b is monotone and bracketed, so the boundary is unconditionally well posed even
for tiny valve areas and choked flow. The same mass / energy rates are applied to the pipe face and to
the volume, so mass and energy of the network are conserved to round-off.

0D volumes: ambient (fixed state), fixed-volume ideal-gas plenum, or engine cylinder during gas exchange
(variable volume from kin.slider_crank tables, frozen-mixture properties from thermo.props, Woschni wall
heat transfer from thermo.woschni). All numbers are MODELLED (there is no physical engine).
"""
from __future__ import annotations

import numpy as np
from numba import njit

from thermo.gasexchange import orifice
from thermo.props import frozen, T_from_u_frozen
from thermo.woschni import coefficient, wall_terms, HP_ON

NVAR = 4                                   # rho, rho u, E, rho Y
# volume state columns: mass, burned mass, energy, then diagnostics integrated with the same RK scheme:
# p dV work, wall heat, net fresh charge through the intake valves, net mass out through the exhaust
# valves, net mass received (cylinders: through the intake valves; plenum / ambient: through all ends)
(S_M, S_MY, S_E, S_W, S_Q, S_FRESH, S_EXH, S_NET) = range(8)
NS = 8
VK_AMBIENT, VK_PLENUM, VK_CYL = 0, 1, 2    # volume kinds
SLOT_CONST, SLOT_INTAKE, SLOT_EXHAUST = -1, 0, 1
LIM_VANLEER, LIM_MINMOD = 0, 1
# numerics vector
(NU_CFL, NU_RHO_FLOOR, NU_P_FLOOR, NU_LIMITER, NU_POS_FRAC, NU_MU_REF, NU_T_SUTH_REF, NU_S_SUTH, NU_PR,
 NU_F_LAM, NU_F_BLASIUS, NU_F_BLASIUS_EXP, NU_RE_MIN, NU_BC_TOL, NU_BC_MAXIT, NU_DT_MAX, NU_NEWTON_TOL,
 NU_NEWTON_MAXIT, NU_TAB_DTHETA, NU_CHAR_BASE_MIN) = range(20)
N_NU = 20
NO_PRESSURE_LINEARISATION = 1.0   # thermo.gasexchange.orifice: pure orifice law (the boundary solve is implicit)
FOUR_PI = 4.0 * np.pi             # one four-stroke cycle in crank radians (ADR 0001)


# ---------------------------------------------------------------------------------------------
# HLLC flux
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def hllc(rl, ul, pl, yl, rr, ur, pr, yr, g, f):
    """HLLC flux of (rho, rho u, E, rho Y) for left / right primitive states; result written to f[0:4]."""
    cl = np.sqrt(g * pl / rl)
    cr = np.sqrt(g * pr / rr)
    # pressure-based wave speeds (Toro eq. 10.59-10.61, PVRS estimate of p*)
    rbar = 0.5 * (rl + rr)
    cbar = 0.5 * (cl + cr)
    ps = 0.5 * (pl + pr) - 0.5 * (ur - ul) * rbar * cbar
    if ps < 0.0:
        ps = 0.0
    gq = (g + 1.0) / (2.0 * g)
    ql = 1.0
    if ps > pl:
        ql = np.sqrt(1.0 + gq * (ps / pl - 1.0))
    qr = 1.0
    if ps > pr:
        qr = np.sqrt(1.0 + gq * (ps / pr - 1.0))
    sl = ul - cl * ql
    sr = ur + cr * qr
    el = pl / (g - 1.0) + 0.5 * rl * ul * ul
    er = pr / (g - 1.0) + 0.5 * rr * ur * ur
    if sl >= 0.0:
        f[0] = rl * ul
        f[1] = rl * ul * ul + pl
        f[2] = ul * (el + pl)
    elif sr <= 0.0:
        f[0] = rr * ur
        f[1] = rr * ur * ur + pr
        f[2] = ur * (er + pr)
    else:
        sm = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur))
        if sm >= 0.0:
            k = rl * (sl - ul) / (sl - sm)
            es = k * (el / rl + (sm - ul) * (sm + pl / (rl * (sl - ul))))
            f[0] = rl * ul + sl * (k - rl)
            f[1] = rl * ul * ul + pl + sl * (k * sm - rl * ul)
            f[2] = ul * (el + pl) + sl * (es - el)
        else:
            k = rr * (sr - ur) / (sr - sm)
            es = k * (er / rr + (sm - ur) * (sm + pr / (rr * (sr - ur))))
            f[0] = rr * ur + sr * (k - rr)
            f[1] = rr * ur * ur + pr + sr * (k * sm - rr * ur)
            f[2] = ur * (er + pr) + sr * (es - er)
    f[3] = f[0] * (yl if f[0] >= 0.0 else yr)


@njit(cache=True)
def _slope(a, b, lim):
    """Limited slope from backward (a) and forward (b) differences."""
    if a * b <= 0.0:
        return 0.0
    if lim == LIM_MINMOD:
        return a if abs(a) < abs(b) else b
    return 2.0 * a * b / (a + b)      # van Leer (harmonic mean)


# ---------------------------------------------------------------------------------------------
# pipe-end boundary: Riemann invariant x orifice law
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def _bc_residual(inflow, pb, A, Ap, rho_i, v_i, p_i, c_i, g, R, pV, TV, gV, RV, res):
    """Residual orifice mass flow - pipe-side mass flow at boundary pressure pb.
    res: [rho_b, |v_b|, T0 of the stream]."""
    beta = (g - 1.0) / (2.0 * g)
    cp = g * R / (g - 1.0)
    v = v_i + 2.0 * c_i / (g - 1.0) * ((pb / p_i) ** beta - 1.0)
    if inflow:
        Tb = TV - 0.5 * v * v / cp
        Tmin = TV * 2.0 / (g + 1.0)          # sonic static temperature: the pipe-side stream stays subsonic
        if Tb < Tmin:
            Tb = Tmin
        rho_b = pb / (R * Tb)
        res[0] = rho_b
        res[1] = v
        res[2] = TV
        return orifice(A, pV, TV, RV, gV, pb, NO_PRESSURE_LINEARISATION) - rho_b * Ap * v
    w = -v
    rho_b = rho_i * (pb / p_i) ** (1.0 / g)
    Tb = pb / (rho_b * R)
    T0 = Tb + 0.5 * w * w / cp
    p0 = pb * (T0 / Tb) ** (g / (g - 1.0))
    res[0] = rho_b
    res[1] = w
    res[2] = T0
    return orifice(A, p0, T0, R, g, pV, NO_PRESSURE_LINEARISATION) - rho_b * Ap * w


@njit(cache=True)
def end_bc(Aeff, Ap, rho_i, v_i, p_i, g, R, pV, TV, gV, RV, tol, maxit, base_min, out):
    """Boundary state of a pipe end attached to a volume at (pV, TV) through effective area Aeff.

    v_i is the interior velocity component directed into the pipe (away from the boundary).
    out: [mdot into the pipe (signed, kg/s), momentum flux p_b + rho_b v_b^2, stagnation temperature of
    the stream, boundary pressure]. Returns the number of iterations.
    """
    c_i = np.sqrt(g * p_i / rho_i)
    beta = (g - 1.0) / (2.0 * g)
    base = 1.0 - 0.5 * (g - 1.0) * v_i / c_i
    if base < base_min:
        base = base_min
    p_w = p_i * base ** (1.0 / beta)          # wall pressure (v_b = 0)
    out[0] = 0.0
    out[1] = p_w
    out[2] = TV
    out[3] = p_w
    if Aeff <= 0.0 or pV == p_w:
        return 0
    A = Aeff if Aeff < Ap else Ap
    res = np.empty(3)
    inflow = pV > p_w
    if inflow:
        lo, hi = p_w, pV                      # residual > 0 at lo, < 0 at hi
    else:
        # sonic exit limit: w = c_b  ->  (p_b / p_i)^beta = (2 - (g - 1) v_i / c_i) / (g + 1)
        r_s = (2.0 - (g - 1.0) * v_i / c_i) / (g + 1.0)
        p_s = p_i * r_s ** (1.0 / beta) if r_s > 0.0 else 0.0
        lo = pV if pV > p_s else p_s          # residual <= 0 at lo unless the exit is choked
        hi = p_w                              # residual > 0 at hi
        if lo >= hi:
            return 0
    flo = _bc_residual(inflow, lo, A, Ap, rho_i, v_i, p_i, c_i, g, R, pV, TV, gV, RV, res)
    fhi = _bc_residual(inflow, hi, A, Ap, rho_i, v_i, p_i, c_i, g, R, pV, TV, gV, RV, res)
    pb = lo
    it = 0
    if flo * fhi < 0.0:
        # Illinois-modified regula falsi on a guaranteed bracket
        side = 0
        scale = tol * p_i
        while it < maxit:
            it += 1
            pb = (lo * fhi - hi * flo) / (fhi - flo)
            if not (lo < pb < hi):
                pb = 0.5 * (lo + hi)
            fb = _bc_residual(inflow, pb, A, Ap, rho_i, v_i, p_i, c_i, g, R, pV, TV, gV, RV, res)
            if fb == 0.0:
                break
            if fb * flo > 0.0:
                lo, flo = pb, fb
                if side == -1:
                    fhi *= 0.5
                side = -1
            else:
                hi, fhi = pb, fb
                if side == 1:
                    flo *= 0.5
                side = 1
            if hi - lo < scale:
                pb = 0.5 * (lo + hi)
                break
    elif inflow:
        pb = hi if fhi >= 0.0 else lo
    else:
        pb = lo if flo >= 0.0 else hi        # choked (sonic) exit sits at the lower bound
    _bc_residual(inflow, pb, A, Ap, rho_i, v_i, p_i, c_i, g, R, pV, TV, gV, RV, res)
    md = res[0] * Ap * res[1]
    out[0] = md if inflow else -md
    out[1] = pb + res[0] * res[1] * res[1]
    out[2] = res[2]
    out[3] = pb
    return it


# ---------------------------------------------------------------------------------------------
# right-hand side
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def _tab(tab, theta, dth):
    """Periodic linear interpolation of a crank-angle table spanning one 720 deg cycle."""
    n = tab.shape[0]
    s = (theta % FOUR_PI) / dth
    i = int(s)
    f = s - i
    if i >= n:
        i -= n
    j = i + 1
    if j >= n:
        j = 0
    return tab[i] + f * (tab[j] - tab[i])


@njit(cache=True)
def rhs(U, S, Tg, theta, omega, pipes, ends, slots, vols, tabs, props, hp, nu, dU, dS, W, sl, vst):
    """Time derivatives dU (cells) and dS (volumes). Returns the largest (|u| + c) / dx."""
    p_off, p_dx, p_area, p_diam, p_gam, p_R, p_Tw, p_fmul, p_hmul = pipes
    e_vol, e_slot = ends
    sl_const, sl_tab, sl_cyl = slots
    v_kind, v_vol, v_gam, v_R, v_p, v_T, v_Y, v_phase, v_act = vols
    tabV, tabdV, tabAi, tabAe = tabs
    Fz, T0f, dTf, R_F, R_P = props
    lim = int(nu[NU_LIMITER])
    pos = nu[NU_POS_FRAC]
    dth = nu[NU_TAB_DTHETA]
    nv = v_kind.shape[0]
    npipe = p_dx.shape[0]
    f = np.empty(NVAR)
    bc = np.empty(4)

    # ---- volume thermodynamic states: vst[v] = (p, T, gamma, R, Y, h, dV/dt, active) -------------
    for v in range(nv):
        for k in range(NS):
            dS[v, k] = 0.0
        kind = v_kind[v]
        vst[v, 7] = 1.0
        vst[v, 6] = 0.0
        if kind == VK_AMBIENT:
            vst[v, 0] = v_p[v]
            vst[v, 1] = v_T[v]
            vst[v, 2] = v_gam[v]
            vst[v, 3] = v_R[v]
            vst[v, 4] = v_Y[v]
            vst[v, 5] = 0.0
        elif kind == VK_PLENUM:
            g = v_gam[v]
            R = v_R[v]
            T = S[v, S_E] * (g - 1.0) / (S[v, S_M] * R)
            vst[v, 0] = S[v, S_M] * R * T / v_vol[v]
            vst[v, 1] = T
            vst[v, 2] = g
            vst[v, 3] = R
            vst[v, 4] = S[v, S_MY] / S[v, S_M]
            vst[v, 5] = 0.0
        else:
            if v_act[v] == 0:
                vst[v, 7] = 0.0
                continue
            m = S[v, S_M]
            x = S[v, S_MY] / m
            if x < 0.0:
                x = 0.0
            elif x > 1.0:
                x = 1.0
            T = T_from_u_frozen(S[v, S_E] / m, x, Tg[v], Fz, T0f, dTf, R_F, R_P, nu[NU_NEWTON_TOL],
                                int(nu[NU_NEWTON_MAXIT]))
            Tg[v] = T
            h, cp, R = frozen(T, x, Fz, T0f, dTf, R_F, R_P)
            th = theta - v_phase[v]
            V = _tab(tabV, th, dth)
            dV = _tab(tabdV, th, dth) * omega
            p = m * R * T / V
            vst[v, 0] = p
            vst[v, 1] = T
            vst[v, 2] = cp / (cp - R)
            vst[v, 3] = R
            vst[v, 4] = x
            vst[v, 5] = h
            vst[v, 6] = dV
            qdot = 0.0
            if hp[HP_ON] > 0.0:
                hc = coefficient(p, T, V, False, hp)
                Aw, ATw = wall_terms(V, hp)
                qdot = hc * (Aw * T - ATw)
            dS[v, S_E] = -p * dV - qdot
            dS[v, S_W] = p * dV
            dS[v, S_Q] = qdot

    smax = 0.0
    for ip in range(npipe):
        i0 = p_off[ip]
        i1 = p_off[ip + 1]
        n = i1 - i0
        g = p_gam[ip]
        R = p_R[ip]
        dx = p_dx[ip]
        A = p_area[ip]
        D = p_diam[ip]
        cp = g * R / (g - 1.0)
        # primitives
        for i in range(i0, i1):
            r = U[i, 0]
            u = U[i, 1] / r
            p = (g - 1.0) * (U[i, 2] - 0.5 * r * u * u)
            W[i, 0] = r
            W[i, 1] = u
            W[i, 2] = p
            W[i, 3] = U[i, 3] / r
            a = (abs(u) + np.sqrt(g * p / r)) / dx
            if a > smax:
                smax = a
        # limited slopes (first order in the end cells and where positivity would be at risk)
        for i in range(i0, i1):
            if i == i0 or i == i1 - 1:
                for k in range(NVAR):
                    sl[i, k] = 0.0
                continue
            for k in range(NVAR):
                sl[i, k] = _slope(W[i, k] - W[i - 1, k], W[i + 1, k] - W[i, k], lim)
            if (0.5 * abs(sl[i, 0]) > (1.0 - pos) * W[i, 0]) or (0.5 * abs(sl[i, 2]) > (1.0 - pos) * W[i, 2]):
                for k in range(NVAR):
                    sl[i, k] = 0.0
        for i in range(i0, i1):
            for k in range(NVAR):
                dU[i, k] = 0.0
        # interior faces
        for i in range(i0, i1 - 1):
            j = i + 1
            hllc(W[i, 0] + 0.5 * sl[i, 0], W[i, 1] + 0.5 * sl[i, 1], W[i, 2] + 0.5 * sl[i, 2],
                 W[i, 3] + 0.5 * sl[i, 3],
                 W[j, 0] - 0.5 * sl[j, 0], W[j, 1] - 0.5 * sl[j, 1], W[j, 2] - 0.5 * sl[j, 2],
                 W[j, 3] - 0.5 * sl[j, 3], g, f)
            for k in range(NVAR):
                dU[i, k] -= f[k] / dx
                dU[j, k] += f[k] / dx
        # ends: side 0 = left (into the pipe is +x), side 1 = right (into the pipe is -x)
        for side in range(2):
            i = i0 if side == 0 else i1 - 1
            sgn = 1.0 if side == 0 else -1.0
            v = e_vol[ip, side]
            Aeff = 0.0
            if v >= 0 and vst[v, 7] > 0.0:
                s = e_slot[ip, side]
                kind = sl_tab[s]
                if kind == SLOT_CONST:
                    Aeff = sl_const[s]
                else:
                    th = theta - v_phase[sl_cyl[s]]
                    Aeff = _tab(tabAi, th, dth) if kind == SLOT_INTAKE else _tab(tabAe, th, dth)
            vv = v if v >= 0 else 0
            end_bc(Aeff, A, W[i, 0], sgn * W[i, 1], W[i, 2], g, R, vst[vv, 0], vst[vv, 1], vst[vv, 2],
                   vst[vv, 3], nu[NU_BC_TOL], int(nu[NU_BC_MAXIT]), nu[NU_CHAR_BASE_MIN], bc)
            md = bc[0]
            dU[i, 1] += sgn * bc[1] / dx
            if md != 0.0:
                if md > 0.0:                      # volume -> pipe
                    Ys = vst[v, 4]
                    Hp = md * cp * bc[2]          # energy rate received by the pipe
                    if v_kind[v] == VK_CYL:
                        Hv = md * vst[v, 5]       # cylinder loses its own (table) enthalpy
                    else:
                        Hv = Hp
                else:                             # pipe -> volume
                    Ys = W[i, 3]
                    Hp = md * cp * bc[2]
                    if v_kind[v] == VK_CYL:
                        hin, cpin, Rin = frozen(bc[2], min(max(Ys, 0.0), 1.0), Fz, T0f, dTf, R_F, R_P)
                        Hv = md * hin
                    else:
                        Hv = Hp
                dU[i, 0] += md / (A * dx)
                dU[i, 2] += Hp / (A * dx)
                dU[i, 3] += md * Ys / (A * dx)
                if v_kind[v] != VK_AMBIENT:
                    dS[v, S_M] -= md
                    dS[v, S_MY] -= md * Ys
                    dS[v, S_E] -= Hv
                s = e_slot[ip, side]
                if sl_tab[s] == SLOT_INTAKE:
                    dS[v, S_FRESH] -= md * (1.0 - Ys)     # net fresh charge through the intake valves
                    dS[v, S_NET] -= md
                elif sl_tab[s] == SLOT_EXHAUST:
                    dS[v, S_EXH] += md                    # net mass through the exhaust valves
                else:
                    dS[v, S_NET] -= md
        # wall friction and heat transfer
        fm = p_fmul[ip]
        hm = p_hmul[ip]
        if fm > 0.0 or hm > 0.0:
            Tw = p_Tw[ip]
            for i in range(i0, i1):
                r = W[i, 0]
                u = W[i, 1]
                au = abs(u)
                T = W[i, 2] / (r * R)
                mu = nu[NU_MU_REF] * (T / nu[NU_T_SUTH_REF]) ** 1.5 * (nu[NU_T_SUTH_REF] + nu[NU_S_SUTH]) / (
                    T + nu[NU_S_SUTH])
                Re = r * au * D / mu
                if Re < nu[NU_RE_MIN]:
                    Re = nu[NU_RE_MIN]
                fl = nu[NU_F_LAM] / Re
                ft = nu[NU_F_BLASIUS] * Re ** nu[NU_F_BLASIUS_EXP]
                fd = fl if fl > ft else ft
                dU[i, 1] -= fm * fd * r * u * au / (2.0 * D)
                # Reynolds / Chilton-Colburn analogy: h = (f / 8) Pr^(-2/3) rho |u| cp, wetted area 4 / D per volume
                dU[i, 2] += hm * fd / 8.0 * nu[NU_PR] ** (-2.0 / 3.0) * r * au * cp * (Tw - T) * 4.0 / D
    return smax


@njit(cache=True)
def _floor(U, pipes, nu):
    """Density / pressure floors after a stage. Returns the number of cells that needed them."""
    p_off, p_dx, p_area, p_diam, p_gam, p_R, p_Tw, p_fmul, p_hmul = pipes
    nbad = 0
    for ip in range(p_dx.shape[0]):
        g = p_gam[ip]
        for i in range(p_off[ip], p_off[ip + 1]):
            bad = False
            if not (U[i, 0] > nu[NU_RHO_FLOOR]):
                U[i, 0] = nu[NU_RHO_FLOOR]
                U[i, 1] = 0.0
                bad = True
            ke = 0.5 * U[i, 1] * U[i, 1] / U[i, 0]
            if not ((g - 1.0) * (U[i, 2] - ke) > nu[NU_P_FLOOR]):
                U[i, 2] = nu[NU_P_FLOOR] / (g - 1.0) + ke
                bad = True
            if bad:
                nbad += 1
    return nbad


@njit(cache=True)
def advance(U, S, Tg, theta0, t_span, omega, pipes, ends, slots, vols, tabs, props, hp, nu, trace, trace_dth):
    """March the network over t_span seconds (crank angle theta = theta0 + omega t; omega = 0 for
    non-engine problems). U, S, Tg are updated in place. `trace` (nbin, n_volumes) receives the volume
    pressures binned by crank angle (pass a zero-row array to skip).
    Returns (number of steps, number of floored cells)."""
    nc = U.shape[0]
    nv = S.shape[0]
    dU = np.empty((nc, NVAR))
    dS = np.empty((nv, NS))
    U1 = np.empty((nc, NVAR))
    S1 = np.empty((nv, NS))
    W = np.empty((nc, NVAR))
    sl = np.empty((nc, NVAR))
    vst = np.zeros((nv, 8))
    t = 0.0
    nstep = 0
    nbad = 0
    nbin = trace.shape[0]
    tiny = 1.0e-14 * (t_span if t_span > 0.0 else 1.0)
    while t_span - t > tiny:
        theta = theta0 + omega * t
        smax = rhs(U, S, Tg, theta, omega, pipes, ends, slots, vols, tabs, props, hp, nu, dU, dS, W, sl, vst)
        if nbin > 0:
            b = int((theta % FOUR_PI) / trace_dth)
            if b >= nbin:
                b = nbin - 1
            for v in range(nv):
                trace[b, v] = vst[v, 0] if vst[v, 7] > 0.0 else np.nan
        dt = nu[NU_CFL] / smax if smax > 0.0 else nu[NU_DT_MAX]
        if dt > nu[NU_DT_MAX]:
            dt = nu[NU_DT_MAX]
        if dt > t_span - t:
            dt = t_span - t
        for i in range(nc):
            for k in range(NVAR):
                U1[i, k] = U[i, k] + dt * dU[i, k]
        for v in range(nv):
            for k in range(NS):
                S1[v, k] = S[v, k] + dt * dS[v, k]
        nbad += _floor(U1, pipes, nu)
        rhs(U1, S1, Tg, theta + omega * dt, omega, pipes, ends, slots, vols, tabs, props, hp, nu, dU, dS, W, sl,
            vst)
        for i in range(nc):
            for k in range(NVAR):
                U[i, k] = 0.5 * (U[i, k] + U1[i, k] + dt * dU[i, k])
        for v in range(nv):
            for k in range(NS):
                S[v, k] = 0.5 * (S[v, k] + S1[v, k] + dt * dS[v, k])
        nbad += _floor(U, pipes, nu)
        t += dt
        nstep += 1
    return nstep, nbad
