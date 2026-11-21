"""Quasi-dimensional turbulent flame model (milestone M7), an alternative burn law to the Wiebe function.

Physics (constants and citations in spec/combustion_ext.yaml `flame`):

* Laminar burning velocity of iso-octane, Metghalchi & Keck (1982) with the Rhodes & Keck residual
  dilution factor (Heywood 1988 eqs. 9.36-9.38, Table 9.3).
* 0D k-epsilon turbulence of the unburned charge from IVC: rapid-distortion compression
  (k ~ rho_u^(2/3), eps ~ rho_u^(4/3)) and decay (dk/dt = -eps, deps/dt = -C_eps2 eps^2/k), operator-split
  per step (Strang) with the exact solution of each part. u' = sqrt(2k/3), L = k^1.5/eps.
* Blizard-Keck / Tabaczynski entrainment and burn-up:
    dm_e/dt = rho_u A_f (C_T u' (1 - exp(-r_f/L)) + S_L)
    dm_b/dt = rho_u A_f S_L + (m_e - m_b) / tau_b,   tau_b = lambda_T / S_L,  lambda_T = sqrt(15 nu u'^2/eps)
* Flame geometry: a sphere centred on the bore axis at the fire deck, clipped by the cylinder wall
  (radius B/2) and the piston (flat-disc chamber of height h = V / A_piston). Volume and front area are
  closed-form (see `sphere_cyl_volume`, `sphere_cyl_area`).

Coupling with the two-zone model (thermo.twozone) is staggered: a burn law gives a two-zone run, whose
p, T_u, T_b trajectories drive the flame ODE to a new burn law, repeated until max |delta x_b| <
coupling_tol. The first burn law is the spec Wiebe shape started at the spark plus the spec ignition
delay. All results are MODELLED (no physical engine).
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from kin.slider_crank import cylinder_volume
from thermo import twozone as tz
from thermo.cycle import _grid, run_closed
from thermo.props import T_adiabatic_p, burned, frozen

SEC_PER_MIN = 60.0

# Metghalchi-Keck parameter vector layout
(MK_PHI_M, MK_BM, MK_BPHI, MK_A0, MK_APHI, MK_B0, MK_BPHI_P, MK_T0, MK_P0, MK_CDIL, MK_EDIL, MK_PHI) = range(12)
N_MK = 12


def mk_params(cspec: dict, phi: float) -> np.ndarray:
    lam = cspec["flame"]["laminar"]
    v = np.zeros(N_MK)
    v[MK_PHI_M], v[MK_BM], v[MK_BPHI] = lam["phi_M"], lam["B_M"], lam["B_phi"]
    v[MK_A0], v[MK_APHI], v[MK_B0], v[MK_BPHI_P] = lam["alpha_0"], lam["alpha_phi"], lam["beta_0"], lam["beta_phi"]
    v[MK_T0], v[MK_P0] = lam["T0"], lam["p0"]
    v[MK_CDIL], v[MK_EDIL] = lam["dilution_coeff"], lam["dilution_exp"]
    v[MK_PHI] = phi
    return v


@njit(cache=True)
def laminar_speed(Tu, p, x_dil, mk):
    """Metghalchi-Keck laminar burning velocity [m/s] with residual dilution x_dil (mass fraction)."""
    phi = mk[MK_PHI]
    s0 = mk[MK_BM] + mk[MK_BPHI] * (phi - mk[MK_PHI_M]) ** 2
    a = mk[MK_A0] + mk[MK_APHI] * (phi - 1.0)
    b = mk[MK_B0] + mk[MK_BPHI_P] * (phi - 1.0)
    f = 1.0 - mk[MK_CDIL] * x_dil ** mk[MK_EDIL] if x_dil > 0.0 else 1.0
    return s0 * (Tu / mk[MK_T0]) ** a * (p / mk[MK_P0]) ** b * f


# ---------------------------------------------------------------------------------------------
# flame geometry: sphere of radius r centred on the axis at the fire deck, inside a cylinder of radius
# R and height h. Slices z in [0, min(r, h)] are discs of radius min(sqrt(r^2 - z^2), R); the spherical
# front lies where sqrt(r^2 - z^2) <= R, i.e. z >= z0 = sqrt(max(r^2 - R^2, 0)), with dA = 2 pi r dz.
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def sphere_cyl_volume(r, h, R):
    z0 = math.sqrt(max(r * r - R * R, 0.0))
    zt = min(r, h)
    za = min(z0, zt)
    return math.pi * R * R * za + math.pi * ((r * r * zt - zt**3 / 3.0) - (r * r * za - za**3 / 3.0))


@njit(cache=True)
def sphere_cyl_area(r, h, R):
    z0 = math.sqrt(max(r * r - R * R, 0.0))
    zt = min(r, h)
    return 2.0 * math.pi * r * max(zt - z0, 0.0)


@njit(cache=True)
def flame_radius(Vf, h, R):
    """Invert sphere_cyl_volume(r, h, R) = Vf by bisection (monotone in r)."""
    lo, hi = 0.0, math.sqrt(R * R + h * h)
    if Vf >= sphere_cyl_volume(hi, h, R):
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if sphere_cyl_volume(mid, h, R) < Vf:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------------------------
# turbulence
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def _decay(k, e, dt, c2):
    s = 1.0 + (c2 - 1.0) * dt * e / k
    return k * s ** (-1.0 / (c2 - 1.0)), e * s ** (-c2 / (c2 - 1.0))


@njit(cache=True)
def turbulence(rho_u, dt, k0, e0, c2):
    """k, eps on the given grid (uniform time step dt) from IVC, Strang split decay / RDT / decay."""
    n = rho_u.size
    k = np.empty(n)
    e = np.empty(n)
    k[0], e[0] = k0, e0
    for i in range(n - 1):
        kk, ee = _decay(k[i], e[i], 0.5 * dt, c2)
        r = rho_u[i + 1] / rho_u[i]
        kk *= r ** (2.0 / 3.0)
        ee *= r ** (4.0 / 3.0)
        k[i + 1], e[i + 1] = _decay(kk, ee, 0.5 * dt, c2)
    return k, e


# ---------------------------------------------------------------------------------------------
# entrainment / burn-up ODE
# ---------------------------------------------------------------------------------------------
@njit(cache=True)
def _rhs(me, mb, m, rho_u, rho_b, V, Tu, p, k, e, nu, xr, Ap, R, CT, dev, mk, out):
    h = V / Ap
    Vf = (me - mb) / rho_u + mb / rho_b
    if Vf > V:
        Vf = V
    r = flame_radius(Vf, h, R)
    A = sphere_cyl_area(r, h, R)
    sl = laminar_speed(Tu, p, xr, mk)
    up = math.sqrt(2.0 * k / 3.0)
    L = k**1.5 / e
    lam_t = math.sqrt(15.0 * nu * up * up / e)
    if me < m:
        fdev = 1.0 - math.exp(-r / (dev * L)) if dev > 0.0 else 1.0
        dme = rho_u * A * (CT * up * fdev + sl)
        dmb = rho_u * A * sl + (me - mb) * sl / lam_t
    else:
        dme = 0.0
        dmb = (me - mb) * sl / lam_t
    out[0] = dme
    out[1] = dmb
    out[2] = r
    out[3] = sl
    out[4] = up
    return out


@njit(cache=True)
def integrate_flame(th_deg, theta_spark, omega, m, xr, rho_u, rho_b, V, Tu, p, k, e, nu, Ap, R, CT, dev, mk,
                    me0):
    """Integrate the flame ODE from the (continuous) spark angle along the half grid th_deg.

    Inputs are sampled on th_deg; stage values between grid points are linear interpolations.
    Returns x_b, dx_b/dtheta [1/rad], flame radius [m], S_L [m/s], u' [m/s] on th_deg and the index of
    the first grid point after the spark."""
    n = th_deg.size
    xb = np.zeros(n)
    dxb = np.zeros(n)
    rf = np.zeros(n)
    slv = np.zeros(n)
    upv = np.zeros(n)
    i0 = 0
    while i0 < n and th_deg[i0] < theta_spark:
        i0 += 1
    if i0 >= n:
        return xb, dxb, rf, slv, upv, i0
    buf = np.empty(5)
    d1 = np.empty(2)
    d2 = np.empty(2)
    d3 = np.empty(2)
    d4 = np.empty(2)
    deg = math.pi / 180.0

    def lerp(a, j, f):
        return a[j] + f * (a[j + 1] - a[j])

    me = me0
    mb = me0
    # march from the spark (between grid points j-1 and j) to each next grid point
    for i in range(i0, n):
        if i == i0:
            ta = theta_spark
            ja = max(i0 - 1, 0)
        else:
            ta = th_deg[i - 1]
            ja = i - 1
        tb = th_deg[i]
        dt = (tb - ta) * deg / omega
        if i0 == 0:
            fa = 0.0
        else:
            fa = (ta - th_deg[ja]) / (th_deg[ja + 1] - th_deg[ja])
        fm = 0.5 * (fa + 1.0)
        if dt > 0.0:
            for st in range(4):
                if st == 0:
                    f, y0, y1 = fa, me, mb
                elif st == 1:
                    f, y0, y1 = fm, me + 0.5 * dt * d1[0], mb + 0.5 * dt * d1[1]
                elif st == 2:
                    f, y0, y1 = fm, me + 0.5 * dt * d2[0], mb + 0.5 * dt * d2[1]
                else:
                    f, y0, y1 = 1.0, me + dt * d3[0], mb + dt * d3[1]
                y0 = min(y0, m)
                y1 = min(y1, y0)
                _rhs(y0, y1, m, lerp(rho_u, ja, f), lerp(rho_b, ja, f), lerp(V, ja, f), lerp(Tu, ja, f),
                     lerp(p, ja, f), lerp(k, ja, f), lerp(e, ja, f), lerp(nu, ja, f), xr, Ap, R, CT, dev, mk, buf)
                d = d1 if st == 0 else (d2 if st == 1 else (d3 if st == 2 else d4))
                d[0] = buf[0]
                d[1] = buf[1]
            me = min(me + dt / 6.0 * (d1[0] + 2.0 * d2[0] + 2.0 * d3[0] + d4[0]), m)
            mb = min(mb + dt / 6.0 * (d1[1] + 2.0 * d2[1] + 2.0 * d3[1] + d4[1]), me)
        _rhs(me, mb, m, rho_u[i], rho_b[i], V[i], Tu[i], p[i], k[i], e[i], nu[i], xr, Ap, R, CT, dev, mk, buf)
        xb[i] = mb / m
        dxb[i] = buf[1] / (omega * m)
        rf[i] = buf[2]
        slv[i] = buf[3]
        upv[i] = buf[4]
    return xb, dxb, rf, slv, upv, i0


@njit(cache=True)
def _burned_density(p, Tu, Tb, xr, Fz, T0f, dTf, R_F, R_P, B, T0b, dTb, lp0, dlp, tol, maxit, Tguess):
    """Burned density on the grid; where the two-zone run has no burned zone yet (NaN T_b) the
    constant-pressure adiabatic flame temperature of the local unburned gas is used."""
    n = p.size
    rb = np.empty(n)
    for i in range(n):
        T = Tb[i]
        if not np.isfinite(T):
            hu, cpu, Ru = frozen(Tu[i], xr, Fz, T0f, dTf, R_F, R_P)
            T = T_adiabatic_p(hu, p[i], Tguess, B, T0b, dTb, lp0, dlp, tol, maxit)
        h, R, hT, hp, RT, Rp = burned(T, p[i], B, T0b, dTb, lp0, dlp)
        rb[i] = p[i] / (R * T)
    return rb


def _to_half(full):
    """Full-grid (N+1) samples -> half grid (2N+1) by midpoint averaging."""
    out = np.empty(2 * full.size - 1)
    out[0::2] = full
    out[1::2] = 0.5 * (full[:-1] + full[1:])
    return out


class FlameModel:
    """Spark-timed quasi-dimensional burn law coupled to the two-zone closed cycle."""

    def __init__(self, spec: dict, tspec: dict, cspec: dict, gspec: dict, tables, geom):
        self.spec, self.tspec, self.cspec, self.tables, self.geom = spec, tspec, cspec, tables, geom
        fl = cspec["flame"]
        self.fl = fl
        self.mk = mk_params(cspec, 1.0 / float(tspec["mixture"]["lambda"]))
        self.tr = gspec["transport"]
        tu = fl["turbulence"]
        self.c2 = float(tu["C_eps2"])
        self.u_ivc_per_sp = float(tu["u_ivc_per_sp"])
        self.L_ivc = float(tu["L_ivc_per_lift"]) * float(spec["valves"]["intake_max_lift"])
        self.CT = float(fl["entrainment_coeff"])
        self.dev = float(fl["kernel_development_scale"])
        self.r_kernel = float(fl["kernel_radius"])
        if fl["plug"] != "centre":
            raise NotImplementedError(fl["plug"])
        comb = spec["combustion"]
        self.ign_delay_deg = float(comb["spark_advance_deg"]) + float(comb["heat_release_start_deg"])

    def viscosity(self, T):
        """Sutherland dynamic viscosity [Pa s] (spec/gasdyn.yaml transport)."""
        t = self.tr
        return t["mu_ref"] * (T / t["T_ref"]) ** 1.5 * (t["T_ref"] + t["sutherland"]) / (T + t["sutherland"])

    def closed(self, m, xr, T0, theta_s, theta_e, step_deg, hp, rpm, spark_adv_deg, num):
        """Converged flame-model closed part. Returns (ClosedResult, info dict)."""
        tables, geom = self.tables, self.geom
        comb = self.spec["combustion"]
        theta_spark = -float(spark_adv_deg)
        n, hdeg, half = _grid(theta_s, theta_e, step_deg)
        Vh = cylinder_volume(np.radians(half), geom)
        omega = rpm * 2.0 * math.pi / SEC_PER_MIN
        sp = 2.0 * geom.stroke * rpm / SEC_PER_MIN
        k0 = 1.5 * (self.u_ivc_per_sp * sp) ** 2      # isotropic: k = 3/2 u'^2
        e0 = k0**1.5 / self.L_ivc
        targs = tables.args()
        Ru = frozen(T0, xr, *targs[:5])[2]
        Ap = geom.piston_area
        Rcyl = 0.5 * geom.bore
        tol, maxit = float(num["newton_tol"]), int(num["newton_max_iter"])
        cl = run_closed(tables, geom, m, xr, T0, theta_s, theta_e, step_deg, hp, "wiebe",
                        theta_spark + self.ign_delay_deg, comb["burn_duration_deg"], comb["wiebe_a"],
                        comb["wiebe_m"], num=num)
        xb_prev = cl.col(tz.O_XB).copy()
        hstep = math.radians(hdeg) / 2.0 / omega
        dxmax = float("inf")
        it = 0
        self.trace_dx = []
        for it in range(int(self.fl["coupling_max_iter"])):
            p = _to_half(cl.col(tz.O_P))
            Tu = _to_half(cl.col(tz.O_TU))
            Tb_full = cl.col(tz.O_TB)
            rho_b = _to_half(_burned_density(cl.col(tz.O_P), cl.col(tz.O_TU), Tb_full, xr, *targs, tol, maxit,
                                             float(num["T_flame_guess"])))
            rho_u = p / (Ru * Tu)
            k, e = turbulence(rho_u, hstep, k0, e0, self.c2)
            nu = self.viscosity(Tu) / rho_u
            # ignition kernel: unburned-density mass of the clipped sphere of kernel radius at the spark
            j = int(np.searchsorted(half, theta_spark))
            j = min(max(j, 0), half.size - 1)
            me0 = rho_u[j] * sphere_cyl_volume(self.r_kernel, Vh[j] / Ap, Rcyl)
            xb, dxb, rf, sl, up, i0 = integrate_flame(half, theta_spark, omega, m, xr, rho_u, rho_b, Vh, Tu, p, k,
                                                      e, nu, Ap, Rcyl, self.CT, self.dev, self.mk, me0)
            if i0 >= half.size:
                raise ValueError("spark after the end of the closed part")
            k_act = (i0 + 1) // 2
            cl = run_closed(tables, geom, m, xr, T0, theta_s, theta_e, step_deg, hp, "profile",
                            xb_half=xb, dxb_half=dxb, act_index=k_act, num=num)
            xb_new = cl.col(tz.O_XB)
            dxmax = float(np.max(np.abs(xb_new - xb_prev)))
            self.trace_dx.append(dxmax)
            xb_prev = xb_new.copy()
            if dxmax < float(self.fl["coupling_tol"]):
                break
        th = half[::2]
        info = {"xb_flame": xb[::2], "k_act": k_act, "i0": i0, "coupling_iter": it + 1, "coupling_dx": dxmax,
                "u_prime": up[::2], "S_L": sl[::2], "r_flame": rf[::2], "theta": th, "k": k[::2], "eps": e[::2],
                "u_prime_tdc": float(np.interp(0.0, half, np.sqrt(2.0 * k / 3.0))),
                "S_L_spark": float(sl[i0]) if i0 < sl.size else float("nan")}
        return cl, info
