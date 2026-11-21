"""Build a pipe / volume network for the kernels in gasdyn.fv and march it.

Generic (no engine knowledge): the analytic tests (shock tube, nozzle, organ pipe, closed-system mass
conservation) use this class directly; gasdyn.engine builds the intake and exhaust systems with it.
"""
from __future__ import annotations

import math

import numpy as np

from gasdyn import fv
from gasdyn.spec import gasdyn_spec
from thermo.woschni import N_HP

DEG = math.pi / 180.0
LIMITERS = {"vanleer": fv.LIM_VANLEER, "minmod": fv.LIM_MINMOD}
WALL = None


def numerics_vector(gspec: dict, tol: float = 1.0e-10, maxit: int = 50) -> np.ndarray:
    """Pack spec/gasdyn.yaml numerics, transport and pipe-friction constants for the numba kernels."""
    n, tr, pf = gspec["numerics"], gspec["transport"], gspec["pipe_friction"]
    nu = np.zeros(fv.N_NU)
    nu[fv.NU_CFL] = n["cfl"]
    nu[fv.NU_RHO_FLOOR] = n["rho_floor"]
    nu[fv.NU_P_FLOOR] = n["p_floor"]
    nu[fv.NU_LIMITER] = LIMITERS[n["limiter"]]
    nu[fv.NU_POS_FRAC] = n["positivity_fraction"]
    nu[fv.NU_MU_REF] = tr["mu_ref"]
    nu[fv.NU_T_SUTH_REF] = tr["T_ref"]
    nu[fv.NU_S_SUTH] = tr["sutherland"]
    nu[fv.NU_PR] = tr["prandtl"]
    nu[fv.NU_F_LAM] = pf["laminar_coeff"]
    nu[fv.NU_F_BLASIUS] = pf["blasius_coeff"]
    nu[fv.NU_F_BLASIUS_EXP] = pf["blasius_exp"]
    nu[fv.NU_RE_MIN] = pf["re_min"]
    nu[fv.NU_BC_TOL] = n["bc_tol"]
    nu[fv.NU_BC_MAXIT] = n["bc_max_iter"]
    nu[fv.NU_DT_MAX] = n["dt_max"]
    nu[fv.NU_NEWTON_TOL] = tol
    nu[fv.NU_NEWTON_MAXIT] = maxit
    nu[fv.NU_TAB_DTHETA] = n["table_step_deg"] * DEG
    nu[fv.NU_CHAR_BASE_MIN] = n["char_base_min"]
    return nu


class Network:
    """Pipes joined to 0D volumes through orifices. Build with add_*, then finalize()."""

    def __init__(self, gspec: dict | None = None, nu: np.ndarray | None = None):
        self.gspec = gspec or gasdyn_spec()
        self.nu = numerics_vector(self.gspec) if nu is None else nu.copy()
        self._vols: list[dict] = []
        self._pipes: list[dict] = []
        self._slots: list[tuple] = []     # (constant area, table kind, cylinder volume index)
        self.tabs = tuple(np.zeros(2) for _ in range(4))
        self.props = (np.zeros((4, 2)), 0.0, 1.0, 1.0, 1.0)
        self.hp = np.zeros(N_HP)
        self.final = False

    # ---- volumes ------------------------------------------------------------------------------
    def add_ambient(self, p, T, gamma, R, Y=0.0) -> int:
        self._vols.append(dict(kind=fv.VK_AMBIENT, vol=0.0, g=gamma, R=R, p=p, T=T, Y=Y, phase=0.0))
        return len(self._vols) - 1

    def add_plenum(self, volume, p, T, gamma, R, Y=0.0) -> int:
        self._vols.append(dict(kind=fv.VK_PLENUM, vol=volume, g=gamma, R=R, p=p, T=T, Y=Y, phase=0.0))
        return len(self._vols) - 1

    def add_cylinder(self, phase_rad: float) -> int:
        """Engine cylinder; local crank angle = theta - phase. State is set by gasdyn.engine."""
        self._vols.append(dict(kind=fv.VK_CYL, vol=0.0, g=0.0, R=0.0, p=0.0, T=0.0, Y=0.0, phase=phase_rad))
        return len(self._vols) - 1

    # ---- pipes --------------------------------------------------------------------------------
    def _end(self, end):
        """end: WALL, (volume, area) for a fixed orifice / open end, or (cylinder, 'intake' | 'exhaust')."""
        if end is WALL:
            return -1, 0
        vol, area = end
        if isinstance(area, str):
            kind = {"intake": fv.SLOT_INTAKE, "exhaust": fv.SLOT_EXHAUST}[area]
            self._slots.append((0.0, kind, vol))
        else:
            self._slots.append((float(area), fv.SLOT_CONST, 0))
        return vol, len(self._slots) - 1

    def add_pipe(self, length, diameter, gamma, R, left=WALL, right=WALL, n_cells=None, dx=None, wall_T=300.0,
                 friction_multiplier=0.0, heat_multiplier=0.0, p=1.0e5, T=300.0, u=0.0, Y=0.0) -> int:
        if n_cells is None:
            dx = dx or self.gspec["numerics"]["dx"]
            n_cells = max(int(self.gspec["numerics"]["min_cells"]), int(round(length / dx)))
        lv, ls = self._end(left)
        rv, rs = self._end(right)
        self._pipes.append(dict(L=length, D=diameter, g=gamma, R=R, n=int(n_cells), Tw=wall_T,
                                fm=friction_multiplier, hm=heat_multiplier, ev=(lv, rv), es=(ls, rs),
                                init=(p, T, u, Y)))
        return len(self._pipes) - 1

    @staticmethod
    def area(diameter):
        return 0.25 * math.pi * diameter**2

    # ---- assembly -----------------------------------------------------------------------------
    def finalize(self):
        P, V = self._pipes, self._vols
        if not self._slots:
            self._slots.append((0.0, fv.SLOT_CONST, 0))
        if not V:
            self.add_ambient(1.0, 1.0, 1.4, 1.0)
        n = np.array([q["n"] for q in P], dtype=np.int64)
        self.p_off = np.concatenate([[0], np.cumsum(n)]).astype(np.int64)
        f = lambda k: np.array([q[k] for q in P], dtype=float)
        self.p_len = f("L")
        self.p_dx = self.p_len / n
        self.p_diam = f("D")
        self.p_area = 0.25 * np.pi * self.p_diam**2
        self.p_gam, self.p_R, self.p_Tw, self.p_fmul, self.p_hmul = f("g"), f("R"), f("Tw"), f("fm"), f("hm")
        self.e_vol = np.array([q["ev"] for q in P], dtype=np.int64).reshape(-1, 2)
        self.e_slot = np.array([q["es"] for q in P], dtype=np.int64).reshape(-1, 2)
        self.sl_const = np.array([s[0] for s in self._slots], dtype=float)
        self.sl_tab = np.array([s[1] for s in self._slots], dtype=np.int64)
        self.sl_cyl = np.array([s[2] for s in self._slots], dtype=np.int64)
        g = lambda k: np.array([q[k] for q in V], dtype=float)
        self.v_kind = np.array([q["kind"] for q in V], dtype=np.int64)
        self.v_vol, self.v_gam, self.v_R = g("vol"), g("g"), g("R")
        self.v_p, self.v_T, self.v_Y, self.v_phase = g("p"), g("T"), g("Y"), g("phase")
        self.v_act = np.zeros(len(V), dtype=np.int64)
        nc = int(self.p_off[-1])
        self.U = np.zeros((nc, fv.NVAR))
        self.S = np.zeros((len(V), fv.NS))
        self.Tg = np.full(len(V), 300.0)
        for ip, q in enumerate(P):
            self.set_pipe(ip, *q["init"])
        for iv, q in enumerate(V):
            if q["kind"] == fv.VK_PLENUM:
                self.set_plenum(iv, q["p"], q["T"], q["Y"])
        self.final = True
        self.nbad = 0
        self.nstep = 0
        return self

    def cells(self, ip) -> slice:
        return slice(int(self.p_off[ip]), int(self.p_off[ip + 1]))

    def x(self, ip):
        n = self.p_off[ip + 1] - self.p_off[ip]
        return (np.arange(n) + 0.5) * self.p_dx[ip]

    def set_pipe(self, ip, p, T, u=0.0, Y=0.0):
        """Uniform or per-cell (arrays) primitive state of one pipe."""
        g, R = self.p_gam[ip], self.p_R[ip]
        c = self.cells(ip)
        rho = np.asarray(p, dtype=float) / (R * np.asarray(T, dtype=float))
        self.U[c, 0] = rho
        self.U[c, 1] = rho * u
        self.U[c, 2] = np.asarray(p, dtype=float) / (g - 1.0) + 0.5 * rho * np.asarray(u, dtype=float) ** 2
        self.U[c, 3] = rho * Y

    def set_plenum(self, iv, p, T, Y=0.0):
        g, R = self.v_gam[iv], self.v_R[iv]
        m = p * self.v_vol[iv] / (R * T)
        self.S[iv, :] = 0.0
        self.S[iv, fv.S_M] = m
        self.S[iv, fv.S_MY] = m * Y
        self.S[iv, fv.S_E] = m * R * T / (g - 1.0)

    def primitives(self, ip):
        """(x, rho, u, p, Y) arrays of one pipe."""
        g = self.p_gam[ip]
        U = self.U[self.cells(ip)]
        rho = U[:, 0]
        u = U[:, 1] / rho
        return self.x(ip), rho.copy(), u, (g - 1.0) * (U[:, 2] - 0.5 * rho * u * u), U[:, 3] / rho

    def plenum_state(self, iv):
        """(p, T) of a fixed-volume plenum."""
        g, R = self.v_gam[iv], self.v_R[iv]
        T = self.S[iv, fv.S_E] * (g - 1.0) / (self.S[iv, fv.S_M] * R)
        return self.S[iv, fv.S_M] * R * T / self.v_vol[iv], T

    def total_mass(self) -> float:
        """Gas mass in all pipes, plenums and active cylinders (ambient volumes excluded)."""
        m = 0.0
        for ip in range(len(self.p_dx)):
            m += float(np.sum(self.U[self.cells(ip), 0])) * self.p_dx[ip] * self.p_area[ip]
        for iv, kind in enumerate(self.v_kind):
            if kind == fv.VK_PLENUM or (kind == fv.VK_CYL and self.v_act[iv]):
                m += self.S[iv, fv.S_M]
        return m

    def kernel_args(self):
        pipes = (self.p_off, self.p_dx, self.p_area, self.p_diam, self.p_gam, self.p_R, self.p_Tw, self.p_fmul,
                 self.p_hmul)
        ends = (self.e_vol, self.e_slot)
        slots = (self.sl_const, self.sl_tab, self.sl_cyl)
        vols = (self.v_kind, self.v_vol, self.v_gam, self.v_R, self.v_p, self.v_T, self.v_Y, self.v_phase,
                self.v_act)
        return pipes, ends, slots, vols, self.tabs, self.props, self.hp, self.nu

    _NO_TRACE = np.zeros((0, 1))

    def advance(self, t_span, theta0=0.0, omega=0.0, trace=None, trace_dth=1.0):
        """March t_span seconds; returns the number of time steps taken."""
        nstep, nbad = fv.advance(self.U, self.S, self.Tg, theta0, t_span, omega, *self.kernel_args(),
                                 self._NO_TRACE if trace is None else trace, trace_dth)
        self.nbad += nbad
        self.nstep += nstep
        return nstep
