"""Big-end and main bearing loads, and a dynamically loaded short-bearing (Ocvirk) film model.

Frames: the cylinder plane is described in (z, y) coordinates (z up the bore, y lateral, ADR 0001)
with angles measured from +z towards +y, which is the crank rotation direction. A throw at local
angle t sits at r (cos t, sin t).

Loads
- Force of the rod on the crankpin (big end):  F_rp = -F_rod (cos b, -sin b) + m_be r w^2 (cos t, sin t)
  where (cos b, -sin b) points from the crankpin to the small end and m_be is the rod big-end mass.
- Crankpin frame: radial (outward along the throw) and tangential (rotation direction) components.
- Rod frame: the load of the pin on the big-end shell (-F_rp) along the rod axis (towards the small
  end) and transverse to it.
- Main bearings: statically determinate split. Each throw (crankpin force + net throw unbalance
  after counterweights) is carried half by each adjacent main journal (5 mains, 4 throws).

Film model (short bearing, Ocvirk, half-Sommerfeld pi-film cavitation), solved for the journal
centre velocity the way the mobility method does (Booker, J. Lubr. Tech. 1965, "Dynamically
loaded journal bearings: mobility method of solution"): for film thickness h = c (1 - eps cos psi),
psi from the eccentricity direction, the short-bearing Reynolds equation gives the film force

    F = -(mu R L^3 A / c^2) * integral over the pressure zone of (-sin(psi - a)) (cos psi, sin psi) / D^3
    with D = 1 - eps cos psi,  u = eps_dot = A sin a,  v = eps (w_bar - phi_dot) = A cos a,
    w_bar = (w_journal + w_bearing) / 2,  pressure zone where sin(psi - a) < 0.

Given the load, the direction fixes the "squeeze angle" a and the magnitude fixes A, which yields
the eccentricity velocity. The orbit is integrated over several cycles. Journal inertia and oil
feed/thermal effects are neglected; this is a design check, not a validated bearing code.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from dyn import params
from dyn.forces import cylinder_forces
from dyn.pressure import PressureTrace


def _rot90(v):
    """Rotate (z, y) vectors by +90 deg in the rotation direction."""
    return np.stack([-v[1], v[0]])


@dataclass
class BigEndLoads:
    theta: np.ndarray            # local crank angle [rad]
    F_pin: np.ndarray            # (2, N) force of rod on crankpin, engine frame (z, y) [N]
    pin_frame: np.ndarray        # (2, N) (radial outward, tangential) of F_pin [N]
    rod_frame: np.ndarray        # (2, N) (axial to small end, transverse) of the load on the shell [N]
    omega_rod: np.ndarray        # absolute angular velocity of the rod [rad/s]


def big_end_loads(theta, trace: PressureTrace, omega: float, spec=None, dspec=None) -> BigEndLoads:
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    g = params.geometry(spec)
    m = params.masses(spec)
    th = np.asarray(theta, dtype=float)
    f = cylinder_forces(th, trace, omega, spec, dspec)
    b = f.beta
    rod_dir = np.stack([np.cos(b), -np.sin(b)])
    e_r = np.stack([np.cos(th), np.sin(th)])
    F_pin = -f.F_rod * rod_dir + m.rod_big_end * g.crank_radius * omega**2 * e_r
    e_t = _rot90(e_r)
    pin_frame = np.stack([(F_pin * e_r).sum(0), (F_pin * e_t).sum(0)])
    shell = -F_pin
    rod_n = _rot90(rod_dir)
    rod_frame = np.stack([(shell * rod_dir).sum(0), (shell * rod_n).sum(0)])
    lam = g.lam
    beta_dot = omega * lam * np.cos(th) / np.cos(b)
    return BigEndLoads(th, F_pin, pin_frame, rod_frame, -beta_dot)


def main_bearing_loads(theta, trace: PressureTrace, omega: float, spec=None, dspec=None) -> np.ndarray:
    """Load of the crank on each main bearing, engine frame (z, y): array (n_cyl + 1, 2, N) [N].

    `theta` is the global crank angle (cylinder 1 local angle)."""
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    g = params.geometry(spec)
    m = params.masses(spec)
    th = np.asarray(theta, dtype=float)
    cr = dspec["crank"]
    net_mr = cr["throw_unbalance_mr"] - cr["counterweight_fraction"] * (
        cr["throw_unbalance_mr"] + m.rod_big_end * g.crank_radius)
    phases = params.firing_offsets(spec)
    n = len(phases)
    loads = np.zeros((n + 1, 2, th.size))
    for k, ph in enumerate(phases):
        t_loc = th - ph
        be = big_end_loads(t_loc, trace, omega, spec, dspec)
        seg = be.F_pin + net_mr * omega**2 * np.stack([np.cos(t_loc), np.sin(t_loc)])
        loads[k] += 0.5 * seg
        loads[k + 1] += 0.5 * seg
    return loads


# ---------------------------------------------------------------------------------------------
# short-bearing film
# ---------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Bearing:
    diameter: float
    length: float
    radial_clearance: float
    viscosity: float

    @property
    def radius(self) -> float:
        return self.diameter / 2

    @property
    def force_scale(self) -> float:
        """mu R L^3 / c^2  [N s]."""
        return self.viscosity * self.radius * self.length**3 / self.radial_clearance**2


def bearing(kind: str, dspec=None) -> Bearing:
    b = (dspec or params.dyn_spec())["bearings"]
    k = b[kind]
    return Bearing(k["diameter"], k["length"], k["radial_clearance"], b["oil_viscosity"])


class ShortBearingFilm:
    """Film-force integrals I(a, eps) of the pi-film short bearing (see module docstring)."""

    def __init__(self, n_psi: int, eps_max: float, n_alpha: int = 721):
        self.psi = np.linspace(0.0, params.FULL_CYCLE, n_psi)   # 0 .. 4 pi covers a + 2 pi for a < 2 pi
        self.eps_max = eps_max
        self.alpha = np.linspace(0.0, 2 * math.pi, n_alpha)
        self._eps = None

    def _cumulative(self, eps: float):
        if self._eps == eps:
            return self._cum
        psi = self.psi
        D3 = (1.0 - eps * np.cos(psi)) ** 3
        s, c = np.sin(psi), np.cos(psi)
        h = np.diff(psi)
        cums = []
        for g in (s * c / D3, c * c / D3, s * s / D3):
            cums.append(np.concatenate([[0.0], np.cumsum(0.5 * h * (g[1:] + g[:-1]))]))
        self._eps, self._cum = eps, cums
        return cums

    def integral(self, alpha, eps: float):
        """I = (I_r, I_t) in the (eccentricity, rotation-direction) frame; alpha scalar or array."""
        eps = min(eps, self.eps_max)
        c1, c2, c3 = self._cumulative(eps)
        a = np.asarray(alpha, dtype=float)
        lo, hi = a + math.pi, a + 2 * math.pi
        C1 = np.interp(hi, self.psi, c1) - np.interp(lo, self.psi, c1)
        C2 = np.interp(hi, self.psi, c2) - np.interp(lo, self.psi, c2)
        C3 = np.interp(hi, self.psi, c3) - np.interp(lo, self.psi, c3)
        sa, ca = np.sin(a), np.cos(a)
        return np.stack([sa * C2 - ca * C1, sa * C1 - ca * C3])

    def solve(self, W_rt, eps: float):
        """Squeeze angle a and magnitude A such that the film force balances load W
        (W = -film force on journal), W given in the (e_r, e_t) frame. Returns (a, A / force_scale)."""
        wr, wt = W_rt
        I = self.integral(self.alpha, eps)
        cross = I[0] * wt - I[1] * wr
        dot = I[0] * wr + I[1] * wt
        idx = np.nonzero((np.sign(cross[:-1]) != np.sign(cross[1:])) & (dot[:-1] > 0))[0]
        if idx.size == 0:
            raise RuntimeError("no squeeze angle found for the bearing load")
        i = idx[0]

        def fcross(a):
            v = self.integral(a, eps)
            return v[0] * wt - v[1] * wr

        a = brentq(fcross, self.alpha[i], self.alpha[i + 1], xtol=1e-12)
        Imag = float(np.hypot(*self.integral(a, eps)))
        return a, math.hypot(wr, wt) / Imag


def ocvirk_static_load(eps, b: Bearing, omega: float):
    """Closed-form Ocvirk short-bearing load capacity (bearing fixed, journal speed omega):
    W = mu U L^3 / (4 c^2) * eps / (1 - eps^2)^2 * sqrt(pi^2 (1 - eps^2) + 16 eps^2), U = R omega."""
    eps = np.asarray(eps, dtype=float)
    U = b.radius * omega
    return (b.viscosity * U * b.length**3 / (4 * b.radial_clearance**2) * eps / (1 - eps**2) ** 2
            * np.sqrt(math.pi**2 * (1 - eps**2) + 16 * eps**2))


def ocvirk_attitude(eps):
    """Ocvirk attitude angle: tan(phi) = pi sqrt(1 - eps^2) / (4 eps)."""
    eps = np.asarray(eps, dtype=float)
    return np.arctan2(math.pi * np.sqrt(1 - eps**2), 4 * eps)


def static_eccentricity(W: float, b: Bearing, omega: float, film: ShortBearingFilm):
    """Steady load W on a fixed bearing: eps where the film model with eps_dot = phi_dot = 0
    carries W (a = 0, A = eps * omega / 2). Returns (eps, attitude angle)."""
    wbar = omega / 2

    def resid(eps):
        I0 = film.integral(0.0, eps)
        return b.force_scale * eps * wbar * float(np.hypot(*I0)) - W

    eps = brentq(resid, 1e-9, film.eps_max, xtol=1e-14)
    I0 = film.integral(0.0, eps)
    # load direction relative to the eccentricity direction is angle(I0); attitude = -that angle
    return eps, -math.atan2(I0[1], I0[0])


@dataclass
class Orbit:
    theta: np.ndarray       # crank angle of the reported (last) cycle [rad]
    eps: np.ndarray         # (2, N) eccentricity ratio vector, engine frame
    h_min: np.ndarray       # minimum film thickness [m]

    @property
    def min_film(self) -> float:
        return float(self.h_min.min())

    @property
    def min_film_angle_deg(self) -> float:
        return float(np.degrees(self.theta[np.argmin(self.h_min)]))


def journal_orbit(theta, film_force, omega_j, omega_b, b: Bearing, omega: float,
                  film: ShortBearingFilm, cycles: int, steps_per_cycle: int) -> Orbit:
    """Integrate the journal-centre orbit for a periodic (720 deg) load.

    film_force: (2, N) force the film must exert on the journal (engine frame) on grid `theta`
    covering [0, 4 pi). omega_j / omega_b: absolute angular velocities (arrays or scalars).
    """
    th = np.asarray(theta, dtype=float)
    th_c = np.append(th, params.FULL_CYCLE)
    Fz = np.append(film_force[0], film_force[0][0])
    Fy = np.append(film_force[1], film_force[1][0])
    wj = np.append(np.broadcast_to(omega_j, th.shape), np.broadcast_to(omega_j, th.shape)[0])
    wb = np.append(np.broadcast_to(omega_b, th.shape), np.broadcast_to(omega_b, th.shape)[0])
    eps_floor = 1e-12

    def rhs(t, y):
        tt = np.mod(t, params.FULL_CYCLE)
        F = np.array([np.interp(tt, th_c, Fz), np.interp(tt, th_c, Fy)])
        W = -F
        wbar = 0.5 * (np.interp(tt, th_c, wj) + np.interp(tt, th_c, wb))
        eps = math.hypot(y[0], y[1])
        if eps > eps_floor:
            e_r = y / eps
        else:
            Wn = math.hypot(*W)
            e_r = W / Wn if Wn > 0 else np.array([1.0, 0.0])
        e_t = np.array([-e_r[1], e_r[0]])
        W_rt = (W @ e_r, W @ e_t)
        if math.hypot(*W_rt) == 0.0:
            u, v = 0.0, eps * wbar
        else:
            a, A_over = film.solve(W_rt, eps)
            A = A_over / b.force_scale
            u, v = A * math.sin(a), A * math.cos(a)
        edot = u * e_r + (eps * wbar - v) * e_t
        return edot / omega

    y0 = np.zeros(2)
    span = (0.0, cycles * params.FULL_CYCLE)
    t_eval = np.linspace((cycles - 1) * params.FULL_CYCLE, cycles * params.FULL_CYCLE,
                         steps_per_cycle, endpoint=False)
    max_step = params.FULL_CYCLE / steps_per_cycle * 2
    sol = solve_ivp(rhs, span, y0, method="RK45", t_eval=t_eval, rtol=1e-6, atol=1e-8,
                    max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    eps_mag = np.hypot(sol.y[0], sol.y[1])
    return Orbit(t_eval - (cycles - 1) * params.FULL_CYCLE, sol.y, b.radial_clearance * (1 - eps_mag))


def film_for(dspec=None) -> ShortBearingFilm:
    bs = (dspec or params.dyn_spec())["bearings"]
    return ShortBearingFilm(bs["psi_grid_points"], bs["eps_max"])
