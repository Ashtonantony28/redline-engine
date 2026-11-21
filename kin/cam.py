"""Cam lift profile design for direct-acting bucket tappets (milestone M3).

All numbers are modelled (CLAUDE.md rule 5). SI units, angles in radians inside the code
(ADR 0001); degrees only where the spec says ``_deg``.

Profile shape (symmetric about the cam nose)
--------------------------------------------
Let ``a = |phi|`` be the cam angle measured from the nose (event centre). Going outwards:

* **main event**, ``0 <= a <= beta_m``: an even "polydyne-style" polynomial in
  ``u = a / beta_m``::

      y = L + sum_k c_k u^k,   k in main_powers (all even, first one 2)

  Even powers make every odd derivative vanish at the nose, so the profile is C-infinity
  there. The ``c_k`` are fixed by four boundary conditions at the ramp junction ``u = 1``:
  lift = ramp height, velocity = ramp velocity, acceleration = 0, jerk = 0.
* **opening / closing ramp**, ``beta_m < a <= beta_m + beta_r``, measured from the base
  circle by ``t = beta_m + beta_r - a``: first a 7th-order polynomial in ``x = t / beta_a``
  with ``y = y' = y'' = y''' = 0`` on the base circle (x = 0) and ``y = h_a, |y'| = v_r,
  y'' = y''' = 0`` at x = 1, then a constant-velocity section at ``v_r`` up to the ramp
  height ``h_r``. The polynomial part has monotone velocity and acceleration only when its
  mean velocity ``h_a / beta_a`` is about half of ``v_r`` (0.42..0.58), so
  ``h_a = ramp_mean_velocity_ratio * v_r * beta_a``. Lash is taken up on the constant-velocity
  section, so the cam meets the tappet at the known speed ``v_r * omega``.

So lift is C3 everywhere (jerk continuous, hence bounded). Valve lift is cam lift minus
lash (direct-acting bucket, rocker ratio 1); the lash is taken up on the low-velocity ramp.

The main-event half-width ``beta_m`` is solved so that the valve lift equals the timing
check lift (``spec/valvetrain.yaml: timing.check_lift``) exactly at the spec opening and
closing angles (``spec/engine.yaml: ivo/ivc/evo/evc``).

Flat-follower geometry: contact offset from the tappet axis ``e = dy/dphi``; cam radius of
curvature ``rho = r_b + y + d2y/dphi2`` (must stay positive: no undercut).
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq

TWO_PI = 2.0 * math.pi
CAM_PER_CRANK = 0.5          # four-stroke: camshaft turns once per two crank turns
N_DERIV = 4                  # lift, velocity, acceleration, jerk


# ----------------------------------------------------------------------------- polynomials
def solve_power_poly(powers, bcs):
    """Coefficients c_k of p(x) = sum c_k x^k matching derivative conditions at points.

    ``bcs`` is a list of ``(x, derivative_order, value)``; there must be as many conditions
    as powers. Returns an array aligned with ``powers``.
    """
    powers = np.asarray(powers, dtype=float)
    a = np.zeros((len(bcs), len(powers)))
    b = np.zeros(len(bcs))
    for row, (x, order, value) in enumerate(bcs):
        fall = np.ones_like(powers)
        for j in range(order):
            fall = fall * (powers - j)
        expo = powers - order
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(fall == 0.0, 0.0, fall * np.power(float(x), np.maximum(expo, 0.0)))
        a[row] = term
        b[row] = value
    return np.linalg.solve(a, b)


def eval_power_poly(coeffs, powers, x, order: int = 0):
    """d^order/dx^order of sum c_k x^k (terms that differentiate to zero are dropped)."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    for c, k in zip(coeffs, powers):
        fall = 1.0
        for j in range(order):
            fall *= k - j
        if fall == 0.0:
            continue
        out = out + c * fall * x ** (k - order)
    return out


def poly345_coefficients():
    """3-4-5 rise polynomial s(x) = 10x^3 - 15x^4 + 6x^5 solved from its six BCs."""
    powers = (0, 1, 2, 3, 4, 5)
    bcs = [(0.0, 0, 0.0), (0.0, 1, 0.0), (0.0, 2, 0.0), (1.0, 0, 1.0), (1.0, 1, 0.0), (1.0, 2, 0.0)]
    return powers, solve_power_poly(powers, bcs)


# ----------------------------------------------------------------------------- harmonic cam
@dataclass(frozen=True)
class HarmonicCam:
    """Symmetric simple-harmonic cam: y = (h/2)(1 + cos(pi a / beta)) for |phi| = a <= beta.

    Closed-form reference for the analytic tests (lift, derivatives, curvature, contour).
    """

    lift: float
    half_width: float
    base_radius: float

    def derivatives(self, phi):
        phi = np.asarray(phi, dtype=float)
        k = math.pi / self.half_width
        inside = np.abs(phi) <= self.half_width
        half = self.lift / 2.0
        y = np.where(inside, half * (1.0 + np.cos(k * phi)), 0.0)
        yp = np.where(inside, -half * k * np.sin(k * phi), 0.0)
        ypp = np.where(inside, -half * k**2 * np.cos(k * phi), 0.0)
        yppp = np.where(inside, half * k**3 * np.sin(k * phi), 0.0)
        return np.stack([y, yp, ypp, yppp])


# ----------------------------------------------------------------------------- bucket cam
@dataclass(frozen=True)
class CamProfile:
    """Designed symmetric cam lobe. ``phi`` is cam angle from the nose [rad]."""

    name: str
    cam_lift: float            # peak cam (tappet) lift = valve max lift + lash [m]
    lash: float                # [m]
    ramp_height: float         # [m]
    ramp_velocity: float       # [m/rad cam]
    ramp_angle: float          # beta_r [rad cam], whole ramp
    ramp_accel_angle: float    # beta_a [rad cam], polynomial part of the ramp
    ramp_accel_height: float   # h_a [m]
    main_half_angle: float     # beta_m [rad cam]
    main_powers: tuple
    main_coeffs: tuple
    ramp_powers: tuple
    ramp_coeffs: tuple
    base_radius: float         # [m]
    centre_crank: float        # crank angle of the nose, ADR 0001 origin [rad]

    @property
    def half_event(self) -> float:
        return self.main_half_angle + self.ramp_angle

    @property
    def valve_max_lift(self) -> float:
        return self.cam_lift - self.lash

    def derivatives(self, phi):
        """Stack [y, dy/dphi, d2y/dphi2, d3y/dphi3] of cam lift at cam angles ``phi``."""
        phi_in = np.asarray(phi, dtype=float)
        phi = np.atleast_1d(phi_in)
        a = np.abs(phi)
        sgn = np.where(phi < 0.0, -1.0, 1.0)
        out = np.zeros((N_DERIV,) + phi.shape)
        bm, br = self.main_half_angle, self.ramp_angle

        main = a <= bm
        if np.any(main):
            u = a[main] / bm
            for n in range(N_DERIV):
                out[n][main] = eval_power_poly(self.main_coeffs, self.main_powers, u, n) / bm**n

        ba = self.ramp_accel_angle
        t = bm + br - a                                  # angle from the base circle
        poly = (t > 0.0) & (t <= ba)
        if np.any(poly):
            x = t[poly] / ba
            for n in range(N_DERIV):
                out[n][poly] = (-1.0) ** n * eval_power_poly(self.ramp_coeffs, self.ramp_powers, x, n) / ba**n
        const = (t > ba) & (a > bm)
        if np.any(const):
            out[0][const] = self.ramp_accel_height + self.ramp_velocity * (t[const] - ba)
            out[1][const] = -self.ramp_velocity
        # f was written in a = |phi|; odd derivatives flip sign on the opening (phi < 0) side
        out[1] *= sgn
        out[3] *= sgn
        return out.reshape((N_DERIV,) + phi_in.shape)

    def lift(self, phi):
        return self.derivatives(phi)[0]

    def valve_lift(self, phi):
        return np.maximum(self.lift(phi) - self.lash, 0.0)

    def radius_of_curvature(self, phi):
        d = self.derivatives(phi)
        return flat_follower_radius(self.base_radius, d[0], d[2])

    def contact_offset(self, phi):
        return self.derivatives(phi)[1]

    # crank-angle helpers ------------------------------------------------------------
    def cam_angle_from_crank(self, theta, advance: float = 0.0):
        """Cam angle from the nose for crank angle ``theta`` [rad], cam advanced by
        ``advance`` crank radians (events happen earlier)."""
        rel = np.asarray(theta, dtype=float) - (self.centre_crank - advance)
        rel = (rel + TWO_PI) % (2.0 * TWO_PI) - TWO_PI      # wrap into one 720 deg cycle
        return CAM_PER_CRANK * rel

    def valve_lift_crank(self, theta, advance: float = 0.0):
        return self.valve_lift(self.cam_angle_from_crank(theta, advance))

    def valve_open_crank(self):
        """Crank angles [rad] where the valve leaves / returns to its seat (lash point)."""
        a = brentq(lambda t: self.lift(np.array(t)) - self.lash, 0.0, self.half_event)
        return self.centre_crank - a / CAM_PER_CRANK, self.centre_crank + a / CAM_PER_CRANK


def flat_follower_radius(base_radius, y, ypp):
    """Radius of curvature of a cam driving a flat-faced translating follower."""
    return base_radius + y + ypp


def flat_follower_contour(base_radius, phi, y, yp):
    """Cam surface points in the cam-fixed frame for a flat translating follower.

    The follower face sits at distance ``r_b + y`` from the cam centre; the contact point is
    offset ``y'`` along the face. Rotating back by the cam angle gives the cam contour.
    """
    rr = base_radius + y
    c, s = np.cos(phi), np.sin(phi)
    x = rr * c - yp * s
    z = rr * s + yp * c
    return x, z


def contour_curvature_radius(x, z, phi):
    """Radius of curvature of a parametric planar curve (x(phi), z(phi)), by finite differences."""
    dx, dz = np.gradient(x, phi), np.gradient(z, phi)
    ddx, ddz = np.gradient(dx, phi), np.gradient(dz, phi)
    return (dx**2 + dz**2) ** 1.5 / np.abs(dx * ddz - dz * ddx)


# ----------------------------------------------------------------------------- design
RAMP_POWERS = (4, 5, 6, 7)   # y, y', y'', y''' = 0 at the base circle -> only x^4..x^7


def _ramp_coeffs(height, velocity, angle):
    # in x, dy/dx = (dy/da)(da/dx) = (-v)(-beta_r): junction velocity v * beta_r
    bcs = [(1.0, 0, height), (1.0, 1, velocity * angle), (1.0, 2, 0.0), (1.0, 3, 0.0)]
    return solve_power_poly(RAMP_POWERS, bcs)


def _main_coeffs(cam_lift, height, velocity, main_half, powers):
    # f(u) = L + sum c u^k; at u = 1: f = h_r, df/da = -v_r, f'' = f''' = 0
    shifted = solve_power_poly(
        powers,
        [(1.0, 0, height - cam_lift), (1.0, 1, -velocity * main_half), (1.0, 2, 0.0), (1.0, 3, 0.0)],
    )
    return shifted


def build_profile(name, cam_lift, lash, ramp_height, ramp_velocity, ramp_accel_angle,
                  ramp_mean_ratio, main_half, powers, base_radius, centre_crank) -> CamProfile:
    powers = tuple(int(p) for p in powers)
    main = _main_coeffs(cam_lift, ramp_height, ramp_velocity, main_half, powers)
    accel_height = ramp_mean_ratio * ramp_velocity * ramp_accel_angle
    if accel_height >= ramp_height:
        raise ValueError("ramp polynomial section would exceed the ramp height")
    ramp_angle = ramp_accel_angle + (ramp_height - accel_height) / ramp_velocity
    ramp = _ramp_coeffs(accel_height, ramp_velocity, ramp_accel_angle)
    # the constant L is carried as an extra power-0 term so eval_power_poly handles it
    return CamProfile(
        name=name, cam_lift=cam_lift, lash=lash, ramp_height=ramp_height,
        ramp_velocity=ramp_velocity, ramp_angle=ramp_angle,
        ramp_accel_angle=ramp_accel_angle, ramp_accel_height=accel_height, main_half_angle=main_half,
        main_powers=(0,) + powers, main_coeffs=(cam_lift,) + tuple(main),
        ramp_powers=RAMP_POWERS, ramp_coeffs=tuple(ramp), base_radius=base_radius,
        centre_crank=centre_crank,
    )


def design_cam(name: str, engine_spec: dict, vt_spec: dict) -> CamProfile:
    """Design the intake or exhaust lobe from spec timing and lift (see module docstring)."""
    v = engine_spec["valves"]
    cam = vt_spec["cam"]
    lobe = cam[name]
    if name == "intake":
        open_deg, close_deg, max_lift = v["ivo_deg"], v["ivc_deg"], v["intake_max_lift"]
    elif name == "exhaust":
        open_deg, close_deg, max_lift = v["evo_deg"], v["evc_deg"], v["exhaust_max_lift"]
    else:
        raise ValueError(name)
    lash = lobe["lash"]
    cam_lift = max_lift + lash
    accel_angle = math.radians(lobe["ramp_accel_angle_deg"])
    centre = math.radians(0.5 * (open_deg + close_deg))
    check_half = CAM_PER_CRANK * math.radians(0.5 * (close_deg - open_deg))
    check_lift = vt_spec["timing"]["check_lift"]

    def make(bm):
        return build_profile(name, cam_lift, lash, lobe["ramp_height"], lobe["ramp_velocity"],
                             accel_angle, cam["ramp_mean_velocity_ratio"], bm, cam["main_powers"], cam["base_circle_radius"], centre)

    def miss(bm):
        return float(make(bm).valve_lift(np.array(check_half))) - check_lift

    # bracket: at lo the check angle is the end of the event (lift 0 < check); at hi it is
    # well inside the main event (lift > check)
    lo = check_half - make(check_half).ramp_angle
    hi = check_half * cam["main_half_bracket_ratio"]
    bm = brentq(miss, lo, hi, xtol=1e-14)
    return make(bm)


def design_cams(engine_spec: dict, vt_spec: dict):
    return {n: design_cam(n, engine_spec, vt_spec) for n in ("intake", "exhaust")}
