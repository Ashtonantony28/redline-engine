"""Valvetrain dynamics: single-mass quasi-static float estimate and a multi-mass model.

All results are modelled (CLAUDE.md rule 5). SI units (ADR 0001).

Multi-mass model (positions measured in the valve-opening direction, which compresses the
spring; origin = valve closed with zero seat load)::

    cam --k_c (unilateral, lash)--> tappet --k_stem (unilateral)--> valve --spring chain--> head
                                                                     |
                                                                    seat (unilateral)

The spring of rate ``k`` and mass ``m_s`` is ``n`` internal lumped masses ``m_s/(n+1)`` joined
by ``n+1`` segments of stiffness ``(n+1) k``; the end half-segment masses go to the head
(fixed) and to the valve. Damping in the chain is stiffness proportional, set so the first
fixed-fixed surge mode has the spec damping ratio. With ``n = 4`` the model has 6 masses.

Integration is semi-implicit (symplectic) Euler at fixed cam-angle steps, batched over
engine speeds, so all speeds share one precomputed cam-lift table.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import eigh

from kin.cam import CAM_PER_CRANK, TWO_PI, CamProfile

SECONDS_PER_MINUTE = 60.0


def cam_omega(engine_rpm):
    """Camshaft angular speed [rad/s] for engine speed [rpm]."""
    return np.asarray(engine_rpm, dtype=float) * TWO_PI / SECONDS_PER_MINUTE * CAM_PER_CRANK


def engine_rpm_from_cam_omega(omega):
    return np.asarray(omega, dtype=float) / CAM_PER_CRANK * SECONDS_PER_MINUTE / TWO_PI


def spring_force(x, preload, rate):
    """Static (Hooke) spring force at valve lift x."""
    return preload + rate * np.asarray(x, dtype=float)


def surge_frequencies_continuous(rate, spring_mass, modes: int = 3):
    """Fixed-fixed helical spring surge: f_j = (j/2) sqrt(k/m_s) [Hz]."""
    j = np.arange(1, modes + 1)
    return 0.5 * j * math.sqrt(rate / spring_mass)


def surge_frequencies_discrete(rate, spring_mass, n, modes: int = 3):
    """Exact natural frequencies [Hz] of the n-mass fixed-fixed chain used here."""
    ke, m = (n + 1) * rate, spring_mass / (n + 1)
    j = np.arange(1, min(modes, n) + 1)
    return 2.0 * math.sqrt(ke / m) * np.sin(j * math.pi / (2 * (n + 1))) / TWO_PI


@dataclass(frozen=True)
class ValvetrainModel:
    name: str
    m_tappet: float
    m_valve: float           # valve + retainer + keepers (spring end mass added in .masses)
    spring_rate: float
    preload: float
    spring_mass: float
    n_spring: int
    k_contact: float
    k_stem: float
    k_seat: float
    zeta_contact: float
    zeta_stem: float
    zeta_seat: float
    zeta_guide: float
    zeta_surge: float

    # --- lumped parameters ------------------------------------------------------------
    @property
    def segment_stiffness(self) -> float:
        return (self.n_spring + 1) * self.spring_rate

    @property
    def node_mass(self) -> float:
        return self.spring_mass / (self.n_spring + 1)

    @property
    def masses(self) -> np.ndarray:
        """[tappet, valve (+ half end segment), spring nodes 1..n (valve side first)]."""
        valve = self.m_valve + 0.5 * self.node_mass
        return np.array([self.m_tappet, valve] + [self.node_mass] * self.n_spring)

    @property
    def single_mass(self) -> float:
        """Equivalent single-DOF mass: tappet + valve + 1/3 of the spring."""
        return self.m_tappet + self.m_valve + self.spring_mass / 3.0

    @property
    def segment_damping(self) -> float:
        omega1 = TWO_PI * surge_frequencies_discrete(self.spring_rate, self.spring_mass, self.n_spring, 1)[0]
        return 2.0 * self.zeta_surge * self.segment_stiffness / omega1

    def _local_damping(self, zeta, k, m):
        return 2.0 * zeta * math.sqrt(k * m)

    @property
    def c_contact(self):
        return self._local_damping(self.zeta_contact, self.k_contact, self.m_tappet)

    @property
    def c_stem(self):
        return self._local_damping(self.zeta_stem, self.k_stem, self.m_tappet)

    @property
    def c_seat(self):
        return self._local_damping(self.zeta_seat, self.k_seat, self.m_valve)

    @property
    def c_guide(self):
        return self._local_damping(self.zeta_guide, self.spring_rate, self.single_mass)

    # --- linear analysis --------------------------------------------------------------
    def stiffness_matrix(self, contact=True, seat=False) -> np.ndarray:
        """Linearised K with chosen unilateral contacts closed."""
        n = 2 + self.n_spring
        k = np.zeros((n, n))

        def link(i, j, kk):
            k[i, i] += kk
            if j is not None:
                k[j, j] += kk
                k[i, j] -= kk
                k[j, i] -= kk

        if contact:
            link(0, None, self.k_contact)
        link(0, 1, self.k_stem)
        if seat:
            link(1, None, self.k_seat)
        # chain: valve(1) - node 2 - ... - node n+1 - head
        ke = self.segment_stiffness
        chain = [1] + list(range(2, n))
        for a, b in zip(chain[:-1], chain[1:]):
            link(a, b, ke)
        link(chain[-1], None, ke)
        return k

    def natural_frequencies(self, contact=True, seat=False) -> np.ndarray:
        w2 = eigh(self.stiffness_matrix(contact, seat), np.diag(self.masses), eigvals_only=True)
        return np.sqrt(np.maximum(w2, 0.0)) / TWO_PI

    def spring_chain_frequencies(self) -> np.ndarray:
        """Surge modes of the spring alone, both ends fixed (valve held), via eigen solve."""
        n = self.n_spring
        ke = self.segment_stiffness
        k = 2.0 * ke * np.eye(n) - ke * (np.eye(n, k=1) + np.eye(n, k=-1))
        w2 = eigh(k, self.node_mass * np.eye(n), eigvals_only=True)
        return np.sqrt(w2) / TWO_PI

    def static_valve_force(self, x_valve: float) -> float:
        """Spring force on the valve with the chain in static equilibrium at valve lift x."""
        n = self.n_spring
        ke = self.segment_stiffness
        # internal nodes: K_ii x_i = coupling to the valve
        kii = 2.0 * ke * np.eye(n) - ke * (np.eye(n, k=1) + np.eye(n, k=-1))
        rhs = np.zeros(n)
        rhs[0] = ke * x_valve
        xi = np.linalg.solve(kii, rhs)
        return self.preload + ke * (x_valve - xi[0])


def model_from_spec(name: str, vt_spec: dict) -> ValvetrainModel:
    va, sp = vt_spec["valve_assembly"], vt_spec["spring"]
    st, dz = vt_spec["stiffness"], vt_spec["damping_ratio"]
    return ValvetrainModel(
        name=name,
        m_tappet=vt_spec["tappet"]["mass"],
        m_valve=va[f"{name}_valve_mass"] + va["retainer_keepers_mass"],
        spring_rate=sp["rate"], preload=sp["installed_preload"], spring_mass=sp["mass"],
        n_spring=int(sp["n_masses"]),
        k_contact=st["cam_contact"], k_stem=st["stem"], k_seat=st["seat"],
        zeta_contact=dz["cam_contact"], zeta_stem=dz["stem"], zeta_seat=dz["seat"],
        zeta_guide=dz["guide"], zeta_surge=sp["surge_damping_ratio"],
    )


# ----------------------------------------------------------------------------- single mass
def quasi_static_contact_force(model: ValvetrainModel, profile: CamProfile, phi, engine_rpm):
    """Single-mass rigid valvetrain: F_c = spring(valve lift) + m_eq * y'' * omega^2."""
    d = profile.derivatives(phi)
    xv = np.maximum(d[0] - profile.lash, 0.0)
    w = cam_omega(engine_rpm)
    f = spring_force(xv, model.preload, model.spring_rate) + model.single_mass * d[2] * w**2
    return np.where(d[0] > profile.lash, f, 0.0)


def quasi_static_float_rpm(model: ValvetrainModel, profile: CamProfile, samples: int):
    """Speed where the rigid single-mass contact force first reaches zero."""
    phi = np.linspace(-profile.main_half_angle, profile.main_half_angle, samples)
    d = profile.derivatives(phi)
    xv = d[0] - profile.lash
    neg = d[2] < 0.0
    w2 = spring_force(xv[neg], model.preload, model.spring_rate) / (model.single_mass * -d[2][neg])
    return float(engine_rpm_from_cam_omega(math.sqrt(w2.min())))


def free_vibration_period(mass, stiffness, x0, steps_per_period: int, periods: int):
    """Period of an undamped single-DOF oscillator integrated with the same semi-implicit
    Euler scheme as the multi-mass model (zero up-crossings of the velocity)."""
    t_ref = TWO_PI * math.sqrt(mass / stiffness)
    dt = t_ref / steps_per_period
    x, v = x0, 0.0
    crossings = []
    for i in range(steps_per_period * periods):
        v_old = v
        v += -stiffness / mass * x * dt
        x += v * dt
        if v_old < 0.0 <= v:
            crossings.append((i + 1 + v_old / (v_old - v)) * dt)
    return float(np.mean(np.diff(crossings)))


# ----------------------------------------------------------------------------- multi mass
@dataclass
class SimResult:
    rpm: np.ndarray
    min_contact_main: np.ndarray     # min cam contact force in the main event, last rev [N]
    max_contact: np.ndarray          # [N]
    max_stress_ratio: np.ndarray     # max F_c / rho over the last rev [N/m]
    max_seat_speed: np.ndarray       # max closing speed at valve seating [m/s]
    phi: np.ndarray                  # recorded cam angles (last rev, decimated)
    contact: np.ndarray              # (len(phi), len(rpm)) contact force [N]
    valve: np.ndarray                # valve lift [m]

    @property
    def floats(self) -> np.ndarray:
        return self.min_contact_main <= 0.0


def simulate(model: ValvetrainModel, profile: CamProfile, engine_rpm, revolutions: int,
             step_safety: float, record_points: int = 721) -> SimResult:
    rpm = np.atleast_1d(np.asarray(engine_rpm, dtype=float))
    omega = cam_omega(rpm)
    w_max = TWO_PI * model.natural_frequencies(contact=True, seat=True).max()
    dt_max = step_safety / w_max
    n_rev = int(math.ceil(TWO_PI / (omega.min() * dt_max)))
    dphi = TWO_PI / n_rev
    dt = dphi / omega                                   # per speed

    phi_tab = -math.pi + dphi * np.arange(n_rev)        # start on the base circle
    d = profile.derivatives(phi_tab)
    y_tab, yp_tab = d[0], d[1]
    rho_tab = profile.radius_of_curvature(phi_tab)
    main_tab = y_tab > profile.ramp_height

    m = model.masses[:, None]
    n = len(model.masses)
    ke, ce = model.segment_stiffness, model.segment_damping
    kc, cc = model.k_contact, model.c_contact
    ks, cs = model.k_stem, model.c_stem
    kz, cz = model.k_seat, model.c_seat
    cg = model.c_guide
    f0 = model.preload
    nspeeds = rpm.size

    # static equilibrium on the seat: chain uniformly preloaded, valve compresses the seat
    x = np.zeros((n, nspeeds))
    x_seat = -f0 / kz
    x[0] = x_seat
    x[1] = x_seat
    for i in range(model.n_spring):                     # node i+2 is (n_spring - i) segments from head
        x[2 + i] = x_seat * (model.n_spring - i) / (model.n_spring + 1)
    v = np.zeros_like(x)

    min_main = np.full(nspeeds, np.inf)
    max_fc = np.zeros(nspeeds)
    max_ratio = np.zeros(nspeeds)
    max_seat = np.zeros(nspeeds)
    rec_every = max(n_rev // (record_points - 1), 1)
    rec_phi, rec_fc, rec_xv = [], [], []
    force = np.zeros_like(x)
    seg_x = np.zeros((model.n_spring + 1, nspeeds))
    seg_v = np.zeros_like(seg_x)
    last = (revolutions - 1) * n_rev

    for step in range(revolutions * n_rev):
        k = step % n_rev
        y, ydot = y_tab[k], yp_tab[k] * omega
        # cam / tappet
        dc = y - profile.lash - x[0]
        fc = np.where(dc > 0.0, np.maximum(kc * dc + cc * (ydot - v[0]), 0.0), 0.0)
        # tappet / valve stem
        ds = x[0] - x[1]
        fs = np.where(ds > 0.0, np.maximum(ks * ds + cs * (v[0] - v[1]), 0.0), 0.0)
        # seat
        fz = np.where(x[1] < 0.0, np.maximum(-kz * x[1] - cz * v[1], 0.0), 0.0)
        # spring chain: segment j between node (valve, n1, ..., nN) and the next toward the head
        seg_x[:-1] = x[1:-1] - x[2:]
        seg_x[-1] = x[-1]
        seg_v[:-1] = v[1:-1] - v[2:]
        seg_v[-1] = v[-1]
        fseg = ke * seg_x + ce * seg_v                  # compression increments (preload separate)
        force[0] = fc - fs
        force[1] = fs + fz - (f0 + fseg[0]) - cg * v[1]
        force[2:] = fseg[:-1] - fseg[1:]
        v += force / m * dt
        x += v * dt

        if step >= last:
            if main_tab[k]:
                min_main = np.minimum(min_main, fc)
            max_fc = np.maximum(max_fc, fc)
            max_ratio = np.maximum(max_ratio, fc / rho_tab[k])
            seating = (fz > 0.0) & (v[1] < 0.0)
            max_seat = np.where(seating, np.maximum(max_seat, -v[1]), max_seat)
            if (step - last) % rec_every == 0:
                rec_phi.append(phi_tab[k])
                rec_fc.append(fc.copy())
                rec_xv.append(np.maximum(x[1], 0.0))

    return SimResult(rpm, min_main, max_fc, max_ratio, max_seat, np.array(rec_phi),
                     np.array(rec_fc), np.array(rec_xv))


def float_speed(model: ValvetrainModel, profile: CamProfile, sweep, refine_step, revolutions,
                step_safety):
    """Lowest engine speed with cam/tappet separation in the main event (multi-mass).

    Coarse sweep then a refine pass; returns (float_rpm or inf, coarse SimResult)."""
    start, stop, step = sweep
    rpm = np.arange(start, stop + step / 2, step, dtype=float)
    coarse = simulate(model, profile, rpm, revolutions, step_safety)
    if not coarse.floats.any():
        return math.inf, coarse
    i = int(np.argmax(coarse.floats))
    if i == 0:
        return float(rpm[0]), coarse
    fine_rpm = np.arange(rpm[i - 1], rpm[i] + refine_step / 2, refine_step)
    fine = simulate(model, profile, fine_rpm, revolutions, step_safety)
    j = int(np.argmax(fine.floats)) if fine.floats.any() else len(fine_rpm) - 1
    return float(fine_rpm[j]), coarse
