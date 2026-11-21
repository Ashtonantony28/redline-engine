"""Crankshaft torsional vibration: lumped chain model, Holzer and eigen solutions, order analysis,
forced response and a tuned rubber (Den Hartog) damper.

Chain (front to rear): [damper ring] - pulley/hub - throw 1..4 - flywheel. Stiffness of each
segment k = G J_p / L_eq for a plain shaft of main-journal diameter and the equivalent length in
spec/dynamics.yaml. Throw inertia = crank share + rod big end m_be r^2 + mean reciprocating
equivalent m_rec r^2 / 2 (the variable-inertia effect of the slider-crank is neglected).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from dyn import params
from dyn.forces import cylinder_forces
from dyn.pressure import PressureTrace


@dataclass
class TorsionalModel:
    J: np.ndarray               # station inertias [kg m^2]
    k: np.ndarray               # segment stiffnesses [N m/rad], len(J) - 1
    c_abs: np.ndarray           # absolute (to ground) damping per station [N m s/rad]
    c_rel: np.ndarray           # relative (segment) damping [N m s/rad]
    labels: list = field(default_factory=list)
    throws: list = field(default_factory=list)    # station indices of throws 1..n (cylinder order)
    hub: int = 0                                  # station index of the pulley hub

    def K(self) -> np.ndarray:
        return chain_matrix(self.k)

    def C(self) -> np.ndarray:
        return chain_matrix(self.c_rel) + np.diag(self.c_abs)


def chain_matrix(k) -> np.ndarray:
    k = np.asarray(k, dtype=float)
    n = k.size + 1
    M = np.zeros((n, n))
    for i, ki in enumerate(k):
        M[i, i] += ki
        M[i + 1, i + 1] += ki
        M[i, i + 1] -= ki
        M[i + 1, i] -= ki
    return M


def two_mass_frequency(J1: float, J2: float, k: float) -> float:
    """Closed-form natural frequency of a free-free two-mass torsional system [Hz]."""
    return math.sqrt(k * (J1 + J2) / (J1 * J2)) / (2 * math.pi)


def build_model(spec=None, dspec=None) -> TorsionalModel:
    """Bare crank chain (no damper ring); damping set from spec/dynamics.yaml."""
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    ts = dspec["torsional"]
    g = params.geometry(spec)
    m = params.masses(spec)
    n = spec["engine"]["cylinders"]
    ci = spec["masses"]["crank_inertia"]
    r2 = g.crank_radius**2
    j_throw = ci * (1 - ts["crank_nose_fraction"] - ts["crank_flange_fraction"]) / n \
        + m.rod_big_end * r2 + 0.5 * m.reciprocating * r2
    J = np.array([ts["pulley_hub_inertia"] + ci * ts["crank_nose_fraction"]]
                 + [j_throw] * n
                 + [spec["masses"]["flywheel_inertia"] + ci * ts["crank_flange_fraction"]])
    G = params.shear_modulus(spec, dspec)
    d = dspec["bearings"]["main"]["diameter"]
    GJp = G * math.pi * d**4 / 32
    L = ts["equivalent_length"]
    lengths = [L["nose_to_throw1"]] + [L["throw_to_throw"]] * (n - 1) + [L["throw4_to_flywheel"]]
    k = GJp / np.array(lengths)
    c_abs = np.zeros(J.size)
    c_abs[1:1 + n] = ts["absolute_damping_per_throw"]
    model = TorsionalModel(J, k, c_abs, np.zeros(k.size),
                           ["hub"] + [f"throw{i + 1}" for i in range(n)] + ["flywheel"],
                           list(range(1, 1 + n)), 0)
    # stiffness-proportional segment damping giving the stated modal damping in mode 1
    w1 = 2 * math.pi * eigen_frequencies(model.J, model.k)[1]
    model.c_rel = 2 * ts["structural_damping_ratio"] / w1 * model.k
    return model


# --- natural frequencies -------------------------------------------------------------------

def holzer_residual(omega: float, J, k) -> float:
    """Holzer table for a free-free chain: residual torque past the last station for unit
    amplitude at station 1. Zero at a natural frequency."""
    theta, torque = 1.0, 0.0
    for i, Ji in enumerate(J):
        torque += Ji * omega**2 * theta
        if i < len(k):
            theta -= torque / k[i]
    return torque


def holzer_mode(omega: float, J, k) -> np.ndarray:
    theta, torque, shape = 1.0, 0.0, []
    for i, Ji in enumerate(J):
        shape.append(theta)
        torque += Ji * omega**2 * theta
        if i < len(k):
            theta -= torque / k[i]
    return np.array(shape)


def holzer_frequencies(J, k, f_max: float, n_scan: int = 20000) -> np.ndarray:
    """Natural frequencies [Hz] below f_max by scanning the Holzer residual and bracketing roots
    (the rigid-body zero at 0 Hz is included as the first entry)."""
    w = np.linspace(0.0, 2 * math.pi * f_max, n_scan)[1:]
    res = np.array([holzer_residual(wi, J, k) for wi in w])
    roots = [0.0]
    for i in np.nonzero(np.sign(res[:-1]) != np.sign(res[1:]))[0]:
        roots.append(brentq(holzer_residual, w[i], w[i + 1], args=(J, k), xtol=1e-12))
    return np.array(roots) / (2 * math.pi)


def eigen_frequencies(J, k, return_modes: bool = False):
    """Natural frequencies [Hz] from the generalized eigenproblem K phi = w^2 J phi."""
    lam, vec = eigh(chain_matrix(k), np.diag(J))
    f = np.sqrt(np.clip(lam, 0.0, None)) / (2 * math.pi)
    return (f, vec) if return_modes else f


# --- excitation ---------------------------------------------------------------------------

def torque_orders(trace: PressureTrace, max_order: float, n: int = 2880, spec=None, dspec=None):
    """Complex order coefficients of single-cylinder gas torque and of the inertia torque per
    unit omega^2: T(t) = sum_q Re(C_q e^{i q t}), q = 0.5, 1, ..., max_order.
    Returns (orders, C_gas, C_inertia_per_w2)."""
    th = np.linspace(0.0, params.FULL_CYCLE, n, endpoint=False)
    f = cylinder_forces(th, trace, 1.0, spec, dspec)
    orders = np.arange(1, int(round(2 * max_order)) + 1) / 2.0
    E = np.exp(-1j * orders[:, None] * th[None, :])
    return orders, 2 * (E * f.T_gas).mean(1), 2 * (E * f.T_inertia).mean(1)


def phase_sum(mode: np.ndarray, throws, order: float, spec=None) -> float:
    """|sum_k phi_k e^{-i q phi_fire,k}| over throws: how strongly order q excites the mode."""
    ph = params.firing_offsets(spec)
    return float(abs(np.sum(mode[throws] * np.exp(-1j * order * ph))))


# --- forced response ----------------------------------------------------------------------

def forced_response(model: TorsionalModel, rpm, orders, C_gas, C_in, spec=None):
    """Steady-state complex angle amplitudes per order: array (n_orders, n_rpm, n_stations)."""
    ph = params.firing_offsets(spec)
    K, C, J = model.K(), model.C(), np.diag(model.J)
    rpm = np.asarray(rpm, dtype=float)
    out = np.zeros((len(orders), rpm.size, model.J.size), dtype=complex)
    for j, n_rpm in enumerate(rpm):
        w = params.omega_from_rpm(n_rpm)
        for i, q in enumerate(orders):
            Om = q * w
            Tq = C_gas[i] + w**2 * C_in[i]
            F = np.zeros(model.J.size, dtype=complex)
            F[model.throws] = Tq * np.exp(-1j * q * ph)
            out[i, j] = np.linalg.solve(K - Om**2 * J + 1j * Om * C, F)
    return out


# --- tuned damper -------------------------------------------------------------------------

@dataclass
class DamperDesign:
    ring_inertia: float
    stiffness: float
    damping: float
    mass_ratio: float
    tuning_ratio: float
    damping_ratio: float
    modal_inertia_at_hub: float
    target_hz: float


def den_hartog_damper(model: TorsionalModel, ring_inertia: float, mode: int = 1) -> DamperDesign:
    """Den Hartog optimum for a damped vibration absorber attached at the hub, applied to the
    equivalent single-DOF system of `mode` (modal inertia referred to the hub station):
    f = 1/(1+mu), zeta = sqrt(3 mu / (8 (1+mu)^3))."""
    f, vec = eigen_frequencies(model.J, model.k, return_modes=True)
    phi = vec[:, mode]
    J_eq = float(phi @ (model.J * phi)) / phi[model.hub] ** 2
    mu = ring_inertia / J_eq
    wn = 2 * math.pi * f[mode]
    tune = 1.0 / (1.0 + mu)
    zeta = math.sqrt(3 * mu / (8 * (1 + mu) ** 3))
    wd = tune * wn
    return DamperDesign(ring_inertia, ring_inertia * wd**2, 2 * zeta * ring_inertia * wd,
                        mu, tune, zeta, J_eq, f[mode])


def with_damper(model: TorsionalModel, d: DamperDesign) -> TorsionalModel:
    """Prepend the damper ring, coupled to the hub through the rubber (k_d, c_d)."""
    return replace(
        model,
        J=np.concatenate([[d.ring_inertia], model.J]),
        k=np.concatenate([[d.stiffness], model.k]),
        c_abs=np.concatenate([[0.0], model.c_abs]),
        c_rel=np.concatenate([[d.damping], model.c_rel]),
        labels=["damper_ring"] + list(model.labels),
        throws=[t + 1 for t in model.throws],
        hub=model.hub + 1,
    )
