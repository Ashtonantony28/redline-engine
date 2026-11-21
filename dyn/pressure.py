"""Cylinder-pressure input interface for the dynamics package, plus a PLACEHOLDER generator.

Interface: `PressureTrace` holds absolute cylinder pressure p(theta) of ONE cylinder over one
four-stroke cycle, theta in [0, 4 pi] rad local to that cylinder (0 = its firing TDC, ADR 0001).
Any source (the thermo/ two-zone model of M2, a measured trace, a fixture) only has to build a
`PressureTrace` (e.g. `PressureTrace.from_degrees(theta_deg, p_pa, source="thermo.run")`).

PLACEHOLDER: `placeholder_wiebe_trace` is a single-zone, constant-gamma, Wiebe heat-release
stand-in with idealised gas exchange. It exists only so M5 can run before M2; it must be
replaced by thermo/ output. It has no heat-transfer model beyond a lumped heat-to-gas fraction,
no blow-by, no residual gas and no real gas properties.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import cumulative_trapezoid

from kin.slider_crank import cylinder_volume
from dyn import params

PLACEHOLDER_SOURCE = "PLACEHOLDER single-zone Wiebe (dyn.pressure); replace with thermo/ output"
DEG = math.pi / 180.0


@dataclass(frozen=True)
class PressureTrace:
    """Absolute cylinder pressure over one 720-deg cycle, local crank angle [rad]."""

    theta: np.ndarray
    p: np.ndarray
    source: str = "unspecified"
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        th = np.asarray(self.theta, dtype=float)
        p = np.asarray(self.p, dtype=float)
        if th.ndim != 1 or th.shape != p.shape:
            raise ValueError("theta and p must be 1-D arrays of equal length")
        if np.any(np.diff(th) <= 0):
            raise ValueError("theta must be strictly increasing")
        tol = 1e-9
        if th[0] > tol or th[-1] < params.FULL_CYCLE - 1e-6:
            raise ValueError("theta must cover 0..4 pi (0..720 deg)")
        if np.any(p <= 0):
            raise ValueError("pressure must be absolute and positive")
        object.__setattr__(self, "theta", th)
        object.__setattr__(self, "p", p)

    @classmethod
    def from_degrees(cls, theta_deg, p, source="unspecified", meta=None):
        return cls(np.asarray(theta_deg, dtype=float) * DEG, p, source, meta or {})

    def at(self, theta_local):
        """p at local crank angle(s), periodic over 720 deg (linear interpolation)."""
        th = np.mod(np.asarray(theta_local, dtype=float), params.FULL_CYCLE)
        return np.interp(th, self.theta, self.p)

    @property
    def peak(self) -> float:
        return float(self.p.max())

    @property
    def peak_angle_deg(self) -> float:
        return float(self.theta[np.argmax(self.p)] / DEG)


def wiebe_fraction(theta, theta0, duration, a, m):
    """Wiebe mass fraction burned x_b(theta) (Heywood eq. 9.23)."""
    y = np.clip((np.asarray(theta, dtype=float) - theta0) / duration, 0.0, None)
    return 1.0 - np.exp(-a * y ** (m + 1))


def wiebe_rate(theta, theta0, duration, a, m):
    """d x_b / d theta [1/rad]."""
    y = np.clip((np.asarray(theta, dtype=float) - theta0) / duration, 0.0, None)
    return a * (m + 1) / duration * y**m * np.exp(-a * y ** (m + 1))


def placeholder_wiebe_trace(op: str = "full_load", spec: dict | None = None,
                            dspec: dict | None = None) -> PressureTrace:
    """PLACEHOLDER pressure trace (see module docstring). `op` names an operating point in
    spec/dynamics.yaml pressure_placeholder.operating_points."""
    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    geom = params.geometry(spec)
    pp = dspec["pressure_placeholder"]
    point = pp["operating_points"][op]
    v = spec["valves"]
    c = spec["combustion"]
    gamma, R = pp["gamma"], pp["gas_constant"]
    p_int, p_exh = point["intake_pressure"], point["exhaust_pressure"]
    cycle_deg = params.FULL_CYCLE / DEG

    ivo, ivc, evo, evc = v["ivo_deg"], v["ivc_deg"], v["evo_deg"], v["evc_deg"]
    bd = pp["blowdown_deg"]
    if not (evo + bd < ivo < evc < ivc):
        raise ValueError("valve events inconsistent with the placeholder gas-exchange model")

    # closed part IVC -> EVO (+ one cycle), integrated on its own fine grid
    th_c = np.linspace(ivc, evo + cycle_deg, int(round((evo + cycle_deg - ivc) * 20)) + 1) * DEG
    V = cylinder_volume(th_c, geom)
    V_ivc = V[0]
    m_air = p_int * V_ivc / (R * pp["intake_temperature"])
    m_fuel = m_air / (c["afr_stoich"] * pp["lambda"])
    q_total = pp["heat_to_gas_fraction"] * m_fuel * float(c["lhv"])  # float(): PyYAML reads "44.3e6" as str
    theta0 = (cycle_deg + point["soc_deg"]) * DEG
    dur = c["burn_duration_deg"] * DEG
    dq = q_total * wiebe_rate(th_c, theta0, dur, c["wiebe_a"], c["wiebe_m"])
    # exact for constant gamma: d(p V^g)/dtheta = (g - 1) V^(g-1) dQ/dtheta
    pvg = p_int * V_ivc**gamma + (gamma - 1) * cumulative_trapezoid(V ** (gamma - 1) * dq, th_c, initial=0.0)
    p_closed = pvg / V**gamma
    p_evo = p_closed[-1]

    n = pp["grid_points"]
    th_deg = np.linspace(0.0, cycle_deg, n)
    p = np.empty(n)
    d_unwrapped = np.where(th_deg < ivc, th_deg + cycle_deg, th_deg)
    closed = (d_unwrapped >= ivc) & (d_unwrapped <= evo + cycle_deg)
    p[closed] = np.interp(d_unwrapped[closed] * DEG, th_c, p_closed)
    blow = (th_deg > evo) & (th_deg <= evo + bd)
    w = 0.5 * (1 - np.cos(math.pi * (th_deg[blow] - evo) / bd))
    p[blow] = p_evo + (p_exh - p_evo) * w
    exh = (th_deg > evo + bd) & (th_deg < ivo)
    p[exh] = p_exh
    ovl = (th_deg >= ivo) & (th_deg <= evc)
    p[ovl] = p_exh + (p_int - p_exh) * (th_deg[ovl] - ivo) / (evc - ivo)
    intake = (th_deg > evc) & (th_deg < ivc)
    p[intake] = p_int

    meta = {"op": op, "q_total_J": q_total, "m_air_kg": m_air, "placeholder": True}
    return PressureTrace.from_degrees(th_deg, p, source=PLACEHOLDER_SOURCE, meta=meta)
