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


# ------------------------------------------------------------------------------------------------
# REAL pressure sources (M10): M4 1D engine at full load, thermo/ cycle at idle
# ------------------------------------------------------------------------------------------------
M4_SOURCE = "M4 1D gas-dynamics engine (gasdyn.engine, WOT, two-zone Wiebe closed part)"
THERMO_IDLE_SOURCE = "M2 quasi-steady two-zone cycle (thermo.cycle), intake pressure solved for zero brake torque"


def _recording_engine():
    """gasdyn EngineModel subclass that keeps each cylinder's latest closed-part result (the
    pressure between IVC and EVO is computed by the thermo two-zone model outside the 1D network,
    so the network trace has NaN there)."""
    from gasdyn.engine import EngineModel

    class RecordingEngine(EngineModel):
        def _closed(self, c, m, x, T, rpm, st):
            cl = super()._closed(c, m, x, T, rpm, st)
            if not hasattr(self, "last_closed"):
                self.last_closed = {}
            self.last_closed[c] = cl
            return cl

    return RecordingEngine


def m4_cylinder_traces(rpm: float, throttle_deg: float | None = None, dspec: dict | None = None):
    """Per-cylinder p(local theta) of the converged M4 run: (theta_deg grid, p [ncyl, n], metrics).

    Open part (valve off its seat): the 1D-network trace binned at spec trace_step_deg, global crank angle
    shifted by each cylinder's firing phase. Closed part (seat IVC .. seat EVO): the thermo two-zone result
    computed inside M4 at the last IVC (the converged periodic state)."""
    from thermo import twozone as tz

    dspec = dspec or params.dyn_spec()
    eng = _recording_engine()()
    op = eng.run(float(rpm), throttle_deg=throttle_deg, trace=True)
    cycle = params.FULL_CYCLE / DEG
    n = dspec["pressure_m4"]["grid_points"]
    grid = np.linspace(0.0, cycle, n)
    p_all = np.empty((eng.ncyl, n))
    for c in range(eng.ncyl):
        raw = op.trace["cylinder_p"][:, c]
        ok = ~np.isnan(raw)
        loc = np.mod(op.trace_deg[ok] - eng.phase_deg[c], cycle)
        order = np.argsort(loc)
        loc, po = loc[order], raw[ok][order]
        # periodic interpolation of the open-part samples
        p_open = np.interp(grid, np.concatenate([loc - cycle, loc, loc + cycle]), np.tile(po, 3))
        cl = eng.last_closed[c]
        th_c = np.asarray(cl.theta, dtype=float)          # IVC - 720 .. EVO, local deg
        p_c = cl.col(tz.O_P)
        g_unwrapped = np.where(grid > eng.ivc_deg, grid - cycle, grid)
        closed = (g_unwrapped >= th_c[0]) & (g_unwrapped <= th_c[-1])
        p_open[closed] = np.interp(g_unwrapped[closed], th_c, p_c)
        p_all[c] = p_open
    return grid, p_all, op.metrics


def m4_trace(rpm: float, throttle_deg: float | None = None, dspec: dict | None = None) -> PressureTrace:
    """REAL full-load (or throttled) pressure trace from the M4 engine at `rpm` (see m4_cylinder_traces).
    With spec pressure_m4.cylinder_average the four aligned cylinders are averaged (the dynamics model
    assumes identical cylinders); otherwise cylinder 1 is used."""
    dspec = dspec or params.dyn_spec()
    grid, p_all, met = m4_cylinder_traces(rpm, throttle_deg, dspec)
    p = p_all.mean(axis=0) if dspec["pressure_m4"]["cylinder_average"] else p_all[0]
    spread = float(np.ptp(p_all.max(axis=1)) / p_all.max(axis=1).mean())
    meta = {"rpm": float(rpm), "placeholder": False, "torque_nm": met["torque_nm"], "imep_net_bar": met["imep_net_bar"],
            "p_max_bar_m4": met["p_max_bar"], "converged": bool(met["converged"]), "peak_spread_between_cylinders": spread}
    return PressureTrace.from_degrees(grid, p, source=M4_SOURCE, meta=meta)


def thermo_trace(rpm: float, p_intake: float, spec: dict | None = None, dspec: dict | None = None):
    """Quasi-steady thermo/ (M2) cycle at an intake pressure: (PressureTrace, metrics)."""
    import copy
    from thermo.cycle import run_cycle
    from thermo.props import thermo_spec

    dspec = dspec or params.dyn_spec()
    idl = dspec["pressure_m4"]["idle"]
    tspec = copy.deepcopy(thermo_spec())
    tspec["loads"]["dyn_idle"] = {"p_intake": float(p_intake), "T_intake": idl["T_intake"],
                                  "p_exhaust": idl["p_exhaust"], "T_exhaust": idl["T_exhaust"]}
    res = run_cycle(float(rpm), "dyn_idle", tspec=tspec, spec=spec)
    cycle = params.FULL_CYCLE / DEG
    th = np.mod(np.asarray(res.theta, dtype=float), cycle)
    order = np.argsort(th, kind="stable")
    th, p = th[order], np.asarray(res.p, dtype=float)[order]
    keep = np.concatenate([[True], np.diff(th) > 0])
    th, p = th[keep], p[keep]
    n = dspec["pressure_m4"]["grid_points"]
    grid = np.linspace(0.0, cycle, n)
    pg = np.interp(grid, np.concatenate([th - cycle, th, th + cycle]), np.tile(p, 3))
    return PressureTrace.from_degrees(grid, pg, source=THERMO_IDLE_SOURCE,
                                      meta={"rpm": float(rpm), "p_intake": float(p_intake), "placeholder": False,
                                            "imep_net_bar": res.metrics["imep_net_bar"]}), res.metrics


def idle_trace(rpm: float | None = None, spec: dict | None = None, dspec: dict | None = None) -> PressureTrace:
    """REAL idle trace: thermo/ cycle at the intake pressure where net IMEP equals the Chen-Flynn FMEP
    (spec/gasdyn.yaml friction), i.e. zero brake torque at idle speed."""
    from scipy.optimize import brentq
    from gasdyn.friction import fmep_total
    from gasdyn.spec import gasdyn_spec

    spec = spec or params.engine_spec()
    dspec = dspec or params.dyn_spec()
    rpm = float(rpm or spec["engine"]["idle_rpm"])
    fr = gasdyn_spec()["friction"]
    sp = 2.0 * spec["engine"]["stroke"] * rpm / 60.0
    idl = dspec["pressure_m4"]["idle"]

    def bmep(p_int):
        tr, met = thermo_trace(rpm, p_int, spec, dspec)
        return met["imep_net_bar"] * 1.0e5 - fmep_total(tr.peak, sp, fr)

    lo, hi = idl["p_intake_bracket"]
    p_int = brentq(bmep, float(lo), float(hi), xtol=float(idl["p_intake_tol"]))
    tr, met = thermo_trace(rpm, p_int, spec, dspec)
    tr.meta["fmep_bar"] = fmep_total(tr.peak, sp, fr) / 1.0e5
    return tr


class PressureProvider:
    """p(theta) per (operating point, rpm). 'placeholder': the single-zone Wiebe stand-in (ignores rpm);
    'm4': REAL traces, M4 at full load per speed and the thermo/ zero-torque cycle at idle. Cached."""

    def __init__(self, kind: str = "placeholder", spec: dict | None = None, dspec: dict | None = None):
        if kind not in ("placeholder", "m4"):
            raise ValueError(kind)
        self.kind = kind
        self.spec = spec or params.engine_spec()
        self.dspec = dspec or params.dyn_spec()
        self._cache: dict = {}

    def __call__(self, op: str, rpm: float) -> PressureTrace:
        key = (op, None if self.kind == "placeholder" else float(rpm))
        if key not in self._cache:
            if self.kind == "placeholder":
                self._cache[key] = placeholder_wiebe_trace(op, self.spec, self.dspec)
            elif op == "idle":
                self._cache[key] = idle_trace(rpm, self.spec, self.dspec)
            elif op == "full_load":
                self._cache[key] = m4_trace(rpm, None, self.dspec)
            else:
                raise ValueError(op)
        return self._cache[key]
