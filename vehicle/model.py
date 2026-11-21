"""Longitudinal vehicle model: road load, rotating inertia, shift strategy, backward drive-cycle
simulation and forward full-throttle simulation. MODELLED (spec/vehicle.yaml, spec/engine.yaml).

Road load at speed v [m/s] on grade alpha ([G&S] eq. 2.2-2.4, [Gil] ch. 4):
    F = m g (c_r cos(alpha) + sin(alpha)) + m g c_r2 v^2 cos(alpha) + 1/2 rho c_d A v^2
Equation of motion at the wheel: m_eq dv/dt = F_traction - F_road, with m_eq from
vehicle.params.equivalent_mass (wheels always, engine side only with the clutch engaged).
Driveline: constant efficiency eta on positive torque (engine to wheel) and on negative torque
(wheel to engine, engine braking), engine speed locked to wheel speed except during launch slip.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from vehicle import params as vp
from vehicle.engine_map import EngineMap, engine_map

RAD = vp.RAD_S_PER_RPM


# ----------------------------------------------------------------------------------------------
# road load
# ----------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RoadLoad:
    mass: float
    g: float
    c_r: float
    c_r2: float
    rho: float
    cd_a: float
    grade: float

    @classmethod
    def from_spec(cls, mass: float, vs: dict | None = None) -> "RoadLoad":
        v = (vs or vp.vehicle_spec())["vehicle"]
        return cls(mass, v["gravity"], v["rolling_resistance"], v["rolling_resistance_speed"], v["air_density"],
                   v["drag_coefficient"] * v["frontal_area"], v["road_grade"])

    def force(self, v):
        v = np.asarray(v, dtype=float)
        n = self.mass * self.g * math.cos(self.grade)
        return (n * (self.c_r + self.c_r2 * v**2) + self.mass * self.g * math.sin(self.grade)
                + 0.5 * self.rho * self.cd_a * v**2)

    def power(self, v):
        return self.force(v) * np.asarray(v, dtype=float)

    @property
    def coast_coefficients(self):
        """(A, C) of F = A + C v^2 on a flat road (c_r2 folded into C)."""
        n = self.mass * self.g * math.cos(self.grade)
        return n * self.c_r + self.mass * self.g * math.sin(self.grade), n * self.c_r2 + 0.5 * self.rho * self.cd_a


def coastdown_time_analytic(A: float, C: float, m_eq: float, v0: float, v1: float) -> float:
    """Time to coast from v0 to v1 under m dv/dt = -(A + C v^2):
    t = sqrt(m^2 / (A C)) [atan(v0 sqrt(C/A)) - atan(v1 sqrt(C/A))]."""
    k = math.sqrt(C / A)
    return m_eq / math.sqrt(A * C) * (math.atan(v0 * k) - math.atan(v1 * k))


def coastdown_numeric(road: RoadLoad, m_eq: float, v0: float, v1: float, dt: float) -> float:
    """RK4 integration of the same coast-down (neutral, no engine), linear crossing interpolation."""
    def acc(v):
        return -float(road.force(v)) / m_eq
    t, v = 0.0, v0
    while v > v1:
        k1 = acc(v)
        k2 = acc(v + 0.5 * dt * k1)
        k3 = acc(v + 0.5 * dt * k2)
        k4 = acc(v + dt * k3)
        vn = v + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if vn <= v1:
            return t + dt * (v - v1) / (v - vn)
        t, v = t + dt, vn
    return t


# ----------------------------------------------------------------------------------------------
# powertrain
# ----------------------------------------------------------------------------------------------
@dataclass
class Powertrain:
    vs: dict = field(default_factory=vp.vehicle_spec)
    es: dict = field(default_factory=vp.engine_spec)
    emap: EngineMap = field(default_factory=engine_map)

    def __post_init__(self):
        d = self.vs["driveline"]
        self.n_gears = len(d["gear_ratios"])
        self.eta = d["efficiency"]
        self.r_w = self.vs["vehicle"]["wheel_radius"]
        self.idle = float(self.es["engine"]["idle_rpm"])
        self.redline = float(self.es["engine"]["redline_rpm"])
        self.p_acc = self.vs["engine_ops"]["idle_accessory_power"]

    def ratio(self, gear: int) -> float:
        return vp.overall_ratio(gear, self.vs)

    def rpm(self, v: float, gear: int) -> float:
        return vp.engine_rpm(v, gear, self.vs)

    def accessory_torque(self, rpm: float) -> float:
        return self.p_acc / (max(rpm, self.idle) * RAD)

    def wheel_force_wot(self, v: float, gear: int) -> float:
        n = self.rpm(v, gear)
        return float(self.emap.max_torque(n)) * self.ratio(gear) * self.eta / self.r_w

    def traction_limit(self, mass: float) -> float:
        """Peak tractive force of the driven axle with longitudinal load transfer (front drive):
        F <= mu (m g s - F h / L)  =>  F <= mu m g s / (1 + mu h / L)  (s = static share)."""
        v = self.vs["vehicle"]
        mu, s, h, L = v["tyre_friction"], v["driven_axle_static_share"], v["cg_height"], v["wheelbase"]
        sign = 1.0 if v["driven_axle"] == "front" else -1.0
        return mu * mass * v["gravity"] * s / (1.0 + sign * mu * h / L)

    # ---- shift strategy -------------------------------------------------------------------
    def select_gear(self, v: float, p_wheel: float, current: int | None) -> int:
        """Economy strategy: highest gear whose engine speed stays above the lugging limit and whose
        WOT power covers the demand with the reserve; hysteresis keeps the current gear a little
        longer on the way down. Falls back to the gear with the most WOT power."""
        sh = self.vs["shift"]
        p_eng = max(p_wheel, 0.0) / self.eta
        best, best_p = 1, -1.0
        for g in range(self.n_gears, 0, -1):
            n = self.rpm(v, g)
            if n > self.redline:
                continue
            n_min = sh["min_rpm_drive"] - (sh["downshift_hysteresis_rpm"] if g == current else 0.0)
            if current is not None and g > current:
                n_min = max(n_min, sh["min_rpm_upshift_target"])
            p_avail = float(self.emap.max_power(n))
            if p_avail > best_p:
                best, best_p = g, p_avail
            if n >= n_min and p_eng + self.p_acc <= sh["power_reserve"] * p_avail:
                return g
        return best


# ----------------------------------------------------------------------------------------------
# backward drive-cycle simulation
# ----------------------------------------------------------------------------------------------
@dataclass
class CycleResult:
    t: np.ndarray            # interval mid times [s]
    v: np.ndarray            # interval mean speed [m/s]
    a: np.ndarray            # interval acceleration [m/s^2]
    gear: np.ndarray         # 0 = clutch open (standstill / idle)
    rpm: np.ndarray
    torque: np.ndarray       # engine brake torque [N m] (incl. accessory load)
    fuel_rate: np.ndarray    # kg/s
    mode: list               # per interval: idle, launch, drive, fuel_cut, overrun_idle
    shortfall: np.ndarray    # N, traction demand the engine could not deliver
    p_wheel: np.ndarray      # W at the wheel (signed)
    dt: float

    @property
    def distance(self) -> float:
        return float(np.sum(self.v) * self.dt)

    @property
    def fuel_mass(self) -> float:
        return float(np.sum(self.fuel_rate) * self.dt)


def simulate_cycle(t_s: np.ndarray, v_ms: np.ndarray, pt: Powertrain | None = None,
                   mass: float | None = None) -> CycleResult:
    """Backward (demand-driven) simulation over each interval [t_k, t_k+1] of a speed trace."""
    pt = pt or Powertrain()
    vs = pt.vs
    mass = mass if mass is not None else vp.wltp_test_mass(vs)
    road = RoadLoad.from_spec(mass, vs)
    sh, ops, lau = vs["shift"], vs["engine_ops"], vs["driveline"]["launch"]
    emap = pt.emap
    dt_arr = np.diff(t_s)
    dt = float(dt_arr[0])
    if not np.allclose(dt_arr, dt):
        raise ValueError("speed trace must be uniformly sampled")
    vm = 0.5 * (v_ms[1:] + v_ms[:-1])
    acc = np.diff(v_ms) / dt
    n = len(vm)
    out = {k: np.zeros(n) for k in ("gear", "rpm", "torque", "fuel", "short", "pw")}
    modes = []
    gear = None
    for k in range(n):
        v, a = vm[k], acc[k]
        if v < sh["standstill_speed"]:
            gear = None
            rpm = pt.idle
            T = pt.accessory_torque(rpm)
            modes.append("idle")
            out["rpm"][k], out["torque"][k], out["fuel"][k] = rpm, T, emap.fuel_rate(rpm, T)
            continue
        f_road = float(road.force(v))
        p_wheel_quasi = (mass * a + f_road) * v
        g_sel = pt.select_gear(v, p_wheel_quasi, gear)
        launch = v < sh["first_gear_max_speed"] and pt.rpm(v, 1) < lau["engine_rpm"] and a > 0.0
        if launch:
            g_sel = 1
        locked_rpm = pt.rpm(v, g_sel)
        if launch:
            rpm = max(lau["engine_rpm"], locked_rpm)
            m_eq = vp.equivalent_mass(mass, None, vs, pt.es)
        else:
            rpm = locked_rpm
            m_eq = vp.equivalent_mass(mass, g_sel, vs, pt.es)
        # below idle in gear (slow creeping without launch flag): clutch slips at idle speed
        if rpm < pt.idle:
            rpm = pt.idle
            m_eq = vp.equivalent_mass(mass, None, vs, pt.es)
        gear = g_sel
        F = m_eq * a + f_road
        T_wheel = F * pt.r_w
        T_acc = pt.accessory_torque(rpm)
        if T_wheel >= 0.0:
            T = T_wheel / (pt.ratio(gear) * pt.eta) + T_acc
            T_max = float(emap.max_torque(rpm))
            if T > T_max:
                out["short"][k] = (T - T_max) * pt.ratio(gear) * pt.eta / pt.r_w
                T = T_max
            fuel = float(emap.fuel_rate(rpm, T))
            mode = "launch" if launch else "drive"
        else:
            T_demand = T_wheel * pt.eta / pt.ratio(gear) + T_acc
            T_motor = float(emap.motoring_torque(rpm))
            if rpm > ops["fuel_cut_min_rpm"] and T_demand < 0.0:
                T, fuel, mode = max(T_demand, T_motor), 0.0, "fuel_cut"      # brakes take the rest
            else:
                T = max(T_demand, T_acc)
                fuel, mode = float(emap.fuel_rate(rpm, max(T, T_acc))), "overrun_idle"
        modes.append(mode)
        out["gear"][k], out["rpm"][k], out["torque"][k], out["fuel"][k], out["pw"][k] = gear, rpm, T, fuel, F * v
    tm = 0.5 * (t_s[1:] + t_s[:-1])
    return CycleResult(tm, vm, acc, out["gear"].astype(int), out["rpm"], out["torque"], out["fuel"], modes,
                       out["short"], out["pw"], dt)


# ----------------------------------------------------------------------------------------------
# forward full-throttle simulation
# ----------------------------------------------------------------------------------------------
@dataclass
class AccelResult:
    t: np.ndarray
    v: np.ndarray
    gear: np.ndarray
    rpm: np.ndarray
    shifts: list

    def time_between(self, v0: float, v1: float) -> float:
        """First time the speed crosses v1 minus first time it crosses v0 (linear interpolation)."""
        def cross(x):
            i = int(np.argmax(self.v >= x))
            if self.v[i] < x:
                return float("nan")
            if i == 0:
                return float(self.t[0])
            return float(self.t[i - 1] + (x - self.v[i - 1]) / (self.v[i] - self.v[i - 1]) * (self.t[i] - self.t[i - 1]))
        return cross(v1) - cross(v0)


def simulate_wot(v0: float, v_end: float, gear0: int, pt: Powertrain | None = None, mass: float | None = None,
                 allow_shifts: bool = True, launch: bool = True) -> AccelResult:
    """Full-throttle run from v0 to v_end, explicit RK2 in time. A standing start slips the clutch at
    the launch engine speed until the wheels catch up; upshifts at the accel shift speed interrupt
    traction for the shift time (the car coasts on road load). Traction is capped by the tyre limit."""
    pt = pt or Powertrain()
    vs = pt.vs
    mass = mass if mass is not None else vs["vehicle"]["mass_running_order"]
    road = RoadLoad.from_spec(mass, vs)
    sim, sh, d = vs["simulation"], vs["shift"], vs["driveline"]
    dt, t_max = sim["accel_dt"], sim["accel_max_time"]
    F_trac = pt.traction_limit(mass)
    gear, v, t = gear0, v0, 0.0
    shift_left = 0.0
    T, V, G, R, shifts = [t], [v], [gear], [pt.rpm(v, gear)], []
    slip_rpm = d["launch"]["wot_engine_rpm"]

    def force(vv, gg, in_shift):
        if in_shift:
            return -float(road.force(vv)), vp.equivalent_mass(mass, None, vs, pt.es), pt.rpm(vv, gg)
        locked = pt.rpm(vv, gg)
        slipping = launch and locked < slip_rpm and gg == 1
        n = slip_rpm if slipping else locked
        Fw = min(float(pt.emap.max_torque(n)) * pt.ratio(gg) * pt.eta / pt.r_w, F_trac)
        m_eq = vp.equivalent_mass(mass, None if slipping else gg, vs, pt.es)
        return Fw - float(road.force(vv)), m_eq, n

    while v < v_end and t < t_max:
        in_shift = shift_left > 0.0
        f1, m1, n1 = force(v, gear, in_shift)
        a1 = f1 / m1
        f2, m2, _ = force(v + dt * a1, gear, in_shift)
        v += 0.5 * dt * (a1 + f2 / m2)
        t += dt
        shift_left = max(shift_left - dt, 0.0)
        if allow_shifts and not in_shift and gear < pt.n_gears and pt.rpm(v, gear) >= sh["accel_upshift_rpm"]:
            gear += 1
            shift_left = d["shift_time"]
            shifts.append((t, gear))
        T.append(t)
        V.append(v)
        G.append(gear)
        R.append(n1)
    return AccelResult(np.array(T), np.array(V), np.array(G), np.array(R), shifts)
