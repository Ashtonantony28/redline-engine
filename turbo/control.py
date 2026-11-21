"""Wastegate PI boost controller and turbocharger shaft dynamics (M10 turbo stretch). MODELLED.

The controller runs in Python, not in the C11 ECU: the ECU has no wastegate output, its MAP sensor
range and table axis end at 1.15 bar (ecu_calib_gen.h ECU_MAP_MAX), and its output struct is shared bit
for bit with the RV64 co-simulation, so adding a boost channel there is not a clean, contained change.
"""
from __future__ import annotations

import math

import numpy as np

SEC_PER_MIN = 60.0


def rpm_to_rad(rpm: float) -> float:
    return rpm * 2.0 * math.pi / SEC_PER_MIN


def rad_to_rpm(omega: float) -> float:
    return omega * SEC_PER_MIN / (2.0 * math.pi)


class TurboShaft:
    """J omega d(omega)/dt = eta_m P_T - P_C, integrated on the rotor kinetic energy (explicit Euler)."""

    def __init__(self, cfg: dict, inertia_factor: float = 1.0):
        self.J = float(cfg["inertia"]) * inertia_factor
        self.eta_m = float(cfg["mechanical_efficiency"])
        self.w_max = rpm_to_rad(float(cfg["max_speed_rpm"]))
        self.w_min = rpm_to_rad(float(cfg["min_speed_rpm"]))
        self.omega = self.w_min
        self.overspeed = False

    def net_power(self, P_T: float, P_C: float) -> float:
        return self.eta_m * P_T - P_C

    def step(self, P_T: float, P_C: float, dt: float) -> float:
        E = 0.5 * self.J * self.omega ** 2 + self.net_power(P_T, P_C) * dt
        w = math.sqrt(max(2.0 * E / self.J, 0.0))
        if w > self.w_max:
            self.overspeed = True
        self.omega = min(max(w, self.w_min), self.w_max)
        return self.omega


class BoostPI:
    """Wastegate opening (0 closed .. 1 open) from the filtered intake-plenum pressure.

    opening_cmd = clamp(i - kp e), di/dt = -ki e, e = (target - MAP_f) / target; the integrator is frozen
    while the command is saturated in the direction the error pushes (conditional integration). With
    target None (no boost demand, part load) the wastegate is held at `no_demand_opening` and the integrator
    tracks it, so enabling the controller is bumpless. A first-order actuator lag follows the command."""

    def __init__(self, cfg: dict, wg_cfg: dict, rpm_axis=None, target_axis=None):
        self.kp, self.ki = float(cfg["kp"]), float(cfg["ki"])
        self.tau_f = float(cfg["map_filter_tau"])
        self.tau_a = float(wg_cfg["actuator_tau"])
        self.rpm_axis = np.asarray(rpm_axis if rpm_axis is not None else cfg["target_rpm"], dtype=float)
        self.map_axis = np.asarray(target_axis if target_axis is not None else cfg["target_map"], dtype=float)
        self.i = float(cfg["opening_init"])
        self.no_demand = float(cfg["no_demand_opening"])
        self.opening = self.i
        self.map_f = None
        self.target = None

    def target_for(self, rpm: float) -> float:
        return float(np.interp(rpm, self.rpm_axis, self.map_axis))

    def step(self, map_pa: float, dt: float) -> float:
        a = dt / (self.tau_f + dt)
        self.map_f = map_pa if self.map_f is None else self.map_f + a * (map_pa - self.map_f)
        if self.target is None:
            cmd = self.no_demand
            self.i = cmd
        else:
            e = (self.target - self.map_f) / self.target
            raw = self.i - self.kp * e
            cmd = min(max(raw, 0.0), 1.0)
            sat_lo = raw <= 0.0 and e > 0.0
            sat_hi = raw >= 1.0 and e < 0.0
            if not (sat_lo or sat_hi):
                self.i = min(max(self.i - self.ki * e * dt, 0.0), 1.0)
        b = dt / (self.tau_a + dt)
        self.opening += b * (cmd - self.opening)
        return self.opening
