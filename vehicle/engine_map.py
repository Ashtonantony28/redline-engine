"""Engine maps for the vehicle model: WOT envelope, fuel flow (Willans fit to the M7 maps), motoring.

* Full-load envelope T_max(N): linear interpolation of the committed fixture of the M4 sim.wot
  sweep (fixtures/wot/wot_coarse.csv, `python -m vehicle.wot_fixture`).
* Fuel flow: the M7 speed-load maps (fixtures/maps/maps_m4.csv) cover 1000-7000 rpm but only
  MAP >= 0.5-0.8 bar (the low-speed part-load points were skipped in M7, see its commit message),
  while a drive cycle lives mostly below that. Each map point gives a fuel mean effective pressure
  p_f = mdot_f LHV / (V_d N / 120). At fixed speed p_f is very nearly affine in BMEP (the Willans
  line, [G&S] sec. 3.2): p_f = a(N) + b(N) BMEP, with b ~ 1 / (indicated efficiency) and the
  x-intercept -a/b = -(friction + pumping) MEP. a(N), b(N) are least-squares polynomials in N over
  every throttled map point (WOT points are knock limited with retarded spark and leave the line).
  The fit residual on every fixture point is reported (`fit_report`), and the line extrapolates to
  light load, idle and motoring:
* Motoring (overrun fuel cut): brake torque at zero fuel = the Willans x-intercept.

All numbers are MODELLED (no physical engine).
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from redline.spec import ROOT, geometry, load
from vehicle.params import RAD_S_PER_RPM, vehicle_spec

PA_PER_BAR = 1.0e5
J_PER_KWH = 3.6e6
G_PER_KG = 1.0e3
SEC_PER_MIN = 60.0
REV_PER_CYCLE = 2.0          # four-stroke


def _read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


@dataclass(frozen=True)
class EngineMap:
    wot_rpm: np.ndarray
    wot_torque: np.ndarray
    a_coef: np.ndarray        # p_f intercept polynomial in n = N / N_redline [Pa]
    b_coef: np.ndarray        # p_f slope polynomial in n [-]
    displacement: float
    lhv: float
    idle_rpm: float
    redline_rpm: float
    map_rows: tuple

    # ---- geometry helpers -------------------------------------------------------------------
    def bmep(self, torque):
        return np.asarray(torque, dtype=float) * 2.0 * math.pi * REV_PER_CYCLE / self.displacement

    def torque_from_bmep(self, bmep):
        return np.asarray(bmep, dtype=float) * self.displacement / (2.0 * math.pi * REV_PER_CYCLE)

    def cycles_per_s(self, rpm):
        return np.asarray(rpm, dtype=float) / (SEC_PER_MIN * REV_PER_CYCLE)

    # ---- envelope ---------------------------------------------------------------------------
    def max_torque(self, rpm):
        """WOT brake torque [N m]; flat below the first fixture speed, zero above redline."""
        rpm = np.asarray(rpm, dtype=float)
        t = np.interp(rpm, self.wot_rpm, self.wot_torque)
        return np.where(rpm > self.redline_rpm, 0.0, t)

    def max_power(self, rpm):
        return self.max_torque(rpm) * np.asarray(rpm, dtype=float) * RAD_S_PER_RPM

    # ---- Willans fuel model -----------------------------------------------------------------
    def _n(self, rpm):
        return np.asarray(rpm, dtype=float) / self.redline_rpm

    def willans(self, rpm):
        """(a [Pa], b [-]) of p_f = a + b BMEP at this speed."""
        n = self._n(rpm)
        return np.polyval(self.a_coef[::-1], n), np.polyval(self.b_coef[::-1], n)

    def fuel_rate(self, rpm, torque):
        """Fuel mass flow [kg/s] at brake torque [N m] (>= 0; zero at or below the motoring torque)."""
        a, b = self.willans(rpm)
        pf = np.maximum(a + b * self.bmep(torque), 0.0)
        return pf * self.displacement * self.cycles_per_s(rpm) / self.lhv

    def motoring_torque(self, rpm):
        """Closed-throttle zero-fuel brake torque [N m] (negative): the Willans x-intercept."""
        a, b = self.willans(rpm)
        return self.torque_from_bmep(-a / b)

    def bsfc(self, rpm, torque):
        """Brake specific fuel consumption [g/kWh] (inf at zero or negative torque)."""
        p = np.asarray(torque, dtype=float) * np.asarray(rpm, dtype=float) * RAD_S_PER_RPM
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(p > 0.0, self.fuel_rate(rpm, torque) / p * J_PER_KWH * G_PER_KG, np.inf)

    # ---- diagnostics ------------------------------------------------------------------------
    def in_map_coverage(self, rpm, torque):
        """True where (rpm, torque) lies inside the speed and BMEP range the M7 fixture actually sampled
        (lowest sampled BMEP at each fixture speed, linear in between); elsewhere the Willans line extrapolates."""
        speeds = sorted({float(r["rpm"]) for r in self.map_rows})
        bmin = [min(float(r["bmep_bar"]) for r in self.map_rows if float(r["rpm"]) == n) * PA_PER_BAR for n in speeds]
        rpm = np.asarray(rpm, dtype=float)
        inside = (rpm >= speeds[0]) & (rpm <= speeds[-1])
        return inside & (self.bmep(torque) >= np.interp(rpm, speeds, bmin)) & (np.asarray(torque) <= self.max_torque(rpm))

    def fit_report(self) -> dict:
        """Relative fuel-flow error of the Willans fit on every M7 fixture point."""
        out = {"throttled": [], "wot": []}
        for r in self.map_rows:
            rpm, T = float(r["rpm"]), float(r["torque_nm"])
            mf = float(r["bsfc_g_per_kwh"]) * float(r["power_kw"]) * 1.0e3 / (J_PER_KWH * G_PER_KG)
            err = float(self.fuel_rate(rpm, T) / mf - 1.0)
            out["wot" if r["level"] == "wot" else "throttled"].append({"rpm": rpm, "level": r["level"], "rel_err": err})
        for k in ("throttled", "wot"):
            e = [abs(x["rel_err"]) for x in out[k]]
            out[f"{k}_max_abs_rel_err"] = max(e) if e else float("nan")
            out[f"{k}_rms_rel_err"] = float(np.sqrt(np.mean(np.square(e)))) if e else float("nan")
        return out


def fuel_mep_points(rows, displacement, lhv):
    """(rpm, bmep [Pa], fuel MEP [Pa]) arrays from M7 map rows."""
    rpm = np.array([float(r["rpm"]) for r in rows])
    bmep = np.array([float(r["bmep_bar"]) for r in rows]) * PA_PER_BAR
    mf = np.array([float(r["bsfc_g_per_kwh"]) * float(r["power_kw"]) for r in rows]) * 1.0e3 / (J_PER_KWH * G_PER_KG)
    pf = mf * lhv / (displacement * rpm / (SEC_PER_MIN * REV_PER_CYCLE))
    return rpm, bmep, pf


def fit_willans(rows, displacement, lhv, redline_rpm, a_order, b_order):
    """Least-squares coefficients (ascending powers of n = N / N_redline) of a(n) and b(n)."""
    rpm, bmep, pf = fuel_mep_points(rows, displacement, lhv)
    n = rpm / redline_rpm
    cols = [n**k for k in range(a_order + 1)] + [bmep * n**k for k in range(b_order + 1)]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, pf, rcond=None)
    return coef[:a_order + 1], coef[a_order + 1:]


@lru_cache(maxsize=None)
def engine_map() -> EngineMap:
    vs = vehicle_spec()
    es = load()
    sim, ops = vs["simulation"], vs["engine_ops"]
    wot = _read_csv(ROOT / sim["wot_fixture"])
    rows = _read_csv(ROOT / sim["maps_fixture"])
    g = geometry(es)
    lhv = float(es["combustion"]["lhv"])
    redline = float(es["engine"]["redline_rpm"])
    fit_rows = [r for r in rows if r["level"] not in ops["willans_exclude_levels"]]
    a, b = fit_willans(fit_rows, g.displacement, lhv, redline, ops["willans_fit"]["intercept_poly_order"],
                       ops["willans_fit"]["slope_poly_order"])
    return EngineMap(np.array([float(r["rpm"]) for r in wot]), np.array([float(r["torque_nm"]) for r in wot]),
                     a, b, g.displacement, lhv, float(es["engine"]["idle_rpm"]), redline, tuple(rows))
