"""Simplex-style supervisor for the learned fueling trim (spec/ml.yaml fuel_ppo).

The classical ECU fueling (speed-density + x-tau + closed-loop lambda) is the safety controller and
is always in the loop; the learned policy only adds a multiplicative trim on top, through three
layers of protection:
  1. authority: |trim| <= trim_limit;
  2. washout: the applied trim is the policy output minus its own low-pass (trim_washout_s), so the
     policy has transient authority only. A constant output decays to zero, the closed-loop lambda
     trims (STFT / LTFT) never have to absorb a learned bias, and switching the policy off is not
     itself a fueling step;
  3. switch: while the measured lambda (wideband sensor) is outside `switch_band` the trim is
     forced to zero (pure classical fueling), and stays zero for `supervisor_hold_s` after it returns.
The safety claim checked in tests: per-event lambda stays inside `safe_band` in every evaluation
cycle, for the learned policy and for adversarial policies (constant, bang-bang, random)."""
from __future__ import annotations


class Supervisor:
    def __init__(self, cfg: dict, switching: bool = True):
        self.lo, self.hi = (float(v) for v in cfg["switch_band"])
        self.hold = float(cfg["supervisor_hold_s"])
        self.limit = float(cfg["trim_limit"])
        self.washout = float(cfg["trim_washout_s"])
        self.switching = switching
        self.reset()

    def reset(self):
        self.timer = 0.0
        self.lowpass = 0.0
        self.switches = 0
        self.time_on_classical = 0.0

    @property
    def classical(self) -> bool:
        return self.timer > 0.0

    def filter(self, trim: float, lam_sensor: float, dt: float) -> float:
        """Returns the trim to apply for the next `dt` seconds."""
        trim = min(max(trim, -self.limit), self.limit)
        self.lowpass += (trim - self.lowpass) * dt / (self.washout + dt)
        trim = min(max(trim - self.lowpass, -self.limit), self.limit)
        if self.switching and not self.lo <= lam_sensor <= self.hi:
            if not self.classical:
                self.switches += 1
            self.timer = self.hold
        if self.classical:
            self.timer = max(self.timer - dt, 0.0)
            self.time_on_classical += dt
            return 0.0
        return trim
