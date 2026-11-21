"""gymnasium environment around the MVEM + host ECU on a dyno, with random pedal steps (tip-in /
tip-out). The action is a fuel trim on top of the classical ECU fueling; the reward penalises the
lambda error of every combustion event in the control period. MODELLED plant, seeded."""
from __future__ import annotations

import math

import gymnasium as gym
import numpy as np

from hil import runner
from hil.session import Session
from ml import ml_spec
from ml.fuel.supervisor import Supervisor
from mvem import params

OBS_NAMES = ["pedal", "pedal_rate", "map", "map_rate", "lambda_sensor_err", "rpm", "throttle_err",
             "stft", "trim_prev"]
# observation scaling to order one (documented normalisers, not calibration)
PEDAL_RATE_SCALE = 0.1          # per s: a full pedal sweep in 0.1 s reads 1
MAP_SCALE = 1.0e-5              # 1 / Pa
MAP_RATE_SCALE = 1.0e-6         # s / Pa
LAMBDA_ERR_SCALE = 10.0
THROTTLE_ERR_SCALE = 1.0 / 30.0     # 1 / deg


def observe(sess: Session, prev: dict, trim_prev: float, period: float, cfg: dict) -> np.ndarray:
    eng, out = sess.eng, sess.out
    redline = float(params.engine_spec()["engine"]["redline_rpm"])
    return np.array([
        sess.pedal,
        (sess.pedal - prev["pedal"]) / period * PEDAL_RATE_SCALE,
        eng.p_man * MAP_SCALE,
        (eng.p_man - prev["map"]) / period * MAP_RATE_SCALE,
        (eng.lam_sensor - 1.0) * LAMBDA_ERR_SCALE,
        eng.rpm() / redline,
        (out.throttle_target_deg - eng.thr) * THROTTLE_ERR_SCALE,
        out.stft * LAMBDA_ERR_SCALE,
        trim_prev / float(cfg["trim_limit"]),
    ], dtype=np.float32)


def random_pedal(rng: np.random.Generator, cfg: dict):
    """Piecewise-linear pedal history: hold, then ramps between random levels."""
    lo, hi = (float(v) for v in cfg["pedal_levels"])
    level = rng.uniform(lo, hi)
    t = float(cfg["settle_s"])
    pts = [[0.0, level], [t, level]]
    while t < float(cfg["episode_s"]):
        new = rng.uniform(lo, hi)
        while abs(new - level) < float(cfg["min_pedal_step"]):
            new = rng.uniform(lo, hi)
        t_ramp = rng.uniform(*(float(v) for v in cfg["pedal_ramp_s"]))
        pts.append([t + t_ramp, new])
        t += t_ramp + rng.uniform(*(float(v) for v in cfg["step_interval_s"]))
        pts.append([t, new])
        level = new
    return pts


class Loop:
    """Shared control loop: `period` of plant steps per decision, used by the environment and by the
    evaluation, so training and evaluation see the same observation and actuation."""

    def __init__(self, scn: runner.Scenario, cfg: dict, ecu=None, supervised: bool = True):
        self.cfg = cfg
        self.sup = Supervisor(cfg, switching=supervised)
        self.sess = Session(scn, ecu)
        self.period = float(cfg["control_period_s"])
        self.n_sub = int(round(self.period / self.sess.dt))
        self.trim = 0.0
        self.sess.step()
        self.prev = {"pedal": self.sess.pedal, "map": self.sess.eng.p_man}

    def obs(self) -> np.ndarray:
        return observe(self.sess, self.prev, self.trim, self.period, self.cfg)

    def advance(self, action: float):
        """Policy output `action` in [-1, 1] for one control period, through the supervisor (limit,
        washout, switch). Returns the lambdas of the events that burned in the period."""
        s = self.sess
        self.prev = {"pedal": s.pedal, "map": s.eng.p_man}
        want = float(np.clip(action, -1.0, 1.0)) * self.sup.limit
        self.trim = s.fuel_trim = self.sup.filter(want, s.eng.lam_sensor, self.period)
        n0 = len(s.eng.events)
        for _ in range(self.n_sub):
            if s.done:
                break
            s.step()
        return [e["lambda"] for e in s.eng.events[n0:] if math.isfinite(e["lambda"])]


class FuelEnv(gym.Env):
    metadata = {"render_modes": []}

    def __init__(self, cfg: dict | None = None, ecu=None):
        self.cfg = cfg or ml_spec()["fuel_ppo"]
        self.ecu = ecu
        self.observation_space = gym.spaces.Box(-np.inf, np.inf, (len(OBS_NAMES),), np.float32)
        self.action_space = gym.spaces.Box(-1.0, 1.0, (1,), np.float32)
        self.loop = None

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        cfg = self.cfg
        rpm = float(self.np_random.choice(cfg["rpm_choices"]))
        scn = runner.Scenario(duration=float(cfg["episode_s"]), pedal=random_pedal(self.np_random, cfg),
                              dyno_rpm=[[0.0, rpm]], name="ppo")
        self.loop = Loop(scn, cfg, self.ecu)
        return self.loop.obs(), {}

    def step(self, action):
        cfg = self.cfg
        lams = self.loop.advance(float(action[0]))
        trim = self.loop.trim
        # start-up (ECU syncing, first fuel) is not the controller's business: no lambda reward there
        scored = lams if self.loop.sess.eng.t >= float(cfg["reward_start_s"]) else []
        clip = float(cfg["reward_error_clip"])
        err2 = float(np.mean([min(abs(lam - 1.0), clip) ** 2 for lam in scored])) if scored else 0.0
        reward = -float(cfg["reward_lambda_weight"]) * err2 \
                 - float(cfg["reward_trim_weight"]) * (trim / float(cfg["trim_limit"])) ** 2 * self.loop.period
        return self.loop.obs(), reward, False, self.loop.sess.done, {"lambdas": lams}
