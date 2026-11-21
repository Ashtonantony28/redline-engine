"""Closed-loop host HIL: the MVEM plant (mvem.plant.Engine) drives the C ECU core through the
ctypes bridge. Crank/cam edges are delivered with exact interpolated timer ticks; the ECU's
output-compare events (injection, dwell, spark) are applied to the plant at their tick times; the
slow loop runs once per plant step. All results are MODELLED: there is no physical engine.
"""
from __future__ import annotations

import heapq
import math
from dataclasses import dataclass, field

import numpy as np
import yaml

from hil.bridge import HostEcu
from mvem import params
from mvem.plant import DEG_PER_RAD, Engine
from redline.spec import ROOT

CYCLES_PATH = ROOT / "hil" / "cycles.yaml"
INJ_END = -1          # runner-internal event: injector closes
LOG_KEYS = ["t", "rpm", "rpm_ecu", "map", "thr", "thr_target", "etc_v", "lambda_sensor", "pedal", "spark0",
            "stft", "ltft", "knock_retard_max", "torque_brake", "torque_ind", "closed_loop", "idle_active",
            "rev_cut", "limp", "sync", "pw0", "load"]


def load_cycles() -> dict:
    with open(CYCLES_PATH) as f:
        return yaml.safe_load(f)


def _hist(points, t):
    pts = np.asarray(points, dtype=float)
    return float(np.interp(t, pts[:, 0], pts[:, 1]))


@dataclass
class Scenario:
    duration: float
    pedal: list
    load: list = field(default_factory=lambda: [[0.0, 0.0]])
    dyno_rpm: list | None = None
    init_rpm: float | None = None
    init_map: float = 0.40e5
    hold_s: float = 0.4
    clt: float = 360.0
    knock_windows: list = field(default_factory=list)     # [[t0, t1], ...] force knock, all cylinders
    knock_offset_deg: float = 0.0                         # lowers the simulated knock limit
    tps2_fault_at: float | None = None                    # throttle sensor 2 jumps (disagreement)
    windows: dict = field(default_factory=dict)
    name: str = "custom"

    def __post_init__(self):
        if self.init_rpm is None:
            if self.dyno_rpm is None:
                raise ValueError("free-running scenario needs init_rpm")
            self.init_rpm = float(self.dyno_rpm[0][1])

    @classmethod
    def from_cycle(cls, name: str, **over) -> "Scenario":
        cy = load_cycles()
        d = dict(cy.get("defaults", {}))
        d.update(cy[name])
        d.update(over)
        acc = float(params.mvem_spec()["loads"]["accessory_step"])
        load = [[float(t), acc if v == "accessory_step" else float(v)] for t, v in d.get("load", [[0.0, 0.0]])]
        dyno = d.get("dyno_rpm")
        init_rpm = d.get("init_rpm", dyno[0][1] if dyno else None)
        return cls(duration=float(d["duration"]), pedal=d["pedal"], load=load, dyno_rpm=dyno,
                   init_rpm=float(init_rpm), init_map=float(d["init_map"]), hold_s=float(d["hold_s"]),
                   clt=float(d["clt"]), knock_windows=d.get("knock_windows", []),
                   tps2_fault_at=d.get("tps2_fault_at"), windows=d.get("windows", {}), name=name)


@dataclass
class Result:
    scenario: Scenario
    log: dict
    events: list
    status: dict

    def window(self, key, t0, t1):
        t = self.log["t"]
        m = (t >= t0) & (t < t1)
        return self.log[key][m]

    def event_lambdas(self, t0, t1, fired_only=True):
        return np.array([e["lambda"] for e in self.events
                         if t0 <= e["t"] < t1 and (e["fired"] or not fired_only) and math.isfinite(e["lambda"])])


def run(scn: Scenario, ecu: HostEcu | None = None) -> Result:
    ecu = ecu or HostEcu()
    ecu.reset()
    k = ecu.consts
    hz = float(k.timer_hz)
    ms = params.mvem_spec()
    dt = float(ms["sim"]["dt"])
    eng = Engine(scn.init_rpm, scn.init_map)
    eng.knock_offset_deg = scn.knock_offset_deg
    if scn.knock_windows:
        wins = [(float(a), float(b)) for a, b in scn.knock_windows]
        eng.knock_hook = lambda t, cyl: any(a <= t < b for a, b in wins)
    n = int(round(scn.duration / dt))
    log = {key: np.zeros(n) for key in LOG_KEYS}
    pending: list = []
    seq = 0
    ev_inj, ev_spark = k.ev_inj, k.ev_spark

    def apply(te, typ, cyl, arg, th0, om_deg, t0):
        nonlocal seq
        if typ == ev_inj:
            # the pulse's fuel reaches the port when the injector closes (at the ECU's EOI angle);
            # with long pulses the opening can precede the previous IVC, the fuel still belongs to
            # this cycle
            seq += 1
            heapq.heappush(pending, (te + arg / hz, seq, INJ_END, cyl, arg))
        elif typ == INJ_END:
            eng.inject(cyl, arg / hz)
        elif typ == ev_spark:
            eng.spark(cyl, th0 + om_deg * (te - t0))

    for step in range(n):
        t = eng.t
        pedal = _hist(scn.pedal, t)
        load = _hist(scn.load, t)
        if scn.dyno_rpm is not None:
            dyno = _hist(scn.dyno_rpm, t)
        else:
            dyno = scn.init_rpm if t < scn.hold_s else None
        tps1, tps2 = eng.tps()
        if scn.tps2_fault_at is not None and t >= scn.tps2_fault_at:
            tps2 = eng.thr_max / 2.0
        out = ecu.slow(int(round(t * hz)), eng.p_man, eng.T_man, scn.clt, tps1, tps2, pedal,
                       eng.lam_sensor, eng.vbatt)
        volts = out.etc_voltage

        th0 = eng.theta
        om_deg = eng.omega * DEG_PER_RAD
        if om_deg > 0.0:
            for thc, kind, idx in eng.crossings(th0, th0 + om_deg * dt):
                tc = t + (thc - th0) / om_deg
                while pending and pending[0][0] <= tc:
                    te, _, typ, cyl, arg = heapq.heappop(pending)
                    apply(te, typ, cyl, arg, th0, om_deg, t)
                if kind == "tooth" or kind == "cam":
                    tick = int(round(tc * hz))
                    if kind == "tooth":
                        ecu.crank(tick)
                    else:
                        ecu.cam(tick)
                    for (tk, typ, cyl, arg) in ecu.drain():
                        seq += 1
                        heapq.heappush(pending, (tk / hz, seq, typ, cyl, arg))
                elif kind == "ivc":
                    eng.ivc(idx, tc)
                else:
                    _, counts = eng.burn(idx, tc)
                    ecu.knock(idx, counts)
        while pending and pending[0][0] <= t + dt:
            te, _, typ, cyl, arg = heapq.heappop(pending)
            apply(te, typ, cyl, arg, th0, om_deg, t)
        eng.advance(dt, etc_volts=volts, load_torque=load, dyno_rpm=dyno)

        row = (t, eng.rpm(), out.rpm, eng.p_man, eng.thr, out.throttle_target_deg, volts, eng.lam_sensor, pedal,
               out.spark_deg[0], out.stft, out.ltft, max(out.knock_retard_deg), eng.T_brake, eng.T_ind,
               out.closed_loop, out.idle_active, out.rev_cut, out.limp, out.sync, out.pw_s[0], load)
        for key, v in zip(LOG_KEYS, row):
            log[key][step] = v
    return Result(scn, log, eng.events, ecu.status())


# --------------------------------------------------------------------------- metrics
def summary(res: Result) -> dict:
    s = res.scenario
    w = s.windows
    L = res.log
    out = {"cycle": s.name, "duration_s": s.duration, "max_rpm": float(L["rpm"].max()),
           "ecu_status": res.status}
    if "steady" in w:
        lam = res.event_lambdas(*w["steady"])
        out["steady_lambda_mean"] = float(lam.mean()) if lam.size else float("nan")
        out["steady_lambda_maxdev"] = float(np.abs(lam - 1.0).max()) if lam.size else float("nan")
    if "transient" in w:
        lam = res.event_lambdas(*w["transient"])
        out["transient_lambda_min"] = float(lam.min())
        out["transient_lambda_max"] = float(lam.max())
        out["transient_lambda_maxdev"] = float(np.abs(lam - 1.0).max())
        sens = res.window("lambda_sensor", *w["transient"])
        out["transient_sensor_maxdev"] = float(np.abs(sens - 1.0).max())
    if "after_settle" in w:
        r = res.window("rpm", *w["after_settle"])
        out["idle_rpm_mean_after"] = float(r.mean())
        out["idle_rpm_maxdev_after"] = float(np.abs(r - float(params.engine_spec()["engine"]["idle_rpm"])).max())
        r0 = res.window("rpm", *w["before_step"])
        out["idle_rpm_mean_before"] = float(r0.mean())
        r_dip = res.window("rpm", w["step_at"], w["after_settle"][0])
        out["idle_rpm_min_after_step"] = float(r_dip.min())
    if s.dyno_rpm is not None:
        out["final_ltft"] = float(L["ltft"][-1])
        out["final_stft"] = float(L["stft"][-1])
    return out
