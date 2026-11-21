"""Step-wise closed-loop HIL session (M9): the loop of hil.runner.run, one plant step per call, so
an outer controller (the gymnasium environment in ml/fuel) can act between steps.

`fuel_trim` scales the pulse width of every injection the ECU commits from now on:
pulse' = round(pulse x (1 + fuel_trim)). With fuel_trim = 0 the session reproduces
hil.runner.run exactly (same event order, same arithmetic; checked in tests/test_m9_hil.py).
All results are MODELLED."""
from __future__ import annotations

import heapq

from hil import runner
from hil.bridge import HostEcu
from mvem import params
from mvem.plant import DEG_PER_RAD, Engine


class Session:
    def __init__(self, scn: runner.Scenario, ecu=None):
        self.scn = scn
        self.ecu = ecu or HostEcu()
        self.ecu.reset()
        self.k = self.ecu.consts
        self.hz = float(self.k.timer_hz)
        self.dt = float(params.mvem_spec()["sim"]["dt"])
        self.eng = Engine(scn.init_rpm, scn.init_map)
        self.eng.knock_offset_deg = scn.knock_offset_deg
        self.pending: list = []
        self.seq = 0
        self.fuel_trim = 0.0
        self.out = None
        self.pedal = runner._hist(scn.pedal, 0.0)
        self.n_steps = int(round(scn.duration / self.dt))
        self.step_count = 0

    @property
    def done(self) -> bool:
        return self.step_count >= self.n_steps

    def _apply(self, te, typ, cyl, arg, th0, om_deg, t0):
        if typ == self.k.ev_inj:
            self.seq += 1
            heapq.heappush(self.pending, (te + arg / self.hz, self.seq, runner.INJ_END, cyl, arg))
        elif typ == runner.INJ_END:
            self.eng.inject(cyl, arg / self.hz)
        elif typ == self.k.ev_spark:
            self.eng.spark(cyl, th0 + om_deg * (te - t0))

    def step(self):
        """One plant step (spec/mvem.yaml sim.dt) = one ECU slow-loop iteration. Returns the ECU outputs."""
        scn, eng, ecu, hz, dt = self.scn, self.eng, self.ecu, self.hz, self.dt
        t = eng.t
        self.pedal = pedal = runner._hist(scn.pedal, t)
        load = runner._hist(scn.load, t)
        if scn.dyno_rpm is not None:
            dyno = runner._hist(scn.dyno_rpm, t)
        else:
            dyno = scn.init_rpm if t < scn.hold_s else None
        tps1, tps2 = eng.tps()
        out = ecu.slow(int(round(t * hz)), eng.p_man, eng.T_man, scn.clt, tps1, tps2, pedal,
                       eng.lam_sensor, eng.vbatt)
        volts = out.etc_voltage
        th0 = eng.theta
        om_deg = eng.omega * DEG_PER_RAD
        if om_deg > 0.0:
            for thc, kind, idx in eng.crossings(th0, th0 + om_deg * dt):
                tc = t + (thc - th0) / om_deg
                while self.pending and self.pending[0][0] <= tc:
                    te, _, typ, cyl, arg = heapq.heappop(self.pending)
                    self._apply(te, typ, cyl, arg, th0, om_deg, t)
                if kind == "tooth" or kind == "cam":
                    tick = int(round(tc * hz))
                    if kind == "tooth":
                        ecu.crank(tick)
                    else:
                        ecu.cam(tick)
                    for (tk, typ, cyl, arg) in ecu.drain():
                        if typ == self.k.ev_inj and self.fuel_trim != 0.0:
                            arg = max(int(round(arg * (1.0 + self.fuel_trim))), 0)
                        self.seq += 1
                        heapq.heappush(self.pending, (tk / hz, self.seq, typ, cyl, arg))
                elif kind == "ivc":
                    eng.ivc(idx, tc)
                else:
                    _, counts = eng.burn(idx, tc)
                    ecu.knock(idx, counts)
        while self.pending and self.pending[0][0] <= t + dt:
            te, _, typ, cyl, arg = heapq.heappop(self.pending)
            self._apply(te, typ, cyl, arg, th0, om_deg, t)
        eng.advance(dt, etc_volts=volts, load_torque=load, dyno_rpm=dyno)
        self.out = out
        self.step_count += 1
        return out
