"""M9.1 / M9.2 real-time co-simulation: `python -m hil.realtime` (also run by `make -C hil cosim`).

Plant time advances with core cycles (ADR 0005). The MVEM plant and the Vortex-V RTL share one
clock: t = cycles / core_clock_hz. The Verilator model plays the ECU's capture/compare timer
(hil/cosim/vortex_sim.cpp, ecu/targets/rv64/hil_mbox.h):
  * a crank / cam edge is stored in the capture FIFO at the cycle it occurs, with its captured tick;
    the firmware sees it when its interrupt handler (or, in the polled mode, its main loop) gets
    there. A full FIFO drops the edge: a missed tooth.
  * an event the firmware arms takes effect at max(its tick, the cycle the store reached RAM): an
    output compare armed late fires late, and the plant ignites / injects at that later time.
  * the slow loop runs on the firmware's own schedule; its outputs reach the plant when committed.
So handler latency and CPU overload appear as timing errors in the events, which the frozen-time
lockstep mode (hil/cosim.py, kept for the bit-exactness tests) cannot show.

M9.2 (ADR 0007, mode "capcmp"): the SoC is built with Vortex-V's capture/compare timer (its ADR
0017) and nothing is played any more. The plant drives the cap_in pins (crank = cap_in[0], cam =
cap_in[1]) at the cycle an edge occurs; the firmware's external interrupt reads the RTL FIFO; its
events are compare slots, and the plant acts on the cmp_out pin edges (injector c = cmp_out[c]:
open while high; coil c = cmp_out[4 + c]: dwell on the rising edge, spark on the falling edge).
Every pin edge is recorded with its cycle and timer tick and matched to the push that armed it.
Knock windows still go through the mailbox (the unit has no knock input).

Reported per run: worst tooth-to-event and tooth-to-serviced latency, spark angle error against the
host build of the same sources run in ideal (frozen) time, lateness of every event, CPU load, slow-
loop wall time and overruns, dropped captures and decoder sync counters. All values MODELLED."""
from __future__ import annotations

import argparse
import bisect
import ctypes as C
import heapq
import json
import struct
import sys
import time
from dataclasses import dataclass, field

import numpy as np

from hil import cosim, runner, rv64
from mvem import params
from mvem.plant import DEG_PER_RAD, Engine
from redline.spec import ROOT

OUT_PATH = ROOT / "out" / "hil" / "realtime.json"
SWEEP_PATH = ROOT / "out" / "hil" / "realtime_sweep.json"
CC_SWEEP_PATH = ROOT / "out" / "hil" / "realtime_capcmp_sweep.json"    # M9.2 (--cc-sweep)
US_PER_S = 1.0e6
DEG_PER_REV = 360.0
U32, U64 = 4, 8
RT_EVENT_RAW = 32                 # RT_EVENT_BYTES_MAX in vortex_sim.cpp
DRAIN_MAX = 4096

EV_TYPES = {"inj": 1, "dwell": 2, "spark": 3}            # ECU_EV_* (ecu_types.h; asserted in run_rt)

MODES = {                         # name -> (firmware build, HIL_RT_F_* flag names)
    "irq": ("fast", ("HIL_RT_F_IRQ", "HIL_RT_F_PREDICT")),          # M9.1 design
    "irq_periodic": ("fast", ("HIL_RT_F_IRQ",)),                    # fixed-period timer tick only
    "polled": ("fast", ()),                                         # M9 structure, M9.1 arithmetic
    "m9": ("m9", ()),                                               # M9 structure and code generation
    "m9_irq": ("m9", ("HIL_RT_F_IRQ", "HIL_RT_F_PREDICT")),         # M9.1 structure, M9 arithmetic
    "capcmp": ("fast", ("HIL_RT_F_IRQ", "HIL_RT_F_CAPCMP")),        # M9.2: SoC capture/compare timer
}
CC_PIN = {"HIL_CAP_CRANK": 0, "HIL_CAP_CAM": 1}         # cap_in bit per capture kind (ADR 0007)
CC_KIND = {1: "HIL_CAP_CRANK", 2: "HIL_CAP_CAM"}        # FIFO kind = channel + 1 (Vortex-V ADR 0017)
CC_ACT = {"pulse": 1, "set": 2, "clear": 3}             # CMP_PUSH action field
CC_PENDING_HIST = 256
CC_LATE_EDGE_CYCLES = 4           # a past-due compare's edge follows its push within this many cycles
CC_EDGE_AFTER_TICK = 3            # cmp_out changes this many cycles after TIME reaches the tick (with
                                  # the vx_cc_align phase: the timer increments P/2 - 2 cycles into a
                                  # plant tick, the registered output one cycle after the compare)
CC_SLOTS = 8                      # vortex_capcmp N_CMP (Vortex-V ADR 0017); hil_rv64.c CC_N_CMP
U32_MASK = 0xFFFFFFFF


class CcEdge(C.Structure):
    _fields_ = [("cycle", C.c_uint64), ("tick", C.c_uint32), ("before", C.c_uint32), ("after", C.c_uint32)]


class CcPush(C.Structure):
    _fields_ = [("cycle", C.c_uint64), ("wdata", C.c_uint64), ("arg", C.c_uint32), ("time", C.c_uint32),
                ("pending", C.c_uint32), ("pad", C.c_uint32)]


class CcCap(C.Structure):
    _fields_ = [("cycle", C.c_uint64), ("tick", C.c_uint32), ("kind", C.c_uint32)]


class RtEvent(C.Structure):
    _fields_ = [("arm_cycle", C.c_uint64), ("cap_rd", C.c_uint32), ("n_raw", C.c_uint32),
                ("raw", C.c_uint8 * RT_EVENT_RAW)]


@dataclass
class RtResult:
    window: str
    mode: str
    clk_hz: float
    rpm: float
    log: dict
    events: list            # dicts: tick, type, cyl, arg, arm_cycle, cap (index of the capture in service)
    caps: list              # dicts: cycle, kind, a0, a1, dropped
    service_cycles: list    # cycle at which capture i was fully serviced
    slow_cycles: list       # cycle of every slow-loop commit
    irq_cycles: list
    stats: dict
    status: dict
    combustion: list
    wall_s: float
    host: dict = field(default_factory=dict)
    cc: dict = field(default_factory=dict)      # M9.2: capture/compare observations


def _bind(lib):
    vp, u64 = C.c_void_p, C.c_uint64
    lib.vx_rt_begin.argtypes = [vp, C.POINTER(u64), C.c_uint32, C.c_uint32, u64, u64]
    lib.vx_rt_end.argtypes = [vp]
    lib.vx_rt_run_to.argtypes = [vp, u64]
    lib.vx_rt_events.argtypes = [vp, C.POINTER(RtEvent), u64]
    lib.vx_rt_events.restype = u64
    for f in (lib.vx_rt_services, lib.vx_rt_slows, lib.vx_rt_irqs):
        f.argtypes, f.restype = [vp, C.POINTER(u64), u64], u64
    lib.vx_cc_align.argtypes, lib.vx_cc_align.restype = [vp], C.c_uint32
    lib.vx_cc_end.argtypes = [vp]
    lib.vx_cc_pin.argtypes = [vp, C.c_uint32, C.c_uint32]
    for f in (lib.vx_cc_time, lib.vx_cc_pending_max, lib.vx_cc_cap_ovf, lib.vx_cc_cmp_drop):
        f.argtypes, f.restype = [vp], C.c_uint32
    lib.vx_cc_pending_hist.argtypes = [vp, C.POINTER(u64), u64]
    for f, t in ((lib.vx_cc_edges, CcEdge), (lib.vx_cc_pushes, CcPush), (lib.vx_cc_caps, CcCap), (lib.vx_cc_pops, CcCap)):
        f.argtypes, f.restype = [vp, C.POINTER(t), u64], u64


def _drain_all(fn, h, ctype):
    buf = (ctype * DRAIN_MAX)()
    out = []
    while True:
        m = fn(h, buf, DRAIN_MAX)
        out += [{k: getattr(buf[i], k) for k, _ in ctype._fields_} for i in range(m)]
        if m < DRAIN_MAX:
            return out


class RecordingEcu:
    """Wraps an ECU bridge and keeps every drained event (the host run in ideal time)."""

    def __init__(self, ecu):
        self.ecu = ecu
        self.consts = ecu.consts
        self.events: list = []

    def __getattr__(self, name):
        return getattr(self.ecu, name)

    def drain(self):
        ev = self.ecu.drain()
        self.events += ev
        return ev


def host_run(window: str) -> dict:
    """The host reference build on the same scenario in ideal time: its events and plant log."""
    rec = RecordingEcu(rv64.RefEcu())
    res = runner.run(cosim.scenario(window), rec)
    return {"events": rec.events, "log": res.log, "status": res.status, "combustion": res.events}


def run_rt(window: str, ecu: rv64.Rv64Ecu, mode: str = "irq", clk_hz: float | None = None) -> RtResult:
    """One window in real time on `ecu` (which must be the firmware build MODES[mode] names)."""
    cfg = rv64.cosim_spec()
    rtc = cfg["realtime"]
    clk = float(clk_hz or cfg["core_clock_hz"])
    E, lib, h = ecu.E, ecu.lib, ecu.h
    _bind(lib)
    scn = cosim.scenario(window)
    k = ecu.consts
    hz = float(k.timer_hz)
    assert (k.ev_inj, k.ev_dwell, k.ev_spark) == (EV_TYPES["inj"], EV_TYPES["dwell"], EV_TYPES["spark"])
    dt = float(params.mvem_spec()["sim"]["dt"])
    flags = 0
    for name in MODES[mode][1]:
        flags |= E[name]
    cc = "HIL_RT_F_CAPCMP" in MODES[mode][1]
    if cc and abs(clk / hz - round(clk / hz)) > 1e-9:
        raise ValueError(f"capture/compare mode needs a whole number of clocks per timer tick ({clk} Hz)")

    ecu.reset()
    off = ecu._l
    rt_cfg = [0] * E["HIL_RT_NCFG"]
    rt_cfg[E["HIL_RT_CLK_HZ"]] = int(round(clk))
    rt_cfg[E["HIL_RT_FLAGS"]] = flags
    rt_cfg[E["HIL_RT_POLL_CYCLES"]] = int(rtc["poll_cycles"])
    rt_cfg[E["HIL_RT_LEAD_CYCLES"]] = int(rtc["lead_cycles"])
    rt_cfg[E["HIL_RT_SPIN_CYCLES"]] = int(rtc["spin_cycles"])
    ecu._wr(off("HIL_L_RT_CFG"), struct.pack(f"<{len(rt_cfg)}I", *rt_cfg))
    ecu._wr(off("HIL_L_RT_STOP"), struct.pack("<I", 0))
    ecu._wr(off("HIL_L_CMD"), struct.pack("<I", E["HIL_CMD_RT"]))
    if lib.vx_run_until(h, ecu.mbox + off("HIL_L_RT_READY"), 1, ecu._max_cmd) < 0:
        raise RuntimeError("firmware did not enter real-time mode")

    eng = Engine(scn.init_rpm, scn.init_map)
    eng.knock_offset_deg = scn.knock_offset_deg
    in_size, out_size = off("HIL_L_IN_SIZE"), off("HIL_L_OUT_SIZE")
    cap_size, n_slots = off("HIL_L_RT_CAP_SIZE"), E["HIL_CAP_SLOTS"]
    ev_size = off("HIL_L_EVENT_SIZE")
    ot, oa, oy, oc = (off(n) for n in ("HIL_L_EVENT_TICK", "HIL_L_EVENT_ARG", "HIL_L_EVENT_TYPE", "HIL_L_EVENT_CYL"))
    inp = ecu.inp
    out = type(ecu.out)()
    in_sel = 0

    def post_inputs(t):
        nonlocal in_sel
        tps1, tps2 = eng.tps()
        inp.now_tick = int(round(t * hz)) & 0xFFFFFFFF       # the firmware reads its own timer instead
        inp.map_pa, inp.iat_k, inp.clt_k, inp.tps1_deg, inp.tps2_deg = eng.p_man, eng.T_man, scn.clt, tps1, tps2
        inp.pedal, inp.lambda_, inp.vbatt = runner._hist(scn.pedal, t), eng.lam_sensor, eng.vbatt
        in_sel ^= 1
        ecu._wr(off("HIL_L_RT_IN") + in_sel * in_size, bytes(inp))
        ecu._wr(off("HIL_L_RT_IN_SEL"), struct.pack("<I", in_sel))

    def u32(name):
        return struct.unpack("<I", ecu._rd(off(name), U32))[0]

    post_inputs(0.0)
    offs = (C.c_uint64 * 5)(*(ecu.mbox + off(n) for n in ("HIL_L_RT_TIMER", "HIL_L_RT_CAP_RD", "HIL_L_RT_EVT_WR",
                                                             "HIL_L_EVENTS", "HIL_L_RT_SLOW_SEQ")))
    lib.vx_rt_begin(h, offs, ev_size, E["HIL_MBOX_EVENTS"], int(round(clk)), int(k.timer_hz))
    t0_hw = 0
    if cc:                                           # plant time 0 on the hardware timer's phase
        t0_hw = int(lib.vx_cc_align(h))
        ecu._wr(off("HIL_L_RT_T0"), struct.pack("<I", t0_hw))
    ecu._wr(off("HIL_L_RT_GO"), struct.pack("<I", 1))

    events, caps, pending = [], [], []
    cc_edges, inj_open = [], {}
    seq = 0
    cap_wr = cap0 = u32("HIL_L_RT_CAP_WR")           # free-running across runs
    ev_buf = (RtEvent * DRAIN_MAX)()
    wall0 = time.time()

    def run_to(t):
        if lib.vx_rt_run_to(h, int(round(t * clk))):
            raise RuntimeError("firmware ended during the real-time run")
        nonlocal seq
        if cc:                                       # output pins: the plant acts on their edges
            for ed in _drain_all(lib.vx_cc_edges, h, CcEdge):
                cc_edges.append(ed)
                te = ed["cycle"] / clk
                for o in range(2 * k.ncyl):
                    bit = 1 << o
                    if (ed["before"] ^ ed["after"]) & bit:
                        seq += 1
                        rising = bool(ed["after"] & bit)
                        heapq.heappush(pending, (te, seq, ("cc", o, rising), o % k.ncyl, 0))
        for i in range(lib.vx_rt_events(h, ev_buf, DRAIN_MAX)):
            e = ev_buf[i]
            raw = bytes(e.raw[:ev_size])
            ev = {"tick": struct.unpack_from("<I", raw, ot)[0], "arg": struct.unpack_from("<I", raw, oa)[0],
                  "type": raw[oy], "cyl": raw[oc], "arm_cycle": int(e.arm_cycle),
                  "cap": (int(e.cap_rd) - cap0) & 0xFFFFFFFF}
            # output compare: at its tick, or at once if that has passed when it is armed
            ev["t_fire"] = max(ev["tick"] / hz, ev["arm_cycle"] / clk)
            events.append(ev)
            seq += 1
            heapq.heappush(pending, (ev["t_fire"], seq, ev["type"], ev["cyl"], ev["arg"]))

    def post(t, kind, a0, a1=0):
        nonlocal cap_wr
        rec = {"cycle": int(round(t * clk)), "kind": kind, "a0": a0, "a1": a1, "dropped": False}
        if cc and kind in CC_PIN:                      # a sensor edge on the capture pin
            lib.vx_cc_pin(h, CC_PIN[kind], int(cfg["capcmp"]["pin_width_cycles"]))
            caps.append(rec)
            return
        if cap_wr - u32("HIL_L_RT_CAP_RD") >= n_slots:
            rec["dropped"] = True                      # FIFO full: the edge is lost
        else:
            tag = (t0_hw + int(round(t * hz))) & U32_MASK if cc else 0     # knock: tick, for the merge
            ecu._wr(off("HIL_L_RT_CAP") + (cap_wr % n_slots) * cap_size,
                    struct.pack("<4I", E[kind], a0 & 0xFFFFFFFF, a1 & 0xFFFFFFFF, tag))
            cap_wr += 1
            ecu._wr(off("HIL_L_RT_CAP_WR"), struct.pack("<I", cap_wr & 0xFFFFFFFF))
        caps.append(rec)

    def inj_credit(cyl, until):
        """Fuel of an open injector up to `until` (pin semantics: an injector held open across two
        pulses, or across its intake valve closing, is credited in pieces; its dead time once)."""
        t_open, first = inj_open[cyl]
        eng.inject(cyl, until - t_open + (0.0 if first else eng.deadtime()))
        inj_open[cyl] = [until, False]

    def apply(te, typ, cyl, arg, th0, om_deg, t0):
        nonlocal seq
        if isinstance(typ, tuple):                     # M9.2: a cmp_out pin edge
            _, o, rising = typ
            if o < k.ncyl:                             # injector: fuel delivered when it closes
                if rising:
                    inj_open[cyl] = [te, True]
                elif cyl in inj_open:
                    inj_credit(cyl, te)
                    del inj_open[cyl]
            elif not rising:                           # coil current cut: the spark
                eng.spark(cyl, th0 + om_deg * (te - t0))
            return
        if typ == k.ev_inj:
            seq += 1
            heapq.heappush(pending, (te + arg / hz, seq, runner.INJ_END, cyl, arg))
        elif typ == runner.INJ_END:
            eng.inject(cyl, arg / hz)
        elif typ == k.ev_spark:
            eng.spark(cyl, th0 + om_deg * (te - t0))

    n = int(round(scn.duration / dt))
    log = {key: np.zeros(n) for key in runner.LOG_KEYS}
    for step in range(n):
        t = eng.t
        pedal = runner._hist(scn.pedal, t)
        load = runner._hist(scn.load, t)
        dyno = runner._hist(scn.dyno_rpm, t) if scn.dyno_rpm is not None else (scn.init_rpm if t < scn.hold_s else None)
        run_to(t)
        if step:
            post_inputs(t)
        th0 = eng.theta
        om_deg = eng.omega * DEG_PER_RAD
        if om_deg > 0.0:
            for thc, kind, idx in eng.crossings(th0, th0 + om_deg * dt):
                tc = t + (thc - th0) / om_deg
                run_to(tc)
                while pending and pending[0][0] <= tc:
                    te, _, typ, cyl, arg = heapq.heappop(pending)
                    apply(te, typ, cyl, arg, th0, om_deg, t)
                if kind == "tooth":
                    post(tc, "HIL_CAP_CRANK", int(round(tc * hz)))
                elif kind == "cam":
                    post(tc, "HIL_CAP_CAM", int(round(tc * hz)))
                elif kind == "ivc":
                    if cc and idx in inj_open:         # still open at IVC: what flowed so far is inducted
                        inj_credit(idx, tc)
                    eng.ivc(idx, tc)
                else:
                    _, counts = eng.burn(idx, tc)
                    post(tc, "HIL_CAP_KNOCK", idx, int(counts))
        run_to(t + dt)
        while pending and pending[0][0] <= t + dt:
            te, _, typ, cyl, arg = heapq.heappop(pending)
            apply(te, typ, cyl, arg, th0, om_deg, t)
        # actuator command: the last slow-loop outputs the firmware has committed by the end of the step
        sel = u32("HIL_L_RT_OUT_SEL") & 1
        C.memmove(C.byref(out), ecu._rd(off("HIL_L_RT_OUT") + sel * out_size, out_size), out_size)
        eng.advance(dt, etc_volts=out.etc_voltage, load_torque=load, dyno_rpm=dyno)
        row = (t, eng.rpm(), out.rpm, eng.p_man, eng.thr, out.throttle_target_deg, out.etc_voltage, eng.lam_sensor,
               pedal, out.spark_deg[0], out.stft, out.ltft, max(out.knock_retard_deg), eng.T_brake, eng.T_ind,
               out.closed_loop, out.idle_active, out.rev_cut, out.limp, out.sync, out.pw_s[0], load)
        for key, v in zip(runner.LOG_KEYS, row):
            log[key][step] = v
    wall = time.time() - wall0

    ecu._wr(off("HIL_L_RT_STOP"), struct.pack("<I", 1))
    if lib.vx_run_until(h, ecu.mbox + off("HIL_L_CMD"), E["HIL_CMD_NONE"], ecu._max_cmd) < 0:
        raise RuntimeError("firmware did not leave real-time mode")
    run_to(lib.vx_rt_now(h) / clk)                       # collect what was armed while stopping

    def stamps(fn):
        buf = (C.c_uint64 * DRAIN_MAX)()
        out_l = []
        while True:
            m = fn(h, buf, DRAIN_MAX)
            out_l += [int(buf[i]) for i in range(m)]
            if m < DRAIN_MAX:
                return out_l
    service, slows, irqs = stamps(lib.vx_rt_services), stamps(lib.vx_rt_slows), stamps(lib.vx_rt_irqs)
    cc_obs = {}
    if cc:
        hist = (C.c_uint64 * CC_PENDING_HIST)()
        lib.vx_cc_pending_hist(h, hist, CC_PENDING_HIST)
        cc_obs = {"t0": t0_hw, "edges": cc_edges, "pushes": _drain_all(lib.vx_cc_pushes, h, CcPush),
                  "caps": _drain_all(lib.vx_cc_caps, h, CcCap), "pops": _drain_all(lib.vx_cc_pops, h, CcCap),
                  "pending_max": int(lib.vx_cc_pending_max(h)), "pending_hist": [int(x) for x in hist],
                  "cap_ovf": int(lib.vx_cc_cap_ovf(h)), "cmp_drop": int(lib.vx_cc_cmp_drop(h)),
                  "mbox_service": service, "end_cycle": int(lib.vx_rt_now(h))}
        events, service = _cc_events(cc_obs, caps, clk, hz, k.ncyl)
        lib.vx_cc_end(h)
    lib.vx_rt_end(h)
    n_stats = E["HIL_RT_NSTATS"]
    vals = struct.unpack(f"<{n_stats}Q", ecu._rd(off("HIL_L_RT_STATS"), n_stats * U64))
    names = {v: key[len("HIL_RT_S_"):].lower() for key, v in E.items() if key.startswith("HIL_RT_S_")}
    stats = {names[i]: int(vals[i]) for i in range(n_stats)}
    return RtResult(window=window, mode=mode, clk_hz=clk, rpm=float(log["rpm"][-1]), log=log, events=events,
                    caps=caps, service_cycles=service, slow_cycles=slows, irq_cycles=irqs, stats=stats,
                    status=ecu.status(), combustion=eng.events, wall_s=wall, cc=cc_obs)


def _cc_events(obs: dict, caps: list, clk: float, hz: float, ncyl: int) -> tuple:
    """M9.2: turn the observed compare pushes and cmp_out edges into the event list metrics() reads.

    Each CMP_PUSH is one event (type and cylinder from its output and action, tick in the host's
    time base = pushed tick - T0). Its fire time is the cmp_out edge it caused. The unit fires a
    compare in the first clock where TIME >= tick (the next clock if that has passed at the push),
    so an on-time event is matched to the edge of its polarity on its output whose TIME equals its
    tick, a late one to the first such edge after the push. A compare that does not change the
    level (a set on a coil that is already charging, a pulse while the injector is still open, a
    clear on a coil that never charged) leaves no edge: it is kept with the unit's fire time and
    flagged `no_edge`, and `edge_missing` counts those whose level should have changed (a fault).
    `cap` names the tooth whose window the event is in (the last crank edge whose captured tick is
    not after the event's tick): the scheduler armed it from that tooth. Also returns the service
    stamps aligned with the posted captures: hardware pops for crank / cam, mailbox pops for knock."""
    t0 = obs["t0"]
    P = int(round(clk / hz))
    edges: dict = {}                                 # (out, rising) -> [edge]
    for ed in obs["edges"]:
        for o in range(2 * ncyl):
            bit = 1 << o
            if (ed["before"] ^ ed["after"]) & bit:
                edges.setdefault((o, bool(ed["after"] & bit)), []).append(ed)
    by_tick = {key: {e["tick"]: e for e in lst} for key, lst in edges.items()}
    cyc = {key: [e["cycle"] for e in lst] for key, lst in edges.items()}

    def level_before(o, cycle):                      # cmp_out[o] just before `cycle`
        lv = 0
        for ed in obs["edges"]:
            if ed["cycle"] >= cycle:
                break
            lv = (ed["after"] >> o) & 1
        return lv

    posted = [c for c in caps if not c["dropped"]]
    crank_idx = [i for i, c in enumerate(posted) if c["kind"] == "HIL_CAP_CRANK"]
    crank_tick = [posted[i]["a0"] for i in crank_idx]
    events = []
    for p in obs["pushes"]:
        w = p["wdata"]
        tick, o, act = w & U32_MASK, (w >> 32) & 0xFF, (w >> 40) & 0x7
        if o < ncyl and act == CC_ACT["pulse"]:
            typ, cyl, rising = EV_TYPES["inj"], o, True
        elif o >= ncyl and act == CC_ACT["set"]:
            typ, cyl, rising = EV_TYPES["dwell"], o - ncyl, True
        elif o >= ncyl and act == CC_ACT["clear"]:
            typ, cyl, rising = EV_TYPES["spark"], o - ncyl, False
        else:
            raise RuntimeError(f"unexpected compare push {w:#x}")
        late = (p["time"] - tick) & U32_MASK < 0x80000000
        key = (o, rising)
        if late:                                     # fires in the clock after the push
            lst = cyc.get(key, [])
            i = bisect.bisect_left(lst, p["cycle"])
            ed = edges[key][i] if i < len(lst) and lst[i] <= p["cycle"] + CC_LATE_EDGE_CYCLES else None
            rule = p["cycle"] + 1
        else:
            ed = by_tick.get(key, {}).get(tick)
            rule = (((tick - t0) & U32_MASK) * P) - P // 2 + CC_EDGE_AFTER_TICK
        j = bisect.bisect_right(crank_tick, (tick - t0) & U32_MASK) - 1   # the tooth whose window it is in
        ev = {"tick": (tick - t0) & U32_MASK, "arg": p["arg"] if typ == EV_TYPES["inj"] else 0, "type": typ,
              "cyl": cyl, "out": o, "arm_cycle": int(p["cycle"]), "cap": crank_idx[j] if j >= 0 else len(posted),
              "pushed_late": late, "pending_before": p["pending"],
              "edge_cycle": ed["cycle"] if ed else None, "edge_tick": ed["tick"] if ed else None,
              "no_edge": ed is None, "edge_missing": False}
        if ed is None:
            ev["edge_missing"] = level_before(o, rule) != (1 if rising else 0) and rule < obs["end_cycle"]
        ev["t_fire"] = (ed["cycle"] if ed else rule) / clk
        if typ == EV_TYPES["inj"] and ed:            # pulse length = the falling edge that follows
            fl = cyc.get((o, False), [])
            i = bisect.bisect_right(fl, ed["cycle"])
            ev["pulse_ticks"] = ((edges[(o, False)][i]["tick"] - ed["tick"]) & U32_MASK) if i < len(fl) else None
        events.append(ev)
    # a pulse is cut short or stretched only by another pulse on the same injector meanwhile
    for e in events:
        if e.get("pulse_ticks") is not None:
            nxt = [f for f in events if f["out"] == e["out"] and f["type"] == e["type"] and
                   0 < (f["tick"] - e["tick"]) & U32_MASK < e["pulse_ticks"] + 1 and f is not e]
            e["pulse_overlap"] = bool(nxt)
    hw = iter(obs["pops"])
    mb = iter(obs["mbox_service"])
    service = []
    for c in posted:
        nxt = next(hw if c["kind"] in CC_PIN else mb, None)
        if nxt is None:
            break
        service.append(nxt["cycle"] if isinstance(nxt, dict) else nxt)
    return events, service


# --------------------------------------------------------------------------- metrics
def _by_cyl(events, typ, key):
    out: dict = {}
    for e in events:
        if e["type"] == typ:
            out.setdefault(e["cyl"], []).append(e[key])
    return out


def metrics(r: RtResult, host: dict) -> dict:
    """Timing figures of one real-time run against the host run in ideal time (same scenario)."""
    cfg = rv64.cosim_spec()
    consts_hz = float(params.ecu_spec()["hardware"]["timer_hz"])
    types = EV_TYPES
    clk = r.clk_hz
    rpm = np.asarray(r.log["rpm"])
    dt = float(params.mvem_spec()["sim"]["dt"])

    def omega_deg(t):                                      # crank speed [deg/s] at time t
        i = min(max(int(t / dt), 0), len(rpm) - 1)
        return rpm[i] * DEG_PER_REV / params.SEC_PER_MIN

    posted = [c for c in r.caps if not c["dropped"]]
    # tooth -> event armed: the capture in service when the event was armed
    lat_event = [e["arm_cycle"] - posted[e["cap"]]["cycle"] for e in r.events if e["cap"] < len(posted)]
    lat_service = [s - c["cycle"] for s, c in zip(r.service_cycles, posted)]
    crank_service = [s - c["cycle"] for s, c in zip(r.service_cycles, posted) if c["kind"] == "HIL_CAP_CRANK"]
    late = [(e["t_fire"] - e["tick"] / consts_hz) for e in r.events]
    late_deg = [l * omega_deg(e["t_fire"]) for l, e in zip(late, r.events)]

    host_ev = [{"tick": tk, "type": ty, "cyl": cy, "arg": ar, "t_fire": tk / consts_hz} for tk, ty, cy, ar in host["events"]]
    out = {"window": r.window, "mode": r.mode, "clk_hz": clk, "rpm": r.rpm,
           "events": len(r.events), "host_events": len(host_ev)}
    # spark angle error vs host: same cylinder, same occurrence; speed imposed, so angle = speed x time
    err_deg, n_match, missing = [], 0, 0
    for name in ("spark", "dwell"):
        a, b = _by_cyl(r.events, types[name], "t_fire"), _by_cyl(host_ev, types[name], "t_fire")
        errs = []
        for cyl in set(a) | set(b):
            ta, tb = a.get(cyl, []), b.get(cyl, [])
            missing += abs(len(ta) - len(tb)) if name == "spark" else 0
            errs += [abs(x - y) * omega_deg(y) for x, y in zip(ta, tb)]
        out[f"{name}_angle_error_max_deg"] = float(max(errs)) if errs else float("nan")
        out[f"{name}_angle_error_mean_deg"] = float(np.mean(errs)) if errs else float("nan")
        if name == "spark":
            err_deg, n_match = errs, len(errs)
    pa, pb = _by_cyl(r.events, types["inj"], "arg"), _by_cyl(host_ev, types["inj"], "arg")
    pw = [abs(x - y) for cyl in set(pa) & set(pb) for x, y in zip(pa[cyl], pb[cyl])]
    s = r.stats
    spec_deg = float(cfg["realtime"]["angle_error_max_deg"])
    om_end = omega_deg(len(rpm) * dt)
    # every event's error: sparks and dwell against the host's fire time, and the lateness of any
    # event (an injection's start moves with its pulse width, so only its lateness is an error)
    observed = max([out["spark_angle_error_max_deg"], out["dwell_angle_error_max_deg"]] + late_deg, default=0.0)
    # an event can fall anywhere in a tooth window, also right at its start: it is then late by the
    # whole tooth-to-armed latency. The worst latency seen bounds that case.
    bound = (max(lat_event) / clk * om_end) if lat_event else 0.0
    idle_cycles = s["idle_iters"] * s["idle_cal_cycles"] / max(s["idle_cal_iters"], 1)
    out.update({
        "sparks_compared": n_match, "spark_count_difference": missing,
        "event_angle_error_max_deg": float(observed),
        "worst_case_angle_error_deg": float(bound),
        "angle_within_spec": bool(err_deg and missing == 0 and max(observed, bound) <= spec_deg),
        "inj_pulse_diff_max_us": float(max(pw)) / consts_hz * US_PER_S if pw else float("nan"),
        "event_late_max_us": float(max(late)) * US_PER_S if late else 0.0,
        "event_late_max_deg": float(max(late_deg)) if late_deg else 0.0,
        "events_late": int(sum(l > 0.0 for l in late)),
        "tooth_to_event_max_cycles": int(max(lat_event)) if lat_event else 0,
        "tooth_to_event_max_us": float(max(lat_event)) / clk * US_PER_S if lat_event else 0.0,
        "tooth_to_serviced_max_cycles": int(max(crank_service)) if crank_service else 0,
        "tooth_to_serviced_max_us": float(max(crank_service)) / clk * US_PER_S if crank_service else 0.0,
        "tooth_to_serviced_mean_us": float(np.mean(crank_service)) / clk * US_PER_S if crank_service else 0.0,
        "capture_to_serviced_max_us": float(max(lat_service)) / clk * US_PER_S if lat_service else 0.0,
        "tooth_period_us": cosim.tooth_period_s(r.rpm) * US_PER_S,
        "cpu_load": 1.0 - idle_cycles / max(s["run_cycles"], 1),
        "slow_steps": s["slow_steps"], "slow_steps_expected": len(rpm),
        "slow_overruns": s["slow_overruns"],
        "slow_wall_max_us": s["slow_wall_max"] / clk * US_PER_S,
        "slow_wall_mean_us": s["slow_wall_sum"] / max(s["slow_steps"], 1) / clk * US_PER_S,
        "irqs": s["irqs"], "isr_max_cycles": s["isr_max"],
        "isr_mean_cycles": s["isr_cycles"] / max(s["irqs"], 1), "spin_timeouts": s["spin_timeouts"],
        "snapshot_retries": s["snap_retries"],
        "captures_posted": len(posted), "captures_dropped": len(r.caps) - len(posted),
        "captures_unserviced": len(posted) - len(r.service_cycles),
        "sync": r.status["sync"], "n_sync": r.status["n_sync"], "n_loss": r.status["n_loss"],
        "n_dropout": r.status["n_dropout"], "n_noise": r.status["n_noise"], "n_cam_err": r.status["n_cam_err"],
        "host_n_sync": host["status"]["n_sync"],
        "wall_s": r.wall_s,
    })
    if r.cc:
        out.update(_cc_metrics(r, host, omega_deg))
    out["teeth_ok"] = bool(out["captures_dropped"] == 0 and out["n_loss"] == 0 and out["n_dropout"] == 0 and
                           out["n_noise"] == 0 and out["sync"] == 2 and out["n_sync"] == host["status"]["n_sync"])
    out["rt_ok"] = bool(out["teeth_ok"] and out["slow_overruns"] == 0 and out["captures_unserviced"] <= 1)
    if r.cc:                                       # the unit must not have lost a capture or a compare
        out["rt_ok"] = bool(out["rt_ok"] and out["slots_ok"] and out["cap_ovf"] == 0)
    out["meets_spec"] = bool(out["angle_within_spec"] and out["rt_ok"])
    return out


def _cc_metrics(r: RtResult, host: dict, omega_deg) -> dict:
    """M9.2 figures from the capture/compare observations (pins, pushes, FIFO)."""
    o = r.cc
    clk = r.clk_hz
    hz = float(params.ecu_spec()["hardware"]["timer_hz"])
    s = r.stats
    pins = [c for c in r.caps if c["kind"] in CC_PIN]
    # every sensor edge on the pins must be in the FIFO, with the tick the host build was given
    got = o["caps"]
    tick_bad = sum(1 for c, g in zip(pins, got)
                   if CC_KIND.get(g["kind"]) != c["kind"] or g["tick"] != (o["t0"] + c["a0"]) & U32_MASK)
    lost = len(pins) - len(got)
    # tooth -> interrupt taken (first dbg_irq stamp at or after the edge), crank edges only
    irqs = np.asarray(r.irq_cycles, dtype=np.int64)
    to_irq = []
    for c in pins:
        if c["kind"] == "HIL_CAP_CRANK" and len(irqs):
            i = int(np.searchsorted(irqs, c["cycle"]))
            if i < len(irqs):
                to_irq.append(int(irqs[i]) - c["cycle"])
    ev = r.events
    no_edge = sum(e["no_edge"] for e in ev)
    missing = sum(e["edge_missing"] for e in ev)
    pulses = [e for e in ev if e.get("pulse_ticks") is not None and not e.get("pulse_overlap")]
    pulse_bad = sum(1 for e in pulses if e["pulse_ticks"] != e["arg"])
    late = [e for e in ev if e["pushed_late"]]
    late_deg = [(e["t_fire"] - e["tick"] / hz) * omega_deg(e["t_fire"]) for e in late]
    # events pushed with exactly the host's (tick, type, cyl, arg), in order per cylinder and type
    host_by = {}
    for tk, ty, cy, ar in host["events"]:
        host_by.setdefault((ty, cy), []).append((tk, ar))
    mine_by = {}
    for e in ev:
        mine_by.setdefault((e["type"], e["cyl"]), []).append((e["tick"], e["arg"]))
    same = sum(sum(a == b for a, b in zip(mine_by.get(key, []), host_by.get(key, []))) for key in set(mine_by) | set(host_by))
    pend_max = int(o["pending_max"])
    return {
        "cc_t0": o["t0"], "cc_pushes": len(o["pushes"]), "cc_pushes_fw": s.get("cc_pushes", 0),
        "cc_edges": len(o["edges"]), "cc_no_edge": no_edge, "cc_edge_missing": missing,
        "cc_pulses_checked": len(pulses), "cc_pulse_mismatches": pulse_bad,
        "cc_pushed_late": len(late), "cc_pushed_late_max_deg": float(max(late_deg)) if late_deg else 0.0,
        "cc_events_identical_to_host": same, "cc_capture_tick_mismatches": tick_bad, "cc_captures_lost": lost,
        "cap_ovf": max(o["cap_ovf"], s.get("cc_cap_ovf", 0)),
        "cc_pending_max": pend_max, "cc_slots": CC_SLOTS, "cc_pending_hist": o["pending_hist"][:CC_SLOTS + 1],
        "cc_cmp_drop": max(o["cmp_drop"], s.get("cc_drops", 0)), "cc_deferred": s.get("cc_deferred", 0),
        "cc_spurious_irqs": s.get("cc_spurious", 0),
        "slots_ok": bool(pend_max <= CC_SLOTS and o["cmp_drop"] == 0 and s.get("cc_drops", 0) == 0 and
                         s.get("cc_deferred", 0) == 0),
        "tooth_to_irq_max_cycles": int(max(to_irq)) if to_irq else 0,
        "tooth_to_irq_mean_cycles": float(np.mean(to_irq)) if to_irq else 0.0,
        "captures_dropped": lost + o["cap_ovf"] + sum(c["dropped"] for c in r.caps),   # pins + mailbox knock
    }


class Bench:
    """Firmware builds on the SoC model (one model instance per build) and cached host runs."""

    def __init__(self):
        self._ecus: dict = {}
        self._host: dict = {}

    def ecu(self, build: str) -> rv64.Rv64Ecu:
        if build not in self._ecus:
            if build == "m9":
                rv64.build("m9")
            self._ecus[build] = rv64.Rv64Ecu(elf=rv64.M9_ELF if build == "m9" else rv64.ELF)
        return self._ecus[build]

    def host(self, window: str) -> dict:
        if window not in self._host:
            self._host[window] = host_run(window)
        return self._host[window]

    def run(self, window: str, mode: str = "irq", clk_hz: float | None = None) -> tuple:
        r = run_rt(window, self.ecu(MODES[mode][0]), mode, clk_hz)
        return r, metrics(r, self.host(window))

    def min_clock(self, window: str, mode: str = "irq", section: str = "realtime") -> dict:
        """Lowest core clock that meets the spec on `window`, by bisection (hil/cosim.yaml
        `<section>.min_clock`; M9.2 uses the capcmp section's whole-MHz grid)."""
        mc = rv64.cosim_spec()[section]["min_clock"]
        lo, hi, res = float(mc["lo_hz"]), float(mc["hi_hz"]), float(mc["resolution_hz"])
        tried = {}

        def ok(f):
            tried[f] = self.run(window, mode, f)[1]
            return tried[f]["meets_spec"]
        if not ok(hi):
            return {"min_clock_hz": None, "tried": tried}
        if ok(lo):
            return {"min_clock_hz": lo, "tried": tried}
        while hi - lo > res:
            mid = round((lo + hi) / 2.0 / res) * res
            if mid in (lo, hi):
                break
            if ok(mid):
                hi = mid
            else:
                lo = mid
        return {"min_clock_hz": hi, "tried": tried}


def line(m: dict) -> str:
    cc = ""
    if "cc_pending_max" in m:
        cc = (f" | capcmp: pending max {m['cc_pending_max']}/{m['cc_slots']}, pushes late {m['cc_pushed_late']}, "
              f"deferred {m['cc_deferred']}, dropped {m['cc_cmp_drop']}, CAP_OVF {m['cap_ovf']}, "
              f"tooth->irq {m['tooth_to_irq_max_cycles']} cyc")
    return (f"{m['window']:10s} {m['mode']:12s} {m['clk_hz'] / 1e6:5.1f} MHz {m['rpm']:5.0f} rpm | "
            f"tooth->event {m['tooth_to_event_max_us']:6.1f} us, ->serviced {m['tooth_to_serviced_max_us']:6.1f} us "
            f"({m['tooth_period_us']:.0f} us/tooth) | angle err: spark {m['spark_angle_error_max_deg']:.3f}, any event "
            f"{m['event_angle_error_max_deg']:.3f}, worst case {m['worst_case_angle_error_deg']:.3f} deg | load {100 * m['cpu_load']:5.1f} % | slow {m['slow_wall_max_us']:5.0f} us "
            f"max, {m['slow_steps']}/{m['slow_steps_expected']} steps, {m['slow_overruns']} overruns | dropped "
            f"{m['captures_dropped']}, loss {m['n_loss']}, sync {m['sync']} | real time {'HELD' if m['rt_ok'] else 'LOST'}, "
            f"angle spec {'MET' if m['meets_spec'] else 'NOT MET'}" + cc)


def main(argv=None) -> int:
    rtc = rv64.cosim_spec()["realtime"]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--windows", nargs="+", default=list(rtc["windows"]))
    ap.add_argument("--modes", nargs="+", default=["m9", "irq", "capcmp"], choices=list(MODES))
    ap.add_argument("--clocks", nargs="+", type=float, default=None, help="core clocks [Hz] (default: the board clock)")
    ap.add_argument("--sweep", action="store_true", help="clock sweep and minimum-clock bisection (hil/cosim.yaml)")
    ap.add_argument("--cc-sweep", action="store_true",
                    help="M9.2: capture/compare mode over capcmp.clocks_hz and its minimum-clock bisection")
    a = ap.parse_args(argv)
    if a.cc_sweep:
        return cc_sweep()
    missing = rv64.toolchain_missing()
    if missing:
        print(f"missing tools: {missing}", file=sys.stderr)
        return 2
    bench = Bench()
    base = float(rv64.cosim_spec()["core_clock_hz"])
    clocks = a.clocks or ([float(c) for c in rtc["clocks_hz"]] if a.sweep else [base])
    print("Redline M9.1 / M9.2 real-time co-simulation: plant time = core cycles. All values MODELLED.")
    print(f"core: vortex_soc, {rv64.VORTEX_V_ROOT}; spec: spark angle error vs host <= "
          f"{rtc['angle_error_max_deg']} deg for every event and for the worst tooth-to-event latency, no missed teeth, no loss of sync, no slow-loop overrun")
    rows = []
    for mode in a.modes:
        for w in a.windows:
            for f in clocks:
                m = bench.run(w, mode, f)[1]
                rows.append(m)
                print(line(m), flush=True)
    result = {"runs": rows}
    if a.sweep:
        mc = bench.min_clock(rtc["windows"][0], "irq")
        result["min_clock"] = {"window": rtc["windows"][0], "min_clock_hz": mc["min_clock_hz"],
                               "tried": {str(f): m for f, m in mc["tried"].items()}}
        for f, m in sorted(mc["tried"].items()):
            print(line(m))
        got = mc["min_clock_hz"]
        print(f"minimum core clock meeting the spec on {rtc['windows'][0]}: "
              + (f"{got / 1e6:.1f} MHz" if got else "none in the searched range"))
    out_path = SWEEP_PATH if a.sweep else OUT_PATH
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"wrote {out_path}")
    ok = True
    for mode, label in (("irq", "M9.1 firmware, mailbox capture/compare"), ("capcmp", "M9.2, SoC capture/compare timer")):
        mine = [m for m in rows if m["mode"] == mode and m["clk_hz"] == base]
        if not mine:
            continue
        held = all(m["rt_ok"] for m in mine)
        spec = all(m["meets_spec"] for m in mine)
        ok = ok and held
        print(f"RESULT at {base / 1e6:.0f} MHz ({label}): real time "
              + ("HELD (no missed tooth, no loss of sync, no slow-loop overrun"
                 + (", no compare slot overflow, no CAP_OVF)" if mode == "capcmp" else ")") if held else "LOST")
              + f"; {rtc['angle_error_max_deg']} deg angle spec " + ("MET" if spec else "NOT MET on every window (see the table)"))
    return 0 if ok else 1


def cc_sweep() -> int:
    """M9.2 clock sweep: capture/compare mode on every capcmp window at capcmp.clocks_hz, then the
    lowest whole-MHz clock that meets the spec on the first window."""
    ccc = rv64.cosim_spec()["capcmp"]
    bench = Bench()
    print("Redline M9.2 clock sweep, SoC capture/compare timer (ADR 0007). All values MODELLED.")
    rows = []
    for w in ccc["windows"]:
        for f in ccc["clocks_hz"]:
            m = bench.run(w, "capcmp", float(f))[1]
            rows.append(m)
            print(line(m), flush=True)
    win = ccc["windows"][0]
    mc = bench.min_clock(win, "capcmp", "capcmp")
    for f, m in sorted(mc["tried"].items()):
        print(line(m))
    got = mc["min_clock_hz"]
    print(f"minimum core clock meeting the spec on {win} (capcmp): "
          + (f"{got / 1e6:.0f} MHz" if got else "none in the searched range"))
    result = {"runs": rows, "min_clock": {"window": win, "min_clock_hz": got,
                                          "tried": {str(f): m for f, m in mc["tried"].items()}}}
    CC_SWEEP_PATH.parent.mkdir(parents=True, exist_ok=True)
    CC_SWEEP_PATH.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"wrote {CC_SWEEP_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
