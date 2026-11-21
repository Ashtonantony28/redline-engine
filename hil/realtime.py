"""M9.1 real-time co-simulation: `python -m hil.realtime` (also run by `make -C hil cosim`).

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

Reported per run: worst tooth-to-event and tooth-to-serviced latency, spark angle error against the
host build of the same sources run in ideal (frozen) time, lateness of every event, CPU load, slow-
loop wall time and overruns, dropped captures and decoder sync counters. All values MODELLED."""
from __future__ import annotations

import argparse
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
}


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


def _bind(lib):
    vp, u64 = C.c_void_p, C.c_uint64
    lib.vx_rt_begin.argtypes = [vp, C.POINTER(u64), C.c_uint32, C.c_uint32, u64, u64]
    lib.vx_rt_end.argtypes = [vp]
    lib.vx_rt_run_to.argtypes = [vp, u64]
    lib.vx_rt_events.argtypes = [vp, C.POINTER(RtEvent), u64]
    lib.vx_rt_events.restype = u64
    for f in (lib.vx_rt_services, lib.vx_rt_slows, lib.vx_rt_irqs):
        f.argtypes, f.restype = [vp, C.POINTER(u64), u64], u64


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
    ecu._wr(off("HIL_L_RT_GO"), struct.pack("<I", 1))

    events, caps, pending = [], [], []
    seq = 0
    cap_wr = cap0 = u32("HIL_L_RT_CAP_WR")           # free-running across runs
    ev_buf = (RtEvent * DRAIN_MAX)()
    wall0 = time.time()

    def run_to(t):
        if lib.vx_rt_run_to(h, int(round(t * clk))):
            raise RuntimeError("firmware ended during the real-time run")
        nonlocal seq
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
        if cap_wr - u32("HIL_L_RT_CAP_RD") >= n_slots:
            rec["dropped"] = True                      # FIFO full: the edge is lost
        else:
            ecu._wr(off("HIL_L_RT_CAP") + (cap_wr % n_slots) * cap_size,
                    struct.pack("<4I", E[kind], a0 & 0xFFFFFFFF, a1 & 0xFFFFFFFF, 0))
            cap_wr += 1
            ecu._wr(off("HIL_L_RT_CAP_WR"), struct.pack("<I", cap_wr & 0xFFFFFFFF))
        caps.append(rec)

    def apply(te, typ, cyl, arg, th0, om_deg, t0):
        nonlocal seq
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
    lib.vx_rt_end(h)
    n_stats = E["HIL_RT_NSTATS"]
    vals = struct.unpack(f"<{n_stats}Q", ecu._rd(off("HIL_L_RT_STATS"), n_stats * U64))
    names = {v: key[len("HIL_RT_S_"):].lower() for key, v in E.items() if key.startswith("HIL_RT_S_")}
    stats = {names[i]: int(vals[i]) for i in range(n_stats)}
    return RtResult(window=window, mode=mode, clk_hz=clk, rpm=float(log["rpm"][-1]), log=log, events=events,
                    caps=caps, service_cycles=service, slow_cycles=slows, irq_cycles=irqs, stats=stats,
                    status=ecu.status(), combustion=eng.events, wall_s=wall)


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
    out["teeth_ok"] = bool(out["captures_dropped"] == 0 and out["n_loss"] == 0 and out["n_dropout"] == 0 and
                           out["n_noise"] == 0 and out["sync"] == 2 and out["n_sync"] == host["status"]["n_sync"])
    out["rt_ok"] = bool(out["teeth_ok"] and out["slow_overruns"] == 0 and out["captures_unserviced"] <= 1)
    out["meets_spec"] = bool(out["angle_within_spec"] and out["rt_ok"])
    return out


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

    def min_clock(self, window: str, mode: str = "irq") -> dict:
        """Lowest core clock that meets the spec on `window`, by bisection (hil/cosim.yaml)."""
        mc = rv64.cosim_spec()["realtime"]["min_clock"]
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
    return (f"{m['window']:10s} {m['mode']:12s} {m['clk_hz'] / 1e6:5.1f} MHz {m['rpm']:5.0f} rpm | "
            f"tooth->event {m['tooth_to_event_max_us']:6.1f} us, ->serviced {m['tooth_to_serviced_max_us']:6.1f} us "
            f"({m['tooth_period_us']:.0f} us/tooth) | angle err: spark {m['spark_angle_error_max_deg']:.3f}, any event "
            f"{m['event_angle_error_max_deg']:.3f}, worst case {m['worst_case_angle_error_deg']:.3f} deg | load {100 * m['cpu_load']:5.1f} % | slow {m['slow_wall_max_us']:5.0f} us "
            f"max, {m['slow_steps']}/{m['slow_steps_expected']} steps, {m['slow_overruns']} overruns | dropped "
            f"{m['captures_dropped']}, loss {m['n_loss']}, sync {m['sync']} | real time {'HELD' if m['rt_ok'] else 'LOST'}, "
            f"angle spec {'MET' if m['meets_spec'] else 'NOT MET'}")


def main(argv=None) -> int:
    rtc = rv64.cosim_spec()["realtime"]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--windows", nargs="+", default=list(rtc["windows"]))
    ap.add_argument("--modes", nargs="+", default=["m9", "irq"], choices=list(MODES))
    ap.add_argument("--clocks", nargs="+", type=float, default=None, help="core clocks [Hz] (default: the board clock)")
    ap.add_argument("--sweep", action="store_true", help="clock sweep and minimum-clock bisection (hil/cosim.yaml)")
    a = ap.parse_args(argv)
    missing = rv64.toolchain_missing()
    if missing:
        print(f"missing tools: {missing}", file=sys.stderr)
        return 2
    bench = Bench()
    base = float(rv64.cosim_spec()["core_clock_hz"])
    clocks = a.clocks or ([float(c) for c in rtc["clocks_hz"]] if a.sweep else [base])
    print("Redline M9.1 real-time co-simulation: plant time = core cycles. All values MODELLED.")
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
    mine = [m for m in rows if m["mode"] == "irq" and m["clk_hz"] == base]
    ok = all(m["rt_ok"] for m in mine)
    spec = all(m["meets_spec"] for m in mine)
    print(f"RESULT at {base / 1e6:.0f} MHz (M9.1 firmware): real time "
          + ("HELD (no missed tooth, no loss of sync, no slow-loop overrun)" if ok else "LOST")
          + f"; {rtc['angle_error_max_deg']} deg angle spec " + ("MET" if spec else "NOT MET on every window (see the table)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
