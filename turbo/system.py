"""Redline I4-2.0T: the M4 1D engine with compressor, intercooler, turbine, wastegate and boost control.

Coupling (operator split, every `steady.substep_deg` crank degrees, gasdyn.engine `substep_fn`):

1. compressor: map flow at the current shaft speed against the charge-air plenum pressure (plus the
   intercooler and air-filter losses) is injected into the charge-air plenum at the intercooler outlet
   temperature (gasdyn.elements.MapCompressor);
2. turbine: the kernel's orifice boundary of nozzle-equivalent area carries the turbine flow; its
   isentropic work times the blade-speed-ratio efficiency is taken out of the post-turbine volume
   (gasdyn.elements.MapTurbine);
3. shaft: J w dw/dt = eta_m P_T - P_C (turbo.control.TurboShaft);
4. wastegate: PI on the filtered intake-plenum pressure sets the bypass orifice (turbo.control.BoostPI).

Steady points iterate whole engine cycles until the IMEP (gasdyn.engine), the cycle-mean turbo speed and
the boost error have settled. To get there in fewer cycles the rotor inertia is scaled by
`steady.inertia_factor` (a steady state does not depend on J); transients use the real inertia.
Spark is knock-limited (thermo.spark / thermo.knock) on the converged trapped state, as in sim.maps.
All results are MODELLED: no physical engine, synthetic turbo maps (CLAUDE.md rule 5).
"""
from __future__ import annotations

import math

import numpy as np

from gasdyn.elements import Intercooler, MapCompressor, MapTurbine, Wastegate, pipe_end_state
from gasdyn.engine import EngineModel
from gasdyn.friction import fmep_chen_flynn
from redline.spec import geometry
from sim.maps import M4Plant
from thermo.m7spec import combustion_spec
from thermo.props import build_tables, thermo_spec
from thermo.spark import SparkModel
from turbo.control import BoostPI, TurboShaft, rad_to_rpm, rpm_to_rad
from turbo.maps import CompressorMap, TurbineMap
from turbo.spec import engine_spec, gasdyn_variant_spec, turbo_spec

SEC_PER_MIN = 60.0
REV_PER_CYCLE = 2.0
CLOSED_EPS = 1.0e-3        # wastegate opening below which it counts as closed (saturated)


class CycleStats:
    """Time-weighted sums of the element outputs over one engine cycle."""

    KEYS = ("mdot_c", "pr_c", "eta_c", "P_C", "T2", "T_ic", "p_charge", "mdot_t", "mdot_wg", "pr_t", "eta_t", "bsr",
            "P_T", "T03", "p03", "p4", "map", "opening", "omega", "m_surge", "m_choke")

    def __init__(self):
        self.reset()

    def reset(self):
        self.t = 0.0
        self.s = dict.fromkeys(self.KEYS, 0.0)
        self.surge_time = 0.0
        self.sm_min = float("inf")
        self.e_c = 0.0          # compressor work over the cycle [J]
        self.e_t = 0.0

    def add(self, dt, vals: dict):
        self.t += dt
        for k in self.KEYS:
            self.s[k] += vals[k] * dt

    def mean(self) -> dict:
        return {k: v / self.t for k, v in self.s.items()} if self.t > 0.0 else dict(self.s)


class TurboEngine:
    def __init__(self, tspec: dict | None = None, tables=None, thermo: dict | None = None,
                 inertia_factor: float | None = None, spark_model: str | None = None):
        self.ts = tspec or turbo_spec()
        self.spec = engine_spec(self.ts)
        self.gspec = gasdyn_variant_spec(self.ts)
        self.thermo = thermo or thermo_spec()
        self.tables = tables or build_tables(self.thermo)
        self.geom = geometry(self.spec)
        self.st = self.ts["steady"]
        layout = dict(self.ts["layout"])
        layout["turbine_area"] = float(self.ts["turbine"]["equivalent_area"])
        layout["wastegate_area"] = 0.0
        self.eng = EngineModel(spec=self.spec, tspec=self.thermo, gspec=self.gspec, tables=self.tables, boost=layout)
        e = self.eng
        # The charge-air path carries the fresh mixture of the M4 intake pipes (port injection downstream;
        # the fuel share is lumped as in M4) so that the plenum energy bookkeeping stays consistent.
        self.cmap = CompressorMap(self.ts["compressor"], e.g_in, e.R_in)
        self.tmap = TurbineMap(self.ts["turbine"], e.g_ex, e.R_ex)
        af = self.ts["air_filter"]
        filt = lambda m, rho: float(af["dp_ref"]) * (m / float(af["m_ref"])) ** 2 * (float(af["rho_ref"]) / rho)
        ic = self.ts["intercooler"]
        self.ic = Intercooler(float(ic["effectiveness"]), float(ic["T_coolant"]), float(ic["dp_ref"]),
                              float(ic["m_ref"]), float(ic["rho_ref"]))
        self.comp = MapCompressor(self.cmap, e.v_charge, e.p_amb, e.T_amb, filt, self.ic)
        self.turb = MapTurbine(self.tmap, e.turbine_pipe, e.v_post_turbine)
        self.wg = Wastegate(e.wastegate_slot, float(self.ts["wastegate"]["max_area"]))
        self.inertia_factor = float(self.st["inertia_factor"]) if inertia_factor is None else inertia_factor
        self.shaft = TurboShaft(self.ts["shaft"], self.inertia_factor)
        self.ctrl = BoostPI(self.ts["boost_control"], self.ts["wastegate"])
        self.cspec = combustion_spec()
        self.sm = SparkModel(self.spec, self.thermo, self.cspec, self.gspec, self.tables, self.geom, spark_model)
        e.substep_fn = self._substep
        e.substep_deg = float(self.st["substep_deg"])
        e.cycle_fn = self._cycle
        self.stats = CycleStats()
        self.cycles: list[dict] = []
        self.log = None
        self.time = 0.0
        self.mode = "steady"
        self._ok = 0
        self._w_prev = None
        self._stop = None
        self.rpm = None

    # ---- state ------------------------------------------------------------------------------------
    def reset_turbo(self, rpm: float):
        self.shaft.omega = rpm_to_rad(float(self.st["init_speed_rpm"]))
        self.ctrl = BoostPI(self.ts["boost_control"], self.ts["wastegate"])
        self.ctrl.target = self.ctrl.target_for(rpm)

    def set_inertia_factor(self, f: float):
        self.shaft.J = float(self.ts["shaft"]["inertia"]) * f
        self.inertia_factor = f

    # ---- hooks ------------------------------------------------------------------------------------
    def _substep(self, net, theta_deg, dt):
        w = self.shaft.omega
        c = self.comp.step(net, w, dt)
        t = self.turb.step(net, w, dt)
        self.shaft.step(t["power"], c["power"], dt)
        p_map, _ = net.plenum_state(self.eng.v_plenum)
        opening = self.ctrl.step(p_map, dt)
        self.wg.set(net, opening)
        m_wg = pipe_end_state(net, self.eng.wastegate_pipe)[0]
        self.time += dt
        vals = {"mdot_c": c["mdot"], "pr_c": c["pr"], "eta_c": c["eta"], "P_C": c["power"], "T2": c["T2"],
                "T_ic": c["T_ic_out"], "p_charge": c["p_charge"], "mdot_t": t["mdot"], "mdot_wg": m_wg,
                "pr_t": t["pr"], "eta_t": t["eta"], "bsr": t["bsr"] if np.isfinite(t["bsr"]) else 0.0,
                "P_T": t["power"], "T03": t["T03"], "p03": t["p03"], "p4": t["p4"], "map": p_map,
                "opening": opening, "omega": self.shaft.omega, "m_surge": c["m_surge"], "m_choke": c["m_choke"]}
        self.stats.add(dt, vals)
        if c["surge"]:
            self.stats.surge_time += dt
        if c["m_surge"] > 0.0:
            self.stats.sm_min = min(self.stats.sm_min, c["surge_margin"])
        if self.log is not None:
            self.log.append((self.time, p_map, self.ctrl.map_f, self.shaft.omega, opening, c["mdot"], c["pr"],
                             c["p_charge"]))

    def _brake(self, rec, rpm):
        imep = float(np.mean([r["imep"] for r in rec]))
        p_max = float(np.mean([r["p_max"] for r in rec]))
        sp = 2.0 * self.geom.stroke * rpm / SEC_PER_MIN
        bmep = imep - fmep_chen_flynn(p_max, sp, self.gspec["friction"])
        return bmep * self.geom.displacement / (2.0 * math.pi * REV_PER_CYCLE)

    def _cycle(self, icyc, rec):
        m = self.stats.mean()
        m.update(cycle=icyc, time_end=self.time, surge_time=self.stats.surge_time, sm_min_inst=self.stats.sm_min,
                 torque_nm=self._brake(rec, self.rpm), speed_rpm=rad_to_rpm(m["omega"]),
                 target=self.ctrl.target if self.ctrl.target is not None else float("nan"))
        self.cycles.append(m)
        self.stats.reset()
        if self.mode != "steady":
            if self._stop is None:
                return False
            level, hold, n = self._stop
            n = n + 1 if m["map"] >= level else 0
            self._stop = (level, hold, n)
            return n >= hold
        w = m["omega"]
        dw = abs(w / self._w_prev - 1.0) if self._w_prev else float("inf")
        self._w_prev = w
        tgt = self.ctrl.target
        if tgt is None:
            ctrl_ok = True
        else:
            err = m["map"] / tgt - 1.0
            ctrl_ok = (abs(err) < float(self.st["map_tol"]) or (m["opening"] < CLOSED_EPS and err < 0.0)
                       or (m["opening"] > 1.0 - CLOSED_EPS and err > 0.0))
        self._ok = self._ok + 1 if (dw < float(self.st["speed_tol"]) and ctrl_ok) else 0
        return self._ok >= int(self.st["consecutive"]) and icyc + 1 >= self._min_cycles

    # ---- runs -------------------------------------------------------------------------------------
    def run(self, rpm, state=None, throttle_deg=None, closed_fn=None, max_cycles=None, mode="steady", trace=False):
        self.rpm = float(rpm)
        self.mode = mode
        self._min_cycles = int(self.st["min_cycles"] if state is None else self.st["min_cycles_warm"])
        self._ok = 0
        self._w_prev = None
        self.stats.reset()
        self.eng.closed_fn = closed_fn
        n = int(max_cycles or self.st["max_cycles"])
        return self.eng.run(float(rpm), state=state, throttle_deg=throttle_deg, max_cycles=n, trace=trace)

    def summary(self, op) -> dict:
        """Operating point metrics: M4 engine metrics + turbo / boost metrics of the last cycle."""
        c = self.cycles[-1]
        out = dict(op.metrics)
        # corrected flow referred to the spec ambient state (compressor inlet T = ambient in this model)
        p1, T1 = self.comp.inlet(c["mdot_c"])
        m_corr = c["mdot_c"] * math.sqrt(T1 / self.eng.T_amb) * (self.eng.p_amb / p1)
        rho_map = c["map"] / (self.eng.R_in * c["T_ic"])
        out.update({
            "map_bar": c["map"] / 1e5, "target_bar": c["target"] / 1e5,
            "boost_error": c["map"] / c["target"] - 1.0 if np.isfinite(c["target"]) else float("nan"),
            "p_charge_bar": c["p_charge"] / 1e5, "turbo_krpm": c["speed_rpm"] / 1e3, "wg_opening": c["opening"],
            "pr_c": c["pr_c"], "mdot_c": c["mdot_c"], "mdot_c_corr": m_corr, "eta_c": c["eta_c"],
            "m_surge": c["m_surge"], "surge_margin": c["mdot_c"] / c["m_surge"] - 1.0 if c["m_surge"] > 0 else float("inf"),
            "surge_margin_inst_min": c["sm_min_inst"], "surge_time": c["surge_time"],
            "P_C_kw": c["P_C"] / 1e3, "P_T_kw": c["P_T"] / 1e3, "T2_K": c["T2"], "T_ic_K": c["T_ic"],
            "T03_K": c["T03"], "p03_bar": c["p03"] / 1e5, "p4_bar": c["p4"] / 1e5, "pr_t": c["pr_t"],
            "eta_t": c["eta_t"], "bsr": c["bsr"], "mdot_t": c["mdot_t"], "mdot_wg": c["mdot_wg"],
            "ve_plenum": out["ve_trapped"] * self.eng.rho_air_amb / rho_map,
            "shaft_balance": (self.shaft.eta_m * c["P_T"] - c["P_C"]) / max(c["P_C"], 1.0),
            "overspeed": self.shaft.overspeed})
        return out

    def steady_point(self, rpm: float, knock_limited: bool = True) -> dict:
        """Converged WOT point with boost control and knock-limited spark (cold start)."""
        self.set_inertia_factor(float(self.st["inertia_factor"]))
        self.reset_turbo(rpm)
        self.cycles = []
        wot = float(self.gspec["intake"]["throttle"]["wot_angle_deg"])
        op = self.run(rpm, None, wot)
        opt = None
        if knock_limited:
            for _ in range(int(self.st["spark_passes"]) + 1):
                stp = M4Plant.trapped(op, self.eng)
                opt = self.sm.optimise(stp, rpm)
                op = self.run(rpm, op.state, wot, self.sm.closed_fn_engine(opt["spark_deg"]))
            stp = M4Plant.trapped(op, self.eng)
            ev = self.sm.evaluate(stp, rpm, opt["spark_deg"])
        out = self.summary(op)
        out["converged"] = bool(op.converged)
        out["cycles_run"] = len(self.cycles)
        if opt is not None:
            out.update({"spark_deg": opt["spark_deg"], "mbt_deg": opt["mbt_deg"], "klsa_deg": opt["klsa_deg"],
                        "klsa_status": opt["klsa_status"], "knock_limited": opt["knock_limited"],
                        "knock_margin": ev["knock"]["margin"], "ca50_deg": ev["ca50"]})
        else:
            out.update({"spark_deg": float(self.spec["combustion"]["spark_advance_deg"]), "mbt_deg": float("nan"),
                        "klsa_deg": float("nan"), "klsa_status": "not evaluated", "knock_limited": False,
                        "knock_margin": float("nan"), "ca50_deg": float("nan")})
        self.eng.closed_fn = None
        return out

    def wot_spark_estimate(self, rpm: float) -> float:
        """Knock-limited spark from one optimisation on the Wiebe-converged steady WOT state (the first pass of
        steady_point, without the flame-model re-runs): a cheaper spark for transients."""
        self.set_inertia_factor(float(self.st["inertia_factor"]))
        self.reset_turbo(rpm)
        self.cycles = []
        wot = float(self.gspec["intake"]["throttle"]["wot_angle_deg"])
        op = self.run(rpm, None, wot)
        return float(self.sm.optimise(M4Plant.trapped(op, self.eng), rpm)["spark_deg"])

    def tipin(self, spark_deg: float | None = None, cfg: dict | None = None) -> dict:
        """Dyno-held tip-in: part throttle (no boost demand) -> WOT with the boost target, real rotor inertia.

        Spark is held at `spark_deg` (the knock-limited spark of the steady WOT point at this speed: the
        most knock-prone state, so conservative while the boost builds); None = spec Wiebe timing. The
        recorded window ends `duration_s` after the step, or earlier once the cycle-mean MAP has stayed
        above both 90 % levels for `stop_hold_cycles` cycles."""
        cfg = cfg or self.ts["tipin"]
        rpm = float(cfg["rpm"])
        fn = self.sm.closed_fn_engine(spark_deg) if spark_deg is not None else None
        # settle at part load (scaled inertia only to reach the part-load steady state)
        self.set_inertia_factor(float(self.st["inertia_factor"]))
        self.reset_turbo(rpm)
        self.ctrl.target = None
        self.cycles = []
        op = self.run(rpm, None, float(cfg["throttle_deg"]), fn, max_cycles=int(cfg["settle_cycles"]))
        # recorded window with the real inertia
        self.set_inertia_factor(1.0)
        self.time = 0.0
        self.log = []
        i0 = len(self.cycles)
        cyc_s = SEC_PER_MIN * REV_PER_CYCLE / rpm
        n_pre = max(1, int(round(float(cfg["pre_s"]) / cyc_s)))
        op = self.run(rpm, op.state, float(cfg["throttle_deg"]), fn, max_cycles=n_pre, mode="transient")
        t_step = self.time
        p0 = float(np.mean([c["map"] for c in self.cycles[i0:]]))
        self.ctrl.target = self.ctrl.target_for(rpm)
        target = self.ctrl.target
        frac = float(cfg["boost_fraction"])
        lvl_rise = p0 + frac * (target - p0)
        lvl_abs = frac * target
        self._stop = (max(lvl_rise, lvl_abs), int(cfg["stop_hold_cycles"]), 0)
        wot = float(self.gspec["intake"]["throttle"]["wot_angle_deg"])
        i1 = len(self.cycles)
        n_post = int(math.ceil(float(cfg["duration_s"]) / cyc_s))
        op = self.run(rpm, op.state, wot, fn, max_cycles=n_post, mode="transient")
        self._stop = None
        log = np.array(self.log)
        self.log = None
        self.eng.closed_fn = None
        cyc = self.cycles[i0:]
        t_mid = np.array([c["time_end"] - cyc_s / 2.0 for c in cyc]) - t_step
        map_c = np.array([c["map"] for c in cyc])
        tq = np.array([c["torque_nm"] for c in cyc])
        post = self.cycles[i1:]
        return {"rpm": rpm, "t": t_mid, "map": map_c, "torque": tq, "map_pre": p0, "target": target,
                "t90_rise": _first_crossing(t_mid, map_c, lvl_rise), "t90_abs": _first_crossing(t_mid, map_c, lvl_abs),
                "map_final": float(map_c[-1]), "torque_pre": float(np.mean(tq[t_mid < 0.0])),
                "torque_final": float(tq[-1]), "log": log, "t_step": t_step, "boost_fraction": frac,
                "spark_deg": spark_deg, "surge_time": float(sum(c["surge_time"] for c in post)),
                "window_s": float(self.time - t_step),
                "surge_margin_cycle_min": float(min(c["mdot_c"] / c["m_surge"] - 1.0 for c in post)),
                "speed_final_krpm": post[-1]["speed_rpm"] / 1e3}


def _first_crossing(t, y, level):
    """First time at which the piecewise-linear y(t) reaches level, clipped at t >= 0; NaN if never."""
    for i in range(1, len(t)):
        if y[i] >= level and y[i - 1] < level:
            tc = float(t[i - 1] + (level - y[i - 1]) / (y[i] - y[i - 1]) * (t[i] - t[i - 1]))
            if tc >= 0.0 or t[i] >= 0.0:
                return max(tc, 0.0)
    return float("nan")
