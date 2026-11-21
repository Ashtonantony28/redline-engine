"""Hendricks-style mean-value engine model (MVEM) with event-resolved fuel and spark. MODELLED.

State and sub-models (parameters: spec/engine.yaml + spec/mvem.yaml, ADR 0001 conventions):
  * throttle: compressible orifice, butterfly area A(a) = A_bore (1 - cos(psi0 + a) / cos psi0) + A_leak
  * intake manifold: isothermal filling, dp/dt = R T / V (mdot_th - mdot_port), implicit step
  * port air flow: speed density, mdot_port = VE(N, p) p / (R T) Vd N / 120 (VE from thermo/ maps)
  * fuel: per-port Aquino x-tau film; injected mass = q (pw - dead time(Vbatt))
  * cylinder events: charge captured at IVC (air inducted since the previous IVC + port fuel),
    work booked after the spark from gross indicated efficiency eta(N, p, SOC) (thermo/ runs) minus
    pumping work PMEP(N, p) Vd; spark sets SOC = ignition delay - advance
  * crank: J domega/dt = T_ind - T_fmep - T_load, J from dyn.flywheel.rotating_inertia
  * electronic throttle: DC motor (L neglected) + gear + preloaded return/limp springs + friction
  * wideband lambda: transport delay (fixed + engine cycles) then first-order lag
  * trigger: 60-2 crank wheel and one cam edge per cycle (crossing angles); knock sensor window
    integral per combustion event with a simulated knock limit and a test hook
"""
from __future__ import annotations

import bisect
import math
import random
from collections import deque
from dataclasses import dataclass

from mvem import maps as mmaps
from mvem import params

TWO_PI = 2.0 * math.pi
DEG_PER_RAD = 180.0 / math.pi


def _seg(xs, x):
    """Index i and fraction f of x in the sorted axis xs (clamped)."""
    n = len(xs)
    if x <= xs[0]:
        return 0, 0.0
    if x >= xs[-1]:
        return n - 2, 1.0
    i = bisect.bisect_right(xs, x) - 1
    return i, (x - xs[i]) / (xs[i + 1] - xs[i])


def bilinear(xs, ys, z, x, y):
    i, fx = _seg(xs, x)
    j, fy = _seg(ys, y)
    a = z[i][j] + fy * (z[i][j + 1] - z[i][j])
    b = z[i + 1][j] + fy * (z[i + 1][j + 1] - z[i + 1][j])
    return a + fx * (b - a)


def interp1(xs, ys, x):
    i, f = _seg(xs, x)
    return ys[i] + f * (ys[i + 1] - ys[i])


def orifice_flow(cd_area: float, p_up: float, T_up: float, p_down: float, R: float, gamma: float) -> float:
    """Isentropic compressible orifice mass flow [kg/s] (Heywood App. C), choked below p*/p0."""
    if p_down > p_up:
        return -orifice_flow(cd_area, p_down, T_up, p_up, R, gamma)
    pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    pr = max(p_down / p_up, pr_crit)
    psi = math.sqrt(2.0 * gamma / (gamma - 1.0) * (pr ** (2.0 / gamma) - pr ** ((gamma + 1.0) / gamma)))
    return cd_area * p_up / math.sqrt(R * T_up) * psi


@dataclass
class Charge:
    air: float = 0.0
    fuel: float = 0.0
    p_map: float = 0.0
    rpm: float = 0.0


class Engine:
    """MVEM state + sub-models. Time advances through `advance()` (continuous states) while the
    caller (hil.runner) delivers crank-angle events in time order: injections, sparks, IVC and
    combustion crossings, and reads trigger crossings from `crossings()`."""

    def __init__(self, rpm0: float, p_map0: float, throttle0_deg: float | None = None,
                 ve_override: float | None = None):
        es = params.engine_spec()
        ms = params.mvem_spec()
        g = params.geometry()
        self.ms = ms
        self.ncyl = g.cylinders
        self.vd_cyl = g.swept_volume
        self.vd = g.displacement
        self.stroke = g.stroke
        amb, thr, man = ms["ambient"], ms["throttle"], ms["manifold"]
        self.p_amb, self.T_amb = float(amb["p"]), float(amb["T"])
        self.R, self.gamma = float(amb["R_air"]), float(amb["gamma_air"])
        self.A_bore = math.pi * float(thr["bore_diameter"]) ** 2 / 4.0
        self.psi0 = math.radians(float(thr["closed_plate_angle_deg"]))
        self.thr_max = float(thr["max_opening_deg"])
        self.cd = float(thr["discharge_coefficient"])
        self.A_leak = float(thr["leak_area"])
        self.V_man = float(man["volume"])
        self.T_man = float(man["temperature"])
        fr = ms["friction"]
        self.fmep_c = (float(fr["fmep_c0"]), float(fr["fmep_c1"]), float(fr["fmep_c2"]))
        cb = ms["combustion"]
        self.lam_lean, self.lam_rich = float(cb["lean_misfire_lambda"]), float(cb["rich_misfire_lambda"])
        self.eval_deg = float(cb["evaluation_deg_atdc"])
        self.lhv = float(es["combustion"]["lhv"])
        self.afr = float(es["combustion"]["afr_stoich"])
        ff = ms["fuel_film"]
        self.film_x = float(ff["x"])
        self.tau_ref, self.tau_rpm_ref, self.tau_exp = (float(ff["tau_ref"]), float(ff["rpm_ref"]),
                                                       float(ff["tau_rpm_exponent"]))
        inj = ms["injector"]
        self.q_inj = float(inj["static_flow"]) * (1.0 + float(inj["flow_error_fraction"]))
        self.dead_v = [float(x) for x in inj["deadtime_vbatt"]]
        self.dead_s = [float(x) for x in inj["deadtime_s"]]
        self.vbatt = float(inj["battery_voltage"])
        et = ms["etc"]
        self.etc = {k: float(v) for k, v in et.items()}
        self.etc_sub = int(et["substeps"])
        ls = ms["lambda_sensor"]
        self.lam_delay_s, self.lam_delay_cyc = float(ls["transport_delay_s"]), float(ls["transport_delay_cycles"])
        self.lam_tau, self.lam_min, self.lam_max = float(ls["time_constant"]), float(ls["min"]), float(ls["max"])
        kn = ms["knock"]
        self.kn = {k: float(v) for k, v in kn.items()}
        self.rng = random.Random(int(kn["seed"]))
        trg = ms["trigger"]
        self.teeth, self.missing = int(trg["teeth_total"]), int(trg["teeth_missing"])
        self.tooth0, self.cam_edge = float(trg["tooth0_angle_deg"]), float(trg["cam_edge_deg"])
        self.tdc = params.firing_tdc_deg(es)
        self.ivc_rel = float(es["valves"]["ivc_deg"]) - params.CYCLE_DEG   # IVC relative to firing TDC
        self.ign_delay = params.ignition_delay_deg(es)
        from dyn.flywheel import rotating_inertia
        self.J = rotating_inertia()

        m = mmaps.load_maps()
        self.ax_rpm = [float(x) for x in m["rpm"]]
        self.ax_map = [float(x) for x in m["map"]]
        self.ax_soc = [float(x) for x in m["soc_deg"]]
        self.ve_tab = m["ve"].tolist()
        self.pmep_tab = m["pmep"].tolist()
        self.eta_tab = m["eta_gross"].tolist()
        self.ve_override = ve_override
        self._build_cycle_table()

        # state
        self.t = 0.0
        self.theta = 0.0                       # unwrapped cycle angle [deg]
        self.omega = rpm0 * params.RPM_TO_RAD_S
        self.p_man = p_map0
        self.thr = self.etc["limp_opening_deg"] if throttle0_deg is None else throttle0_deg
        self.thr_rate = 0.0
        self.film = [0.0] * self.ncyl
        self.port = [0.0] * self.ncyl
        self.air_acc = 0.0
        self.charge = [Charge() for _ in range(self.ncyl)]
        self.spark_adv: list[float | None] = [None] * self.ncyl
        self.T_ind = 0.0
        self.lam_queue: deque = deque()
        self.lam_sensed_target = 1.0
        self.lam_sensor = 1.0
        self.last_lambda = 1.0
        self.knock_hook = None                 # callable(t, cyl) -> bool: force a knocking event
        self.knock_offset_deg = 0.0            # lowers the simulated knock limit by this much (test hook)
        self.t_flow = 0.0
        self.mdot_port = self.port_flow(rpm0, p_map0)
        self.T_brake = 0.0
        self.events: list[dict] = []           # combustion log

    # ------------------------------------------------------------------ sub-models
    def rpm(self) -> float:
        return self.omega / params.RPM_TO_RAD_S

    def throttle_area(self, alpha_deg: float) -> float:
        a = min(max(alpha_deg, 0.0), self.thr_max)
        return self.A_bore * (1.0 - math.cos(self.psi0 + math.radians(a)) / math.cos(self.psi0)) + self.A_leak

    def throttle_flow(self, alpha_deg: float, p_man: float) -> float:
        return orifice_flow(self.cd * self.throttle_area(alpha_deg), self.p_amb, self.T_amb, p_man, self.R, self.gamma)

    def ve(self, rpm: float, p_man: float) -> float:
        if self.ve_override is not None:
            return self.ve_override
        return bilinear(self.ax_rpm, self.ax_map, self.ve_tab, rpm, p_man)

    def port_flow(self, rpm: float, p_man: float) -> float:
        return (self.ve(rpm, p_man) * p_man / (self.R * self.T_man) * self.vd * rpm
                / (params.SEC_PER_MIN * params.REV_PER_CYCLE))

    def fmep(self, rpm: float) -> float:
        sp = 2.0 * self.stroke * rpm / params.SEC_PER_MIN
        c0, c1, c2 = self.fmep_c
        return c0 + c1 * sp + c2 * sp * sp

    def friction_torque(self, rpm: float) -> float:
        return self.fmep(rpm) * self.vd / (2.0 * TWO_PI)

    def film_tau(self, rpm: float) -> float:
        return self.tau_ref * (max(rpm, 1.0) / self.tau_rpm_ref) ** self.tau_exp

    def deadtime(self) -> float:
        return interp1(self.dead_v, self.dead_s, self.vbatt)

    def klsa(self, rpm: float, p_man: float) -> float:
        k = self.kn
        return (k["klsa_wot_deg"] + k["klsa_per_unit_vacuum_deg"] * (1.0 - p_man / self.p_amb)
                + k["klsa_per_krpm_deg"] * (rpm - k["klsa_rpm_ref"]) / 1000.0)

    def eta_gross(self, rpm: float, p_man: float, soc: float) -> float:
        i, fr = _seg(self.ax_rpm, rpm)
        j, fm = _seg(self.ax_map, p_man)
        k, fs = _seg(self.ax_soc, soc)
        e = self.eta_tab

        def lin(a, b):
            return a + fs * (b - a)
        c00 = lin(e[i][j][k], e[i][j][k + 1])
        c01 = lin(e[i][j + 1][k], e[i][j + 1][k + 1])
        c10 = lin(e[i + 1][j][k], e[i + 1][j][k + 1])
        c11 = lin(e[i + 1][j + 1][k], e[i + 1][j + 1][k + 1])
        a = c00 + fm * (c01 - c00)
        b = c10 + fm * (c11 - c10)
        return max(a + fr * (b - a), 0.0)       # residual-dominated corner of the map: no negative gross work

    def pmep(self, rpm: float, p_man: float) -> float:
        return bilinear(self.ax_rpm, self.ax_map, self.pmep_tab, rpm, p_man)

    # ------------------------------------------------------------------ trigger / event angles
    def _build_cycle_table(self):
        pitch = params.REV_DEG / self.teeth
        ev = []
        for rev in range(int(params.REV_PER_CYCLE)):
            for k in range(self.teeth - self.missing):
                ev.append(((self.tooth0 + k * pitch) % params.REV_DEG + rev * params.REV_DEG, "tooth", k))
        ev.append((self.cam_edge % params.CYCLE_DEG, "cam", 0))
        for c in range(self.ncyl):
            ev.append(((self.tdc[c] + self.ivc_rel) % params.CYCLE_DEG, "ivc", c))
            ev.append(((self.tdc[c] + self.eval_deg) % params.CYCLE_DEG, "burn", c))
        ev.sort(key=lambda e: e[0])
        self.cyc_angles = [e[0] for e in ev]
        self.cyc_events = ev

    def crossings(self, th0: float, th1: float):
        """Event angles in (th0, th1] (unwrapped cycle angle): list of (theta, kind, index)."""
        out = []
        base = math.floor(th0 / params.CYCLE_DEG) * params.CYCLE_DEG
        a0 = th0 - base
        i = bisect.bisect_right(self.cyc_angles, a0)
        n = len(self.cyc_angles)
        while True:
            if i >= n:
                i = 0
                base += params.CYCLE_DEG
            th = base + self.cyc_angles[i]
            if th > th1:
                break
            e = self.cyc_events[i]
            out.append((th, e[1], e[2]))
            i += 1
        return out

    # ------------------------------------------------------------------ events
    def advance_flows(self, t: float):
        """Integrate inducted air and film evaporation from t_flow to t (mdot_port held)."""
        dt = t - self.t_flow
        if dt <= 0.0:
            return
        self.air_acc += self.mdot_port * dt
        decay = math.exp(-dt / self.film_tau(self.rpm()))
        for c in range(self.ncyl):
            rel = self.film[c] * (1.0 - decay)
            self.film[c] -= rel
            self.port[c] += rel
        self.t_flow = t

    def inject(self, cyl: int, pw_s: float):
        m = self.q_inj * max(pw_s - self.deadtime(), 0.0)
        self.film[cyl] += self.film_x * m
        self.port[cyl] += (1.0 - self.film_x) * m

    def spark(self, cyl: int, theta_spark: float):
        adv = (self.tdc[cyl] - theta_spark) % params.CYCLE_DEG
        if adv > params.REV_DEG:
            adv -= params.CYCLE_DEG
        self.spark_adv[cyl] = adv

    def ivc(self, cyl: int, t: float):
        self.advance_flows(t)
        self.charge[cyl] = Charge(self.air_acc, self.port[cyl], self.p_man, self.rpm())
        self.air_acc = 0.0
        self.port[cyl] = 0.0

    def burn(self, cyl: int, t: float):
        """Book the work of cylinder `cyl`; returns (lambda, knock window counts)."""
        ch = self.charge[cyl]
        adv = self.spark_adv[cyl]
        self.spark_adv[cyl] = None
        lam = ch.air / (ch.fuel * self.afr) if ch.fuel > 0.0 else math.inf
        fired = adv is not None and self.lam_rich <= lam <= self.lam_lean and ch.air > 0.0
        W = -self.pmep(ch.rpm, ch.p_map) * self.vd_cyl
        counts = self.kn["noise_mean_counts"] + self.kn["noise_counts"] * self.rng.gauss(0.0, 1.0)
        if fired:
            burned = min(ch.fuel, ch.air / self.afr)
            soc = self.ign_delay - adv
            W += self.eta_gross(ch.rpm, ch.p_map, soc) * burned * self.lhv
            excess = adv - (self.klsa(ch.rpm, ch.p_map) - self.knock_offset_deg)
            if excess > 0.0:
                counts += self.kn["gain_counts_per_deg"] * excess
            if self.knock_hook is not None and self.knock_hook(t, cyl):
                counts += self.kn["injected_counts"]
        self.T_ind = W / (params.CYCLE_DEG / self.ncyl / DEG_PER_RAD)
        self.last_lambda = lam
        rpm = max(self.rpm(), 1.0)
        delay = self.lam_delay_s + self.lam_delay_cyc * params.SEC_PER_MIN * params.REV_PER_CYCLE / rpm
        self.lam_queue.append((t + delay, min(max(lam, self.lam_min), self.lam_max)))
        self.events.append({"t": t, "cyl": cyl, "lambda": lam, "adv": adv, "W": W, "fired": fired,
                            "air": ch.air, "fuel": ch.fuel})
        return lam, max(int(round(counts)), 0)

    # ------------------------------------------------------------------ continuous states
    def etc_step(self, volts: float, dt: float):
        e = self.etc
        h = dt / self.etc_sub
        J, R, kt, ke = e["inertia"], e["motor_resistance"], e["torque_constant"], e["backemf_constant"]
        b = e["viscous"] + kt * ke / R
        lim = e["limp_opening_deg"]
        for _ in range(self.etc_sub):
            off = math.radians(self.thr - lim)
            spring = -(e["spring_rate"] * off + e["spring_preload"] * ((off > 0) - (off < 0)))
            F = kt * volts / R + spring
            w = self.thr_rate
            if w == 0.0 and abs(F) <= e["static_friction"]:
                continue                       # stuck: breakaway torque not reached
            mom = J * w + h * F                # momentum after the step without Coulomb friction
            if abs(mom) <= h * e["coulomb_friction"]:
                w_new = 0.0                    # Coulomb friction stops the plate within the substep
            else:
                sgn = 1.0 if mom > 0.0 else -1.0
                w_new = (mom - h * e["coulomb_friction"] * sgn) / (J + h * b)
            self.thr += math.degrees(h * w_new)
            self.thr_rate = w_new
            if self.thr <= 0.0 or self.thr >= self.thr_max:
                self.thr = min(max(self.thr, 0.0), self.thr_max)
                self.thr_rate = 0.0

    def manifold_step(self, dt: float, rpm: float):
        """Implicit Euler on p_man (throttle flow is stiff near pressure ratio 1)."""
        k = self.R * self.T_man / self.V_man
        p0 = self.p_man
        area = self.cd * self.throttle_area(self.thr)

        def g(p):
            return p - p0 - dt * k * (orifice_flow(area, self.p_amb, self.T_amb, p, self.R, self.gamma)
                                      - self.port_flow(rpm, p))
        lo, hi = 1.0, self.p_amb
        glo, ghi = g(lo), g(hi)
        if ghi <= 0.0:
            self.p_man = hi
            return
        # Illinois regula falsi: g is monotone increasing in p
        side = 0
        p = p0
        for _ in range(60):
            p = (lo * ghi - hi * glo) / (ghi - glo)
            gp = g(p)
            if abs(gp) < 1.0e-3 or hi - lo < 1.0e-3:
                break
            if gp > 0.0:
                hi, ghi = p, gp
                if side == -1:
                    glo *= 0.5
                side = -1
            else:
                lo, glo = p, gp
                if side == 1:
                    ghi *= 0.5
                side = 1
        self.p_man = p

    def lambda_step(self, dt: float):
        while self.lam_queue and self.lam_queue[0][0] <= self.t:
            self.lam_sensed_target = self.lam_queue.popleft()[1]
        self.lam_sensor += (self.lam_sensed_target - self.lam_sensor) * (1.0 - math.exp(-dt / self.lam_tau))

    def advance(self, dt: float, etc_volts: float | None = None, load_torque: float = 0.0,
                dyno_rpm: float | None = None, throttle_deg: float | None = None):
        """Continuous states over [t, t + dt]. Event crossings inside the step must already have
        been delivered by the caller (they only use theta and the flows)."""
        rpm = self.rpm()
        self.advance_flows(self.t + dt)
        self.t += dt
        if throttle_deg is not None:
            self.thr = throttle_deg
        elif etc_volts is not None:
            self.etc_step(etc_volts, dt)
        self.manifold_step(dt, rpm)
        self.mdot_port = self.port_flow(rpm, self.p_man)
        self.lambda_step(dt)
        self.T_brake = self.T_ind - self.friction_torque(rpm) - load_torque
        self.theta += self.omega * DEG_PER_RAD * dt
        if dyno_rpm is not None:
            self.omega = dyno_rpm * params.RPM_TO_RAD_S
        else:
            self.omega = max(self.omega + dt * self.T_brake / self.J, 0.0)

    def tps(self):
        """Two complementary throttle position sensors [deg]."""
        return self.thr, self.thr_max - self.thr
