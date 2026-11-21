"""Four-cylinder engine on the 1D intake / exhaust network (milestone M4).

Intake: ambient -> inlet duct -> throttle (variable-area orifice) -> plenum (0D) -> 4 runners -> valves.
Exhaust: valves -> 4 primaries -> 4-2-1 (or 4-1) junctions as 0D volumes -> tailpipe -> silencer -> ambient.

Each cylinder (ADR 0001: firing order from spec, 180 deg apart, local angle 0 = its firing TDC) is

* a 0D volume inside the gas-dynamics kernel while a valve is off its seat (kin.cam seat-to-seat events,
  spec timing interpreted at the check lift), exchanging mass with its runner and primary through the
  valve boundary of gasdyn.fv with the real kin.cam lift and the spec discharge coefficients, and
* the thermo two-zone closed model (thermo.cycle.run_closed, Wiebe burn, equilibrium burned gas, Woschni)
  from seat IVC to seat EVO. The closed part does not interact with the pipes, so it is evaluated in one
  call at IVC and its EVO state is handed back to the kernel when the crank reaches EVO.

Cycles are repeated until the 4-cylinder mean net IMEP changes by less than `cycle.imep_tol` between
consecutive cycles. All results are MODELLED: there is no physical engine (CLAUDE.md rule 5).

Turbocharger hooks (M10 stretch, turbo/; all optional, the defaults build and run exactly the M4 engine):

* `boost` (constructor): layout of a turbocharged variant. The intake ambient becomes a 0D charge-air
  plenum (compressor outlet + intercooler + ducting to the throttle) and the exhaust collector feeds a
  turbine-housing pipe and a wastegate pipe, both ending in variable orifices into a 0D post-turbine
  volume that drains through the catalyst / silencer line. The compressor (mass / enthalpy source), the
  turbine (work sink) and the wastegate area are applied by gasdyn.elements between kernel calls.
* `substep_fn(net, theta_deg, dt)` with `substep_deg` (attributes): called after every kernel call; the
  marching is additionally split every `substep_deg` crank degrees.
* `cycle_fn(icyc, rec)` (attribute): called after every complete cycle; returning False vetoes the IMEP
  convergence test (used to iterate to a turbo / boost-control steady state or to run a fixed-length
  transient).
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field

import numpy as np

from gasdyn import fv
from gasdyn.friction import fmep_jets, fmep_total
from gasdyn.network import Network, numerics_vector
from gasdyn.spec import gasdyn_spec, valvetrain_spec, with_overrides
from kin.cam import design_cams
from kin.slider_crank import cylinder_volume, dV_dtheta
from redline.spec import geometry, load
from thermo import twozone as tz
from thermo.cycle import CYCLE_DEG, effective_valve_area, merge_zones, run_closed
from thermo.props import Tables, build_tables, frozen, thermo_spec, T_from_u_frozen
from thermo.woschni import heat_params

DEG = math.pi / 180.0
FIRING_INTERVAL_DEG = 180.0     # even-firing four-stroke inline-4 (ADR 0001)
REV_PER_CYCLE = 2.0             # four-stroke
SEC_PER_MIN = 60.0


def throttle_area(angle_deg: float, bore: float, closed_angle_deg: float) -> float:
    """Projected open area of a butterfly plate: (pi D^2 / 4) (1 - cos(a) / cos(a0)), the leading term
    of Heywood eq. 7.18 (throttle shaft neglected; its blockage is lumped into the throttle Cd)."""
    a = min(max(angle_deg, closed_angle_deg), 90.0) * DEG
    return 0.25 * math.pi * bore**2 * (1.0 - math.cos(a) / math.cos(closed_angle_deg * DEG))


@dataclass
class EngineState:
    """Everything needed to continue a run (used to warm-start the next speed of a sweep)."""
    U: np.ndarray
    S: np.ndarray
    Tg: np.ndarray
    act: np.ndarray
    pend: np.ndarray          # per cylinder: m, m x, U, T at EVO (from the last closed-part run)
    w_closed: np.ndarray      # per cylinder: work of the last closed part [J]
    p_max: np.ndarray
    q_closed: np.ndarray

    def copy(self):
        return copy.deepcopy(self)


@dataclass
class OperatingPoint:
    rpm: float
    throttle_deg: float
    converged: bool
    cycles: int
    metrics: dict
    per_cylinder: list
    history: list
    state: EngineState
    trace_deg: np.ndarray | None = None
    trace: dict = field(default_factory=dict)


class EngineModel:
    def __init__(self, spec: dict | None = None, tspec: dict | None = None, vspec: dict | None = None,
                 gspec: dict | None = None, tables: Tables | None = None, overrides: dict | None = None,
                 vvt_advance_deg: float = 0.0, boost: dict | None = None):
        self.spec = spec or load()
        self.tspec = tspec or thermo_spec()
        self.vspec = vspec or valvetrain_spec()
        self.gspec = with_overrides(gspec or gasdyn_spec(), overrides)
        self.tables = tables or build_tables(self.tspec)
        self.geom = geometry(self.spec)
        self.targs = self.tables.args()
        self.props = self.targs[:5]
        num_t = self.tspec["numerics"]
        self.tol, self.maxit = float(num_t["newton_tol"]), int(num_t["newton_max_iter"])
        self.vvt_advance_deg = vvt_advance_deg
        self.boost = copy.deepcopy(boost) if boost else None
        # turbo hooks (module docstring); None keeps the M4 behaviour bit for bit
        self.substep_fn = None
        self.substep_deg = None
        self.cycle_fn = None
        # M7 hook: replacement for the spec Wiebe closed part, callable(m, x, T_ivc, rpm, hp, theta_s, theta_e,
        # step_deg) -> thermo.cycle.ClosedResult (e.g. a spark-timed flame model). None = spec Wiebe (M4).
        self.closed_fn = None
        self._valve_tables()
        self._gas()
        self._build()

    # ------------------------------------------------------------------------------------------
    def _valve_tables(self):
        """Cylinder volume and effective valve areas on a uniform local crank-angle grid (one cycle)."""
        val, gex = self.spec["valves"], self.tspec["gas_exchange"]
        step = float(self.gspec["numerics"]["table_step_deg"])
        n = int(round(CYCLE_DEG / step))
        th = np.arange(n) * (CYCLE_DEG / n) * DEG
        self.cams = design_cams(self.spec, self.vspec)
        adv = self.vvt_advance_deg * DEG
        nv = val["per_cylinder"]
        lift_i = self.cams["intake"].valve_lift_crank(th, advance=adv)
        lift_e = self.cams["exhaust"].valve_lift_crank(th)
        Ai = effective_valve_area(lift_i, nv["intake"], val["intake_head_diameter"], gex["cd_intake"],
                                  gex["lift_to_diameter_cap"])
        Ae = effective_valve_area(lift_e, nv["exhaust"], val["exhaust_head_diameter"], gex["cd_exhaust"],
                                  gex["lift_to_diameter_cap"])
        self.tabs = (cylinder_volume(th, self.geom), dV_dtheta(th, self.geom), Ai, Ae)
        self.table_theta_deg = th / DEG
        self.lift_tables = (lift_i, lift_e)
        # seat-to-seat events in local crank degrees
        self.ivo_deg, self.ivc_deg = (a / DEG - self.vvt_advance_deg for a in self.cams["intake"].valve_open_crank())
        self.evo_deg, self.evc_deg = (a / DEG for a in self.cams["exhaust"].valve_open_crank())
        if not (0.0 < self.evo_deg < self.ivc_deg - CYCLE_DEG / 2 and self.ivc_deg < CYCLE_DEG):
            raise ValueError("valve events outside the layout assumed by the closed / open split")

    def _gas(self):
        gs = self.gspec["gas"]
        _, cpF, RF = frozen(float(gs["T_gamma_intake"]), 0.0, *self.props)
        _, cpP, RP = frozen(float(gs["T_gamma_exhaust"]), 1.0, *self.props)
        self.g_in, self.R_in = cpF / (cpF - RF), RF
        self.g_ex, self.R_ex = cpP / (cpP - RP), RP
        amb = self.gspec["ambient"]
        self.p_amb, self.T_amb = float(amb["p"]), float(amb["T"])
        self.rho_air_amb = self.p_amb / (self.tables.R_air * self.T_amb)
        self.c_intake = math.sqrt(self.g_in * self.R_in * self.T_amb)

    def _pipe(self, net, cfg, gamma, R, left, right, T, Y, dx=None):
        return net.add_pipe(cfg["length"], cfg["diameter"], gamma, R, left=left, right=right, dx=dx,
                            wall_T=cfg["wall_T"], friction_multiplier=cfg["friction_multiplier"],
                            heat_multiplier=cfg["heat_multiplier"], p=self.p_amb, T=T, Y=Y)

    def _build(self):
        g = self.gspec
        it, ex, num = g["intake"], g["exhaust"], g["numerics"]
        order = list(self.spec["engine"]["firing_order"])
        self.ncyl = len(order)
        net = Network(g, numerics_vector(g, self.tol, self.maxit))
        net.tabs, net.props = self.tabs, self.props
        A = Network.area
        # cylinders: volume index = cylinder number - 1
        self.phase_deg = np.array([FIRING_INTERVAL_DEG * order.index(c + 1) for c in range(self.ncyl)])
        for c in range(self.ncyl):
            net.add_cylinder(self.phase_deg[c] * DEG)
        gi, Ri, ge, Re = self.g_in, self.R_in, self.g_ex, self.R_ex
        # intake
        bst = self.boost
        if bst is None:
            self.v_amb_in = net.add_ambient(self.p_amb, self.T_amb, gi, Ri, Y=0.0)
            self.v_charge = None
        else:   # turbo: the inlet duct starts at the charge-air plenum (fed by gasdyn.elements.MapCompressor)
            self.v_amb_in = None
            self.v_charge = net.add_plenum(bst["charge_air_volume"], bst["p_charge_init"], self.T_amb, gi, Ri)
        v_src = self.v_amb_in if bst is None else self.v_charge
        self.v_plenum = net.add_plenum(it["plenum"]["volume"], self.p_amb, self.T_amb, gi, Ri)
        duct = it["inlet_duct"]
        self.throttle_slot_pipe = self._pipe(net, duct, gi, Ri, (v_src, duct["entry_cd"] * A(duct["diameter"])),
                                             (self.v_plenum, 0.0), self.T_amb, 0.0)
        run = it["runner"]
        self.runner_pipes = [self._pipe(net, run, gi, Ri, (self.v_plenum, run["entry_cd"] * A(run["diameter"])),
                                        (c, "intake"), self.T_amb, 0.0) for c in range(self.ncyl)]
        # exhaust
        Tx = float(ex["T_init"])
        self.v_amb_ex = net.add_ambient(self.p_amb, self.T_amb, ge, Re, Y=1.0)
        prim, cdj = ex["primary"], ex["junction_cd"]
        dxt = num["dx_tailpipe"]
        plen = lambda vol: net.add_plenum(vol, self.p_amb, Tx, ge, Re, Y=1.0)
        self.primary_pipes = []
        if ex["layout"] == "4-2-1":
            sec = ex["secondary"]
            v_col = plen(ex["collector_volume"])
            self.v_collector = v_col
            for pair in ((1, 4), (2, 3)):
                vj = plen(ex["junction_volume"])
                for cyl in pair:
                    self.primary_pipes.append(self._pipe(net, prim, ge, Re, (cyl - 1, "exhaust"),
                                                         (vj, cdj * A(prim["diameter"])), Tx, 1.0))
                self._pipe(net, sec, ge, Re, (vj, cdj * A(sec["diameter"])), (v_col, cdj * A(sec["diameter"])), Tx, 1.0)
        elif ex["layout"] == "4-1":
            v_col = plen(ex["junction_volume"])
            self.v_collector = v_col
            for c in range(self.ncyl):
                self.primary_pipes.append(self._pipe(net, prim, ge, Re, (c, "exhaust"),
                                                     (v_col, cdj * A(prim["diameter"])), Tx, 1.0))
        else:
            raise ValueError(ex["layout"])
        v_sil = plen(ex["silencer"]["volume"])
        tf, trr = ex["tailpipe_front"], ex["tailpipe_rear"]
        v_down = v_col
        self.turbine_pipe = self.wastegate_pipe = self.v_post_turbine = None
        if bst is not None:
            # collector -> turbine housing pipe -> turbine (nozzle-equivalent orifice) -> post-turbine volume
            # collector -> wastegate pipe -> wastegate (variable orifice)              -> post-turbine volume
            tp, wp = bst["turbine_pipe"], bst["wastegate_pipe"]
            self.v_post_turbine = plen(bst["post_turbine_volume"])
            self.turbine_pipe = self._pipe(net, tp, ge, Re, (v_col, cdj * A(tp["diameter"])),
                                           (self.v_post_turbine, bst["turbine_area"]), Tx, 1.0)
            self.wastegate_pipe = self._pipe(net, wp, ge, Re, (v_col, cdj * A(wp["diameter"])),
                                             (self.v_post_turbine, bst.get("wastegate_area", 0.0)), Tx, 1.0)
            v_down = self.v_post_turbine
        self._pipe(net, tf, ge, Re, (v_down, cdj * A(tf["diameter"])), (v_sil, cdj * A(tf["diameter"])), Tx, 1.0, dx=dxt)
        self._pipe(net, trr, ge, Re, (v_sil, cdj * A(trr["diameter"])), (self.v_amb_ex, ex["exit_cd"] * A(trr["diameter"])),
                   Tx, 1.0, dx=dxt)
        net.finalize()
        self.net = net
        self._cold = (net.U.copy(), net.S.copy(), net.Tg.copy())   # pristine pipes and plenums
        self.throttle_slot = int(net.e_slot[self.throttle_slot_pipe, 1])
        if bst is not None:
            self.turbine_slot = int(net.e_slot[self.turbine_pipe, 1])
            self.wastegate_slot = int(net.e_slot[self.wastegate_pipe, 1])
        self.set_throttle(float(it["throttle"]["wot_angle_deg"]))

    def set_throttle(self, angle_deg: float):
        th = self.gspec["intake"]["throttle"]
        self.throttle_deg = angle_deg
        self.net.sl_const[self.throttle_slot] = th["cd"] * throttle_area(angle_deg, th["bore"], th["closed_angle_deg"])

    # ------------------------------------------------------------------------------------------
    def _closed(self, c, m, x, T, rpm, st: EngineState):
        """Two-zone closed part from seat IVC to seat EVO for cylinder c; stores its EVO state."""
        comb = self.spec["combustion"]
        V_ivc = float(cylinder_volume(self.ivc_deg * DEG, self.geom))
        _, cpu, Ru = frozen(T, x, *self.props)
        p_ivc = m * Ru * T / V_ivc
        hp = heat_params(self.tspec, self.geom, rpm, "closed", on=True, ref=(p_ivc, V_ivc, T, cpu / (cpu - Ru)))
        step = float(self.gspec["numerics"]["closed_step_deg"])
        if self.closed_fn is not None:
            cl = self.closed_fn(m, x, T, rpm, hp, self.ivc_deg - CYCLE_DEG, self.evo_deg, step)
        else:
            cl = run_closed(self.tables, self.geom, m, x, T, self.ivc_deg - CYCLE_DEG, self.evo_deg, step, hp,
                            "wiebe", comb["heat_release_start_deg"], comb["burn_duration_deg"], comb["wiebe_a"],
                            comb["wiebe_m"], num=self.tspec["numerics"])
        x_open, U, T_open = merge_zones(cl, self.tables, self.tol, self.maxit)
        st.pend[c] = (m, m * x_open, U, T_open)
        st.w_closed[c] = cl.col(tz.O_W)[-1]
        st.q_closed[c] = cl.col(tz.O_QU)[-1] + cl.col(tz.O_QB)[-1]
        st.p_max[c] = float(np.max(cl.col(tz.O_P)))
        return cl

    def initial_state(self, rpm: float) -> EngineState:
        """Cold start: pipes at ambient pressure, cylinders filled according to their stroke at theta = 0."""
        net = self.net
        num = self.tspec["numerics"]
        st = EngineState(*(a.copy() for a in self._cold), np.zeros(len(net.v_kind), dtype=np.int64),
                         np.zeros((self.ncyl, 4)), np.zeros(self.ncyl), np.zeros(self.ncyl), np.zeros(self.ncyl))
        xr, T_ivc = float(num["xr_guess"]), float(num["T_ivc_guess"])
        Tx = float(self.gspec["exhaust"]["T_init"])
        V_ivc = float(cylinder_volume(self.ivc_deg * DEG, self.geom))
        _, _, Ru = frozen(T_ivc, xr, *self.props)
        m_guess = self.p_amb * V_ivc / (Ru * T_ivc)
        for c in range(self.ncyl):
            self._closed(c, m_guess, xr, T_ivc, rpm, st)
            local = (-self.phase_deg[c]) % CYCLE_DEG
            if self.evo_deg <= local < self.ivc_deg:        # gas exchange under way at theta = 0
                x, T = (1.0, Tx) if local < CYCLE_DEG / 2 else (xr, T_ivc)
                h, _, R = frozen(T, x, *self.props)
                m = self.p_amb * float(cylinder_volume(local * DEG, self.geom)) / (R * T)
                st.S[c, :] = 0.0
                st.S[c, fv.S_M], st.S[c, fv.S_MY], st.S[c, fv.S_E] = m, m * x, m * (h - R * T)
                st.Tg[c] = T
                st.act[c] = 1
        return st

    def _events(self):
        ev = []
        for c in range(self.ncyl):
            ev.append(((self.ivc_deg + self.phase_deg[c]) % CYCLE_DEG, c, "ivc"))
            ev.append(((self.evo_deg + self.phase_deg[c]) % CYCLE_DEG, c, "evo"))
        return sorted(ev)

    # ------------------------------------------------------------------------------------------
    def run(self, rpm: float, state: EngineState | None = None, throttle_deg: float | None = None,
            max_cycles: int | None = None, trace: bool = False, probe_deg: float | None = None) -> OperatingPoint:
        """Iterate engine cycles at constant speed until the mean net IMEP is periodic.

        probe_deg (M7, optional): also split the time marching every probe_deg crank degrees and sample the
        exhaust collector gas temperature there; the last cycle's mean is reported as T_exhaust_K. The
        default (None) marches event to event exactly as in M4."""
        cyc = self.gspec["cycle"]
        num = self.gspec["numerics"]
        net = self.net
        if throttle_deg is not None:
            self.set_throttle(throttle_deg)
        st = state.copy() if state is not None else self.initial_state(rpm)
        net.U, net.S, net.Tg = st.U, st.S, st.Tg
        net.v_act[:] = st.act
        net.hp = heat_params(self.tspec, self.geom, rpm, "gas_exchange", on=True)
        omega = rpm * 2.0 * math.pi / SEC_PER_MIN
        net.nu[fv.NU_DT_MAX] = float(num["max_step_deg"]) * DEG / omega
        net.nbad = net.nstep = 0
        Vd = self.geom.swept_volume
        yf = self.tables.fuel_mass_fraction
        events = self._events()
        tr_step = float(num["trace_step_deg"])
        nbin = int(round(CYCLE_DEG / tr_step))
        ntr = len(net.v_kind)
        max_cycles = int(max_cycles or cyc["max_cycles"])
        sub_deg = float(self.substep_deg) if (self.substep_fn is not None and self.substep_deg) else None
        rec = [dict() for _ in range(self.ncyl)]
        history = []
        ok_run = 0
        converged = False
        imep_prev = None
        tr = None
        for icyc in range(max_cycles):
            tr = np.full((nbin, ntr), np.nan) if trace else None
            theta = 0.0
            probes = []
            for th_ev, c, kind in events + [(CYCLE_DEG, -1, "end")]:
                while th_ev > theta:
                    nxt = th_ev
                    if probe_deg:
                        nxt = min(th_ev, (math.floor(theta / probe_deg + 1e-9) + 1.0) * probe_deg)
                    if sub_deg:
                        nxt = min(nxt, (math.floor(theta / sub_deg + 1e-9) + 1.0) * sub_deg)
                    dt = (nxt - theta) * DEG / omega
                    net.advance(dt, theta0=theta * DEG, omega=omega, trace=tr, trace_dth=tr_step * DEG)
                    theta = nxt
                    if self.substep_fn is not None:
                        self.substep_fn(net, theta, dt)
                    if probe_deg and abs(theta / probe_deg - round(theta / probe_deg)) < 1e-9:
                        probes.append(net.plenum_state(self.v_collector)[1])
                if kind == "ivc":
                    S = net.S[c]
                    m = S[fv.S_M]
                    x = min(max(S[fv.S_MY] / m, 0.0), 1.0)
                    T = T_from_u_frozen(S[fv.S_E] / m, x, net.Tg[c], *self.props, self.tol, self.maxit)
                    rec[c] = dict(w_net=st.w_closed[c] + S[fv.S_W], imep=(st.w_closed[c] + S[fv.S_W]) / Vd,
                                  w_open=S[fv.S_W], p_max=st.p_max[c], q_wall=st.q_closed[c] + S[fv.S_Q],
                                  m_trapped=m, xr=x, T_ivc=T, fresh_in=S[fv.S_FRESH], intake_net=S[fv.S_NET],
                                  exhaust_net=S[fv.S_EXH], m_evo=st.pend[c, 0])
                    self._closed(c, m, x, T, rpm, st)
                    net.v_act[c] = 0
                elif kind == "evo":
                    net.S[c, :] = 0.0
                    net.S[c, fv.S_M], net.S[c, fv.S_MY], net.S[c, fv.S_E] = st.pend[c, :3]
                    net.Tg[c] = st.pend[c, 3]
                    net.v_act[c] = 1
            if all(rec):
                imep = float(np.mean([r["imep"] for r in rec]))
                d = abs(imep / imep_prev - 1.0) if imep_prev else float("inf")
                history.append({"cycle": icyc, "imep_net_bar": imep / 1e5, "delta": d,
                                "m_trapped": float(np.mean([r["m_trapped"] for r in rec]))})
                ok_run = ok_run + 1 if d < float(cyc["imep_tol"]) else 0
                imep_prev = imep
                if self.cycle_fn is not None and self.cycle_fn(icyc, rec) is False:
                    ok_run = 0
                if ok_run >= int(cyc["consecutive"]) and icyc + 1 >= int(cyc["min_cycles"]):
                    converged = True
                    break
        st.act = net.v_act.copy()
        net.U, net.S, net.Tg = (a.copy() for a in self._cold)    # the returned state owns the final arrays

        # ---- cycle-averaged metrics (mean over the cylinders) ------------------------------------
        mean = lambda k: float(np.mean([r[k] for r in rec]))
        imep = mean("imep")
        p_max = mean("p_max")
        sp = 2.0 * self.geom.stroke * rpm / SEC_PER_MIN
        fmep = fmep_total(p_max, sp, self.gspec["friction"])
        bmep = imep - fmep
        disp = self.geom.displacement
        torque = bmep * disp / (2.0 * math.pi * REV_PER_CYCLE)
        power = torque * omega
        cyc_per_s = rpm / (SEC_PER_MIN * REV_PER_CYCLE)
        fresh = mean("fresh_in")
        fuel_rate = fresh * yf * self.ncyl * cyc_per_s
        air_rate = fresh * (1.0 - yf) * self.ncyl * cyc_per_s
        m_air_trapped = float(np.mean([r["m_trapped"] * (1.0 - r["xr"]) * (1.0 - yf) for r in rec]))
        lhv = float(self.spec["combustion"]["lhv"])
        metrics = {
            "rpm": rpm,
            "imep_net_bar": imep / 1e5,
            "fmep_bar": fmep / 1e5,
            "fmep_jets_bar": fmep_jets(sp, self.gspec["friction"]) / 1e5,
            "bmep_bar": bmep / 1e5,
            "torque_nm": torque,
            "power_kw": power / 1e3,
            "ve_trapped": m_air_trapped / (self.rho_air_amb * Vd),
            "ve_delivered": fresh * (1.0 - yf) / (self.rho_air_amb * Vd),
            "residual_fraction": mean("xr"),
            "bsfc_g_per_kwh": fuel_rate / power * 3.6e9 if power > 0.0 else float("nan"),
            "isfc_g_per_kwh": fuel_rate / (imep * disp * cyc_per_s) * 3.6e9,
            "eta_ind_net": imep * disp * cyc_per_s / (fuel_rate * lhv),
            "eta_brake": power / (fuel_rate * lhv),
            "air_kg_per_h": air_rate * 3600.0,
            "fuel_kg_per_h": fuel_rate * 3600.0,
            "p_max_bar": p_max / 1e5,
            "T_ivc_K": mean("T_ivc"),
            "mean_piston_speed": sp,
            "imep_cyl_spread": float(np.ptp([r["imep"] for r in rec]) / imep),
            "cycles_run": icyc + 1,
            "converged": converged,
            "imep_delta": history[-1]["delta"] if history else float("nan"),
            "floor_hits": int(net.nbad),
            "steps": int(net.nstep),
        }
        if probe_deg:
            metrics["T_exhaust_K"] = float(np.mean(probes))
        op = OperatingPoint(rpm, self.throttle_deg, converged, icyc + 1, metrics, rec, history, st)
        if trace and tr is not None:
            op.trace_deg = (np.arange(nbin) + 0.5) * tr_step
            op.trace = {"cylinder_p": tr[:, :self.ncyl], "plenum_p": tr[:, self.v_plenum]}
        return op
