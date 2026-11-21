"""Load cases for the M6 FEA, taken from the existing models (nothing invented here):

- Cylinder pressure p(theta): thermo.cycle.run_cycle (M2 two-zone model) at WOT over the rpm sweep of
  spec/fea.yaml; the worst speed is picked per part (largest rod compression, largest crankpin fatigue
  demand, largest main-bearing load, largest peak pressure / crown heat flux).
- Forces: dyn.forces / dyn.bearings (M5) driven by those traces, with the reciprocating and rotating
  masses taken from the CAD solids (cad/ M1) instead of the spec budget.
- Piston crown heat flux: Woschni (thermo.woschni) evaluated along the converged cycle at the spec
  piston wall temperature (thermo.yaml heat_transfer.T_piston), cycle averaged.
Crank angles follow ADR 0001 (0 = firing TDC of the cylinder, 720-deg cycle).
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from cad.conrod import conrod
from cad.massprops import mass_props
from cad.params import Params
from cad.piston import piston_parts
from dyn import params as dparams
from dyn.bearings import big_end_loads, main_bearing_loads
from dyn.forces import cylinder_forces
from dyn.pressure import PressureTrace
from fea.config import fea_spec
from thermo.cycle import CYCLE_DEG, run_cycle
from thermo.props import build_tables, frozen, thermo_spec
from thermo.woschni import coefficient, heat_params
from redline.spec import geometry

DEG = math.pi / 180.0


# ---------------------------------------------------------------------------------------------
# masses from CAD
# ---------------------------------------------------------------------------------------------
def cad_masses(p: Params) -> dict:
    """Piston assembly and conrod masses and the conrod small-end share from the CAD solids."""
    pa = sum(mass_props(s, p.density(m), p.k).mass for s, m in piston_parts(p).values())
    rod = mass_props(conrod(p), p.density("conrod"), p.k)
    frac = float(rod.com[2] / p.geom.conrod_length)
    return {"piston_assembly": pa, "conrod": rod.mass, "conrod_small_end_fraction": frac,
            "reciprocating": pa + frac * rod.mass}


def spec_with_cad_masses(p: Params, masses: dict | None = None) -> dict:
    spec = copy.deepcopy(p.engine)
    m = masses or cad_masses(p)
    for k in ("piston_assembly", "conrod", "conrod_small_end_fraction"):
        spec["masses"][k] = m[k]
    return spec


# ---------------------------------------------------------------------------------------------
# thermo sweep
# ---------------------------------------------------------------------------------------------
@dataclass
class WotPoint:
    rpm: float
    trace: PressureTrace
    p_max: float
    theta_p_max_deg: float
    q_piston_mean: float          # W/m^2 on the piston (per wetted area) at the spec wall temperature
    h_piston_mean: float          # W/(m^2 K)
    metrics: dict = field(default_factory=dict)


def _woschni_flux(res, tspec, spec, tables):
    """Cycle-mean Woschni h and q onto the piston at T_piston (time average = crank-angle average)."""
    geom = geometry(spec)
    ht = tspec["heat_transfer"]
    num_cl = len(res.closed.theta)
    m = res.metrics
    xr = m["residual_fraction"]
    _hu, cpu, Ru = frozen(m["T_ivc_K"], xr, *tables.args()[:5])
    V_ivc = res.closed.V[0]
    hp_c = heat_params(tspec, geom, res.rpm, "closed", ref=(res.p[0], V_ivc, m["T_ivc_K"], cpu / (cpu - Ru)))
    hp_o = heat_params(tspec, geom, res.rpm, "gas_exchange")
    soc = float(spec["combustion"]["heat_release_start_deg"])
    h = np.empty(len(res.theta))
    for i, (th, p, T, V) in enumerate(zip(res.theta, res.p, res.T_bulk, res.V)):
        if i < num_cl:
            h[i] = coefficient(p, T, V, th >= soc, hp_c)
        else:
            h[i] = coefficient(p, T, V, False, hp_o)
    q = h * (res.T_bulk - ht["T_piston"])
    # uniform grid (closed and open parts share the step), duplicate EVO point already dropped
    return float(h.mean()), float(q.mean())


def _trace(res, step_deg: float) -> PressureTrace:
    th = np.mod(res.theta, CYCLE_DEG)
    order = np.argsort(th, kind="stable")
    th, p = th[order], res.p[order]
    grid = np.arange(0.0, CYCLE_DEG + 0.5 * step_deg, step_deg)
    pp = np.interp(grid, th, p, period=CYCLE_DEG)
    return PressureTrace.from_degrees(grid, pp, source=f"thermo.run_cycle {res.rpm:.0f} rpm {res.load}",
                                      meta={"rpm": res.rpm})


@lru_cache(maxsize=4)
def _wot_sweep_cached(spec_key: str, rpms: tuple, load_name: str, step_deg: float):
    tspec = thermo_spec()
    tables = build_tables(tspec)
    out = []
    for rpm in rpms:
        res = run_cycle(rpm, load_name, tables=tables)
        hm, qm = _woschni_flux(res, tspec, res_spec(), tables)
        m = res.metrics
        out.append(WotPoint(rpm, _trace(res, step_deg), m["p_max_bar"] * 1.0e5, m["theta_p_max_deg"], qm, hm,
                            {k: m[k] for k in ("imep_net_bar", "p_max_bar", "theta_p_max_deg", "converged")}))
    return out


def res_spec():
    from redline.spec import load
    return load()


def rpm_sweep(fspec: dict | None = None) -> tuple:
    s = (fspec or fea_spec())["loads"]["rpm_sweep"]
    n = int(round((s["max"] - s["min"]) / s["step"])) + 1
    return tuple(float(v) for v in np.linspace(s["min"], s["max"], n))


def wot_sweep(fspec: dict | None = None) -> list[WotPoint]:
    fspec = fspec or fea_spec()
    ld = fspec["loads"]
    return _wot_sweep_cached("spec", rpm_sweep(fspec), ld["load"], float(ld["theta_step_deg"]))


def theta_grid(fspec: dict | None = None) -> np.ndarray:
    step = float((fspec or fea_spec())["loads"]["theta_step_deg"])
    return np.arange(0.0, CYCLE_DEG, step) * DEG


def overrun_trace(fspec=None, dspec=None) -> PressureTrace:
    """Closed-throttle overrun: cylinder pressure = crankcase pressure (zero net gas force)."""
    dspec = dspec or dparams.dyn_spec()
    p_cc = dspec["layout"]["crankcase_pressure"]
    grid = np.array([0.0, CYCLE_DEG])
    return PressureTrace.from_degrees(grid, np.full(2, p_cc), source="overrun: p_cyl = p_crankcase")


# ---------------------------------------------------------------------------------------------
# conrod
# ---------------------------------------------------------------------------------------------
def conrod_loads(p: Params, fspec=None) -> dict:
    """Max rod compression (firing, worst WOT speed) and max tension (redline overrun TDC) [N]."""
    fspec = fspec or fea_spec()
    masses = cad_masses(p)
    spec = spec_with_cad_masses(p, masses)
    dspec = dparams.dyn_spec()
    th = theta_grid(fspec)
    best = None
    per_rpm = []
    for pt in wot_sweep(fspec):
        f = cylinder_forces(th, pt.trace, dparams.omega_from_rpm(pt.rpm), spec, dspec)
        i = int(np.argmax(f.F_rod))
        per_rpm.append({"rpm": pt.rpm, "F_rod_max_N": float(f.F_rod[i]), "theta_deg": float(th[i] / DEG),
                        "F_gas_N": float(f.F_gas[i]), "F_inertia_N": float(f.F_inertia[i]),
                        "F_rod_min_N": float(f.F_rod.min())})
        if best is None or f.F_rod[i] > best["F_rod_max_N"]:
            best = per_rpm[-1]
    red = float(p.engine["engine"]["redline_rpm"])
    fo = cylinder_forces(th, overrun_trace(fspec, dspec), dparams.omega_from_rpm(red), spec, dspec)
    j = int(np.argmin(fo.F_rod))
    wot_red = [r for r in per_rpm if r["rpm"] == red]
    return {"compression_N": best["F_rod_max_N"], "compression_case": best,
            "tension_N": float(-fo.F_rod[j]), "tension_theta_deg": float(th[j] / DEG), "redline_rpm": red,
            "wot_redline_min_F_rod_N": wot_red[0]["F_rod_min_N"] if wot_red else None,
            "per_rpm": per_rpm, "masses": masses}


# ---------------------------------------------------------------------------------------------
# crankshaft
# ---------------------------------------------------------------------------------------------
def crank_histories(p: Params, fspec=None) -> list[dict]:
    """Per WOT speed: for each throw k (front to rear) the crankpin force in the rotating throw frame
    (radial outward, tangential in the rotation direction) and the torque entering the throw from the
    front (sum of the torques of the cylinders ahead of it), over one 720-deg cycle of cylinder 1."""
    fspec = fspec or fea_spec()
    spec = spec_with_cad_masses(p)
    dspec = dparams.dyn_spec()
    th = theta_grid(fspec)
    phases = dparams.firing_offsets(spec)
    out = []
    for pt in wot_sweep(fspec):
        w = dparams.omega_from_rpm(pt.rpm)
        Fr, Ft, T = [], [], []
        for ph in phases:
            t_loc = th - ph
            be = big_end_loads(t_loc, pt.trace, w, spec, dspec)
            Fr.append(be.pin_frame[0])
            Ft.append(be.pin_frame[1])
            T.append(cylinder_forces(t_loc, pt.trace, w, spec, dspec).T)
        T = np.array(T)
        T_in = np.vstack([np.zeros(th.size), np.cumsum(T, axis=0)[:-1]])
        out.append({"rpm": pt.rpm, "omega": w, "theta": th, "Fr": np.array(Fr), "Ft": np.array(Ft),
                    "T_in": T_in, "T_cyl": T})
    return out


# ---------------------------------------------------------------------------------------------
# main bearings
# ---------------------------------------------------------------------------------------------
def main_bearing_peak(p: Params, fspec=None) -> dict:
    """Largest main-bearing load directed into the cap (crank pushing down, -z), over all mains and speeds."""
    fspec = fspec or fea_spec()
    spec = spec_with_cad_masses(p)
    dspec = dparams.dyn_spec()
    th = theta_grid(fspec)
    best = None
    for pt in wot_sweep(fspec):
        L = main_bearing_loads(th, pt.trace, dparams.omega_from_rpm(pt.rpm), spec, dspec)  # (n,2,N) (z,y)
        down = -L[:, 0, :]
        k, i = np.unravel_index(np.argmax(down), down.shape)
        if best is None or down[k, i] > best["F_down_N"]:
            Fz, Fy = L[k, 0, i], L[k, 1, i]
            best = {"rpm": pt.rpm, "main": int(k), "theta_deg": float(th[i] / DEG), "F_down_N": float(down[k, i]),
                    "Fz_N": float(Fz), "Fy_N": float(Fy), "F_mag_N": float(math.hypot(Fz, Fy))}
    return best


# ---------------------------------------------------------------------------------------------
# piston
# ---------------------------------------------------------------------------------------------
def piston_loads(fspec=None) -> dict:
    fspec = fspec or fea_spec()
    pts = wot_sweep(fspec)
    ip = int(np.argmax([pt.p_max for pt in pts]))
    iq = int(np.argmax([pt.q_piston_mean for pt in pts]))
    tspec = thermo_spec()
    return {"p_max_Pa": pts[ip].p_max, "p_max_rpm": pts[ip].rpm, "p_max_theta_deg": pts[ip].theta_p_max_deg,
            "q_mean_Wm2": pts[iq].q_piston_mean, "h_mean_Wm2K": pts[iq].h_piston_mean, "q_rpm": pts[iq].rpm,
            "T_piston_spec_K": tspec["heat_transfer"]["T_piston"],
            "T_liner_K": tspec["heat_transfer"]["T_liner"],
            "piston_area_factor": tspec["heat_transfer"]["piston_area_factor"],
            "crankcase_pressure_Pa": dparams.dyn_spec()["layout"]["crankcase_pressure"],
            "sweep": [{"rpm": pt.rpm, "p_max_bar": pt.p_max / 1.0e5, "q_mean_Wm2": pt.q_piston_mean} for pt in pts]}
