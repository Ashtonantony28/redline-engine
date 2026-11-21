"""Speed-load maps (milestone M7): BSFC, VE, MBT / knock-limited spark, NOx, CA50, exhaust temperature,
brake torque.

    python -m sim.maps                 # default: M4 1D gas-dynamics engine, coarse grid (spec maps.rpm x map_levels)
    python -m sim.maps --quasi-steady  # fast option for dense grids: thermo.cycle filling/emptying gas exchange

writes out/sim/maps/<plant>.csv and contour PNGs; `--write-fixture` also refreshes fixtures/maps/maps_m4.csv
(a reviewed change only, CLAUDE.md). Each operating point is computed independently from a cold start:

1. load: the throttle plate angle is iterated (secant on effective area, warm-started inside the point)
   until the cycle-mean plenum pressure matches the map level (spec maps.map_tol); "wot" = throttle wide open.
   This step uses the spec Wiebe burn.
2. spark: on the cylinder-mean trapped state (mass, residual, T at IVC) the flame model (thermo.flame) gives
   MBT and the knock-limited advance (thermo.spark, thermo.knock); spark = min(MBT, KLSA - margin).
3. the point is re-run on the M4 engine with the flame model at that spark (warm start); the spark search is
   repeated `spark_reconverge_passes` times on the new trapped state and the engine re-run again.
4. knock, burn angles and emissions (thermo.emissions) are evaluated on the final trapped state and spark;
   torque, BSFC, VE and exhaust (collector) temperature come from the final M4 run.

The quasi-steady plant replaces the 1D network by constant intake (plenum = map level) and exhaust states
(spec maps.quasi_steady); it is ~an order of magnitude faster but has no wave dynamics, so its VE and
pumping work differ from the reported M4 maps. All numbers are MODELLED: there is no physical engine.
"""
from __future__ import annotations

import argparse
import copy
import csv
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from gasdyn.engine import EngineModel, throttle_area
from gasdyn.friction import fmep_total
from redline.spec import ROOT, geometry, load
from thermo import gasexchange as gx
from thermo.cycle import run_cycle
from thermo.emissions import Emissions, specific
from thermo.m7spec import combustion_spec
from thermo.props import build_tables, frozen, thermo_spec
from thermo.spark import SparkModel, Trapped
from gasdyn.spec import gasdyn_spec

OUT_DIR = ROOT / "out" / "sim" / "maps"
FIXTURE_DIR = ROOT / "fixtures" / "maps"
FIXTURE = FIXTURE_DIR / "maps_m4.csv"
PA_PER_BAR = 1.0e5
SEC_PER_MIN = 60.0
REV_PER_CYCLE = 2.0
G_PER_KG = 1.0e3
J_PER_KWH = 3.6e6
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"

COLUMNS = ["rpm", "level", "throttle_deg", "map_bar", "torque_nm", "bmep_bar", "imep_net_bar", "power_kw",
           "bsfc_g_per_kwh", "ve_trapped", "ve_delivered", "residual_fraction", "spark_deg", "mbt_deg", "klsa_deg",
           "klsa_status", "knock_limited", "ca10_deg", "ca50_deg", "ca90_deg", "knock_margin", "knock_onset_deg",
           "p_max_bar", "T_exhaust_K", "nox_g_per_kwh", "no_ppm", "co_g_per_kwh", "co_pct_evo", "co_pct_peak",
           "T_b_max_K", "converged", "cycles"]


def level_name(level) -> str:
    return "wot" if level == "wot" else f"{float(level) / PA_PER_BAR:.2f}bar"


def fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{float(v):.6g}"


# ---------------------------------------------------------------------------------------------
# shared post-processing
# ---------------------------------------------------------------------------------------------
class _Common:
    def __init__(self, model: str | None = None):
        self.spec = load()
        self.tspec = thermo_spec()
        self.cspec = combustion_spec()
        self.gspec = gasdyn_spec()
        self.tables = build_tables(self.tspec)
        self.geom = geometry(self.spec)
        self.sm = SparkModel(self.spec, self.tspec, self.cspec, self.gspec, self.tables, self.geom, model)
        self.em = Emissions(self.tspec, self.cspec, self.tables)
        self.ncyl = int(self.spec["engine"]["cylinders"])

    def finish(self, row: dict, st: Trapped, rpm: float, opt: dict, power_w: float):
        ev = self.sm.evaluate(st, rpm, opt["spark_deg"])
        e = self.em.closed(ev["cl"], rpm)
        row.update({"spark_deg": opt["spark_deg"], "mbt_deg": opt["mbt_deg"], "klsa_deg": opt["klsa_deg"],
                    "klsa_status": opt["klsa_status"], "knock_limited": opt["knock_limited"],
                    "ca10_deg": ev["ca10"], "ca50_deg": ev["ca50"], "ca90_deg": ev["ca90"],
                    "knock_margin": ev["knock"]["margin"], "knock_onset_deg": ev["knock"]["onset_deg"],
                    "nox_g_per_kwh": specific(e["m_NOx_kg"], self.ncyl, rpm, power_w), "no_ppm": e["NO_ppm"],
                    "co_g_per_kwh": specific(e["m_CO_kg"], self.ncyl, rpm, power_w), "co_pct_evo": e["CO_pct_evo"],
                    "co_pct_peak": e["CO_pct_peak"], "T_b_max_K": e["T_b_max_K"]})
        return row


# ---------------------------------------------------------------------------------------------
# M4 plant (default)
# ---------------------------------------------------------------------------------------------
class M4Plant(_Common):
    def __init__(self, model: str | None = None):
        super().__init__(model)
        self.eng = EngineModel(spec=self.spec, tspec=self.tspec, gspec=self.gspec, tables=self.tables)
        self.mp = self.cspec["maps"]
        th = self.gspec["intake"]["throttle"]
        self.th = th
        self.cd_thr = float(th["cd"])
        self.a_wot = float(th["wot_angle_deg"])
        self.a_min = float(self.mp["throttle_min_deg"])
        self.probe = float(self.mp["probe_deg"])

    def _area(self, angle):
        return throttle_area(angle, self.th["bore"], self.th["closed_angle_deg"])

    def _angle(self, area):
        lo, hi = self.a_min, self.a_wot
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self._area(mid) < area:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def _flow_per_area(self, p_down):
        e = self.eng
        return gx.orifice(self.cd_thr, e.p_amb, e.T_amb, e.R_in, e.g_in, p_down,
                          float(self.tspec["gas_exchange"]["linear_flow_pr"]))

    def _run(self, rpm, angle, state, closed_fn):
        self.eng.closed_fn = closed_fn
        op = self.eng.run(rpm, state=state, throttle_deg=angle, trace=True, probe_deg=self.probe)
        op.metrics["map_pa"] = float(np.nanmean(op.trace["plenum_p"]))
        return op

    @staticmethod
    def trapped(op, eng) -> Trapped:
        rec = op.per_cylinder
        mean = lambda k: float(np.mean([r[k] for r in rec]))
        return Trapped(mean("m_trapped"), mean("xr"), mean("T_ivc"), eng.ivc_deg, eng.evo_deg,
                       float(eng.gspec["numerics"]["closed_step_deg"]))

    def load_point(self, rpm, level):
        """Throttle angle and converged Wiebe operating point at the map level (cold start)."""
        if level == "wot":
            return self._run(rpm, self.a_wot, None, None), True
        target = float(level)
        tol = float(self.mp["map_tol"])
        e = self.eng
        # first guess: delivered air at a nominal VE through the throttle orifice from ambient
        vd = self.geom.displacement
        mdot = float(self.mp["ve_guess"]) * target / (e.R_in * e.T_amb) * vd * rpm / SEC_PER_MIN / REV_PER_CYCLE
        A = min(max(mdot / self._flow_per_area(target), self._area(self.a_min)), self._area(self.a_wot))
        op = self._run(rpm, self._angle(A), None, None)
        hist = [(A, op.metrics["map_pa"])]
        xlim = float(self.mp["dilution_limit"])
        for _ in range(int(self.mp["throttle_max_iter"])):
            A0, p0 = hist[-1]
            if abs(p0 / target - 1.0) < tol:
                return op, True
            if p0 > target and op.metrics["residual_fraction"] > xlim:
                # residual only grows as the plenum pressure falls further: beyond the dilution limit
                return op, True
            if len(hist) == 1:
                # same delivered-air model, corrected by the achieved pressure
                A1 = A0 * (target / p0) * self._flow_per_area(p0) / self._flow_per_area(target)
            else:
                Am, pm = hist[-2]
                A1 = A0 + (target - p0) * (A0 - Am) / (p0 - pm) if p0 != pm else A0
            A1 = min(max(A1, self._area(self.a_min)), self._area(self.a_wot))
            if A1 == A0:
                break
            op = self._run(rpm, self._angle(A1), op.state, None)
            hist.append((A1, op.metrics["map_pa"]))
        A0, p0 = hist[-1]
        return op, abs(p0 / target - 1.0) < tol

    def point(self, rpm, level):
        op, reached = self.load_point(rpm, level)
        if not reached:
            return {"skip": "map level above the wide-open-throttle plenum pressure"}
        if op.metrics["residual_fraction"] > float(self.mp["dilution_limit"]):
            return {"skip": f"residual fraction {op.metrics['residual_fraction']:.3f} beyond the dilution limit"}
        angle = op.throttle_deg
        st = self.trapped(op, self.eng)
        opt = self.sm.optimise(st, rpm)
        for _ in range(int(self.mp["spark_reconverge_passes"]) + 1):
            op = self._run(rpm, angle, op.state, self.sm.closed_fn_engine(opt["spark_deg"]))
            st = self.trapped(op, self.eng)
            if _ < int(self.mp["spark_reconverge_passes"]):
                opt = self.sm.optimise(st, rpm)
        m = op.metrics
        row = {"rpm": rpm, "level": level_name(level), "throttle_deg": angle, "map_bar": m["map_pa"] / PA_PER_BAR,
               "torque_nm": m["torque_nm"], "bmep_bar": m["bmep_bar"], "imep_net_bar": m["imep_net_bar"],
               "power_kw": m["power_kw"], "bsfc_g_per_kwh": m["bsfc_g_per_kwh"], "ve_trapped": m["ve_trapped"],
               "ve_delivered": m["ve_delivered"], "residual_fraction": m["residual_fraction"],
               "p_max_bar": m["p_max_bar"], "T_exhaust_K": m["T_exhaust_K"], "converged": m["converged"],
               "cycles": m["cycles_run"]}
        self.eng.closed_fn = None
        return self.finish(row, st, rpm, opt, m["power_kw"] * 1.0e3)


# ---------------------------------------------------------------------------------------------
# quasi-steady plant (fast option)
# ---------------------------------------------------------------------------------------------
class QuasiSteadyPlant(_Common):
    LOAD = "m7_map"

    def __init__(self, model: str | None = None):
        super().__init__(model)
        self.qs = self.cspec["maps"]["quasi_steady"]
        self.step = float(self.qs["step_deg"])
        val = self.spec["valves"]
        self.ivc, self.evo = float(val["ivc_deg"]), float(val["evo_deg"])

    def _tspec(self, p_map):
        t = copy.deepcopy(self.tspec)
        t["loads"][self.LOAD] = {"p_intake": float(p_map), "T_intake": float(self.qs["T_intake"]),
                                 "p_exhaust": float(self.qs["p_exhaust"]), "T_exhaust": float(self.qs["T_exhaust"])}
        return t

    def _cycle(self, rpm, p_map, closed_fn):
        return run_cycle(rpm, self.LOAD, self.step, tables=self.tables, spec=self.spec, tspec=self._tspec(p_map),
                         closed_fn=closed_fn)

    def point(self, rpm, level):
        p_map = float(level)
        res = self._cycle(rpm, p_map, None)
        mk = lambda r: Trapped(r.metrics["m_trapped_kg"], r.metrics["residual_fraction"], r.metrics["T_ivc_K"],
                               self.ivc, self.evo, self.step)
        st = mk(res)
        opt = self.sm.optimise(st, rpm)
        passes = int(self.cspec["maps"]["spark_reconverge_passes"])
        for i in range(passes + 1):
            res = self._cycle(rpm, p_map, self.sm.closed_fn_cycle(opt["spark_deg"], rpm))
            st = mk(res)
            if i < passes:
                opt = self.sm.optimise(st, rpm)
        m = res.metrics
        sp = 2.0 * self.geom.stroke * rpm / SEC_PER_MIN
        fmep = fmep_total(m["p_max_bar"] * PA_PER_BAR, sp, self.gspec["friction"])
        bmep = m["imep_net_bar"] * PA_PER_BAR - fmep
        vd = self.geom.swept_volume
        torque = bmep * self.geom.displacement / (2.0 * math.pi * REV_PER_CYCLE)
        power = torque * rpm * 2.0 * math.pi / SEC_PER_MIN
        oo = res.open_out
        out = np.clip(oo[:, gx.GO_DME], 0.0, None)
        T_exh = float(np.sum(out * oo[:, gx.GO_T]) / np.sum(out))
        amb = self.gspec["ambient"]
        rho_amb = float(amb["p"]) / (self.tables.R_air * float(amb["T"]))
        m_air_in = m["fresh_in_kg"] * (1.0 - self.tables.fuel_mass_fraction)
        m_air_trapped = m["m_trapped_kg"] * (1.0 - m["residual_fraction"]) * (1.0 - self.tables.fuel_mass_fraction)
        row = {"rpm": rpm, "level": level_name(level), "throttle_deg": float("nan"), "map_bar": p_map / PA_PER_BAR,
               "torque_nm": torque, "bmep_bar": bmep / PA_PER_BAR, "imep_net_bar": m["imep_net_bar"],
               "power_kw": power / 1.0e3,
               "bsfc_g_per_kwh": m["m_fuel_kg"] * G_PER_KG / (bmep * vd / J_PER_KWH) if bmep > 0.0 else float("nan"),
               "ve_trapped": m_air_trapped / (rho_amb * vd), "ve_delivered": m_air_in / (rho_amb * vd),
               "residual_fraction": m["residual_fraction"], "p_max_bar": m["p_max_bar"], "T_exhaust_K": T_exh,
               "converged": m["converged"], "cycles": m["cycles"]}
        return self.finish(row, st, rpm, opt, power)


# ---------------------------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------------------------
_PLANT = None


def _task(args):
    kind, rpm, level = args
    global _PLANT
    if _PLANT is None or _PLANT.__class__.__name__ != kind:
        _PLANT = M4Plant() if kind == "M4Plant" else QuasiSteadyPlant()
    return _PLANT.point(float(rpm), level)


def grid(cspec: dict, quasi_steady: bool):
    mp = cspec["maps"]
    if quasi_steady:
        a, b, s = mp["quasi_steady"]["rpm"]
        rpms = [float(a + i * s) for i in range(int(round((b - a) / s)) + 1)]
        return [(r, float(l)) for r in rpms for l in mp["quasi_steady"]["map_levels"]]
    return [(float(r), l if l == "wot" else float(l)) for r in mp["rpm"] for l in mp["map_levels"]]


def run_points(points, quasi_steady=False, jobs=1):
    kind = "QuasiSteadyPlant" if quasi_steady else "M4Plant"
    tasks = [(kind, r, l) for r, l in points]
    if jobs <= 1 or len(tasks) == 1:
        rows = [_task(t) for t in tasks]
    else:
        from sim.wot import _warm_kernels
        _warm_kernels()
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            rows = list(pool.map(_task, tasks))
    return ([r for r in rows if "skip" not in r],
            [(p, r["skip"]) for p, r in zip(points, rows) if "skip" in r])


def csv_rows(rows):
    return [[fmt(r[c]) for c in COLUMNS] for r in rows]


def write_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(csv_rows(rows))


def read_csv(path: Path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


PANELS = [("bsfc_g_per_kwh", "BSFC [g/kWh]", "viridis_r"), ("ve_trapped", "Volumetric efficiency (trapped) [-]", "viridis"),
          ("spark_deg", "Spark advance [deg BTDC]", "viridis"), ("nox_g_per_kwh", "NOx as NO2 [g/kWh]", "magma_r"),
          ("ca50_deg", "CA50 [deg ATDC]", "viridis"), ("T_exhaust_K", "Exhaust temperature [K]", "inferno"),
          ("torque_nm", "Brake torque [Nm]", "viridis")]


def plot(rows, out: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    rpm = np.array([r["rpm"] for r in rows], dtype=float)
    bmep = np.array([r["bmep_bar"] for r in rows], dtype=float)
    kl = np.array([bool(r["knock_limited"]) for r in rows])
    ok = np.isfinite(bmep) & (bmep > 0.0)
    tri = mtri.Triangulation(rpm[ok], bmep[ok])
    for key, label, cmap in PANELS:
        z = np.array([r[key] for r in rows], dtype=float)[ok]
        fig, ax = plt.subplots(figsize=(7.5, 5), facecolor=SURFACE)
        ax.set_facecolor(SURFACE)
        good = np.isfinite(z)
        if key == "bsfc_g_per_kwh":
            lv = np.percentile(z[good], np.linspace(0, 90, 10))
            cs = ax.tricontourf(tri, np.where(good, z, np.nanmax(z)), levels=np.unique(lv), cmap=cmap, extend="max")
        else:
            cs = ax.tricontourf(tri, np.where(good, z, np.nanmean(z)), levels=10, cmap=cmap)
        ax.tricontour(tri, np.where(good, z, np.nanmean(z)), levels=cs.levels, colors=INK_2, linewidths=0.4)
        ax.plot(rpm[ok & ~kl], bmep[ok & ~kl], "o", ms=3, mfc="white", mec=INK, mew=0.6, label="MBT")
        ax.plot(rpm[ok & kl], bmep[ok & kl], "x", ms=5, color=INK, label="knock-limited")
        fig.colorbar(cs, ax=ax, label=label)
        ax.set_xlabel("Engine speed [rpm]", color=INK_2)
        ax.set_ylabel("BMEP [bar]", color=INK_2)
        ax.set_title(f"{label}: {title}\nMODELLED, no physical engine", color=INK, loc="left", fontsize=9)
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        for s in ax.spines.values():
            s.set_color(GRID)
        fig.tight_layout()
        fig.savefig(out / f"{key}.png", dpi=120, facecolor=SURFACE)
        plt.close(fig)


def summary(rows) -> dict:
    ok = [r for r in rows if np.isfinite(r["bsfc_g_per_kwh"]) and r["bsfc_g_per_kwh"] > 0.0]
    best = min(ok, key=lambda r: r["bsfc_g_per_kwh"])
    kl = [(r["rpm"], r["level"]) for r in rows if r["knock_limited"]]
    nox = [r["nox_g_per_kwh"] for r in ok]
    return {"bsfc_min": best["bsfc_g_per_kwh"], "bsfc_min_rpm": best["rpm"], "bsfc_min_level": best["level"],
            "bsfc_min_bmep": best["bmep_bar"], "knock_limited": kl, "nox_min": min(nox), "nox_max": max(nox)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quasi-steady", action="store_true", help="fast quasi-steady thermo gas exchange (dense grid)")
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--out", default=str(OUT_DIR))
    ap.add_argument("--write-fixture", action="store_true", help="refresh fixtures/maps/maps_m4.csv (reviewed change)")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args(argv)
    cspec = combustion_spec()
    pts = grid(cspec, a.quasi_steady)
    rows, skipped = run_points(pts, a.quasi_steady, a.jobs)
    tag = "quasi_steady" if a.quasi_steady else "m4"
    out = Path(a.out) / tag
    write_csv(rows, out / f"maps_{tag}.csv")
    if not a.no_plot:
        plot(rows, out, "M4 1D engine" if tag == "m4" else "quasi-steady gas exchange")
    if a.write_fixture and tag == "m4":
        write_csv(rows, FIXTURE)
    print(f"MODELLED speed-load maps ({tag}; no physical engine)")
    print(f"{'rpm':>5} {'level':>8} {'MAP':>5} {'BMEP':>5} {'T[Nm]':>6} {'BSFC':>6} {'VE':>5} {'spark':>6} {'MBT':>5} "
          f"{'KLSA':>5} {'CA50':>5} {'NOx':>5} {'Texh':>5}")
    for r in rows:
        print(f"{r['rpm']:5.0f} {r['level']:>8} {r['map_bar']:5.2f} {r['bmep_bar']:5.2f} {r['torque_nm']:6.1f} "
              f"{r['bsfc_g_per_kwh']:6.1f} {r['ve_trapped']:5.3f} {r['spark_deg']:6.1f} {r['mbt_deg']:5.1f} "
              f"{r['klsa_deg']:5.1f}{'*' if r['knock_limited'] else ' '}{r['ca50_deg']:5.1f} {r['nox_g_per_kwh']:5.1f} "
              f"{r['T_exhaust_K']:5.0f}")
    with open(out / f"skipped_{tag}.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["rpm", "level", "reason"])
        for (r, l), why in skipped:
            w.writerow([fmt(r), level_name(l), why])
            print(f"skipped {r:.0f} rpm {level_name(l)}: {why}")
    s = summary(rows)
    print(f"min BSFC {s['bsfc_min']:.1f} g/kWh at {s['bsfc_min_rpm']:.0f} rpm, {s['bsfc_min_level']} "
          f"(BMEP {s['bsfc_min_bmep']:.2f} bar); knock-limited points: {s['knock_limited']}; "
          f"NOx {s['nox_min']:.1f}..{s['nox_max']:.1f} g/kWh")
    print(f"wrote {out}" + (f" and {FIXTURE}" if a.write_fixture and tag == "m4" else ""))
    return rows


if __name__ == "__main__":
    main()
