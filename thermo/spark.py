"""Spark timing (milestone M7): MBT search capped by the knock-limited spark advance (KLSA).

Spark is optimised on the closed two-zone part for a fixed trapped state (mass, residual fraction and
temperature at IVC, from the M4 1D engine or the quasi-steady thermo cycle). The objective is the closed-
part work minus the Chen-Flynn friction work, which depends on the peak pressure (gasdyn/friction.py);
gas-exchange work does not depend on spark for a fixed trapped state, so this maximises brake torque.

    MBT    = argmax_adv [W_closed(adv) - FMEP(p_max(adv), S_p) V_d]   coarse scan + golden section
    KLSA   = max adv with knock margin < 1                            bisection (thermo.knock)
    spark  = min(MBT, KLSA - klsa_margin_deg)

Spark advance is in deg BTDC (positive before firing TDC, ADR 0001 angle = -advance). MODELLED only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from gasdyn.friction import fmep_total
from kin.slider_crank import cylinder_volume
from thermo import twozone as tz
from thermo.cycle import CYCLE_DEG, run_closed
from thermo.flame import FlameModel
from thermo.knock import klsa_search, knock_analysis
from thermo.props import frozen
from thermo.woschni import heat_params

GOLDEN = (math.sqrt(5.0) - 1.0) / 2.0
SEC_PER_MIN = 60.0


@dataclass(frozen=True)
class Trapped:
    """Trapped state at IVC and the closed-part window [deg, ADR 0001 local angles]."""
    m: float
    xr: float
    T: float
    ivc_deg: float
    evo_deg: float
    step_deg: float


def burn_points(cl):
    xb = cl.col(tz.O_XB)
    ok = xb < 1.0 - 1.0e-4
    f = lambda v: float(np.interp(v, xb[ok], cl.theta[ok])) if xb[-1] > v else float("nan")
    return f(0.1), f(0.5), f(0.9)


class SparkModel:
    def __init__(self, spec, tspec, cspec, gspec, tables, geom, model: str | None = None):
        self.spec, self.tspec, self.cspec, self.gspec, self.tables, self.geom = spec, tspec, cspec, gspec, tables, geom
        self.model = model or cspec["flame"]["model"]
        self.flame = FlameModel(spec, tspec, cspec, gspec, tables, geom)
        self.kspec = cspec["knock"]
        self.sspec = cspec["spark"]
        comb = spec["combustion"]
        self.ign_delay_deg = float(comb["spark_advance_deg"]) + float(comb["heat_release_start_deg"])
        self.num = tspec["numerics"]

    # ------------------------------------------------------------------------------------------
    def closed(self, m, xr, T, rpm, hp, theta_s, theta_e, step_deg, adv):
        """Closed part with spark advance adv [deg BTDC] (flame model or spark-shifted spec Wiebe)."""
        if self.model == "flame":
            cl, _ = self.flame.closed(m, xr, T, theta_s, theta_e, step_deg, hp, rpm, adv, self.num)
            return cl
        if self.model == "wiebe":
            comb = self.spec["combustion"]
            return run_closed(self.tables, self.geom, m, xr, T, theta_s, theta_e, step_deg, hp, "wiebe",
                              -adv + self.ign_delay_deg, comb["burn_duration_deg"], comb["wiebe_a"], comb["wiebe_m"],
                              num=self.num)
        raise ValueError(self.model)

    def closed_fn_engine(self, adv):
        """Hook for gasdyn.engine.EngineModel.closed_fn at a fixed spark advance."""
        def fn(m, x, T, rpm, hp, theta_s, theta_e, step):
            return self.closed(m, x, T, rpm, hp, theta_s, theta_e, step, adv)
        return fn

    def closed_fn_cycle(self, adv, rpm):
        """Hook for thermo.cycle.run_cycle(closed_fn=...) at a fixed spark advance."""
        def fn(tables, geom, m, xr, T, theta_s, theta_e, step, hp, num):
            return self.closed(m, xr, T, rpm, hp, theta_s, theta_e, step, adv)
        return fn

    def hp_closed(self, st: Trapped, rpm):
        V_ivc = float(cylinder_volume(math.radians(st.ivc_deg), self.geom))
        _, cpu, Ru = frozen(st.T, st.xr, *self.tables.args()[:5])
        p_ivc = st.m * Ru * st.T / V_ivc
        return heat_params(self.tspec, self.geom, rpm, "closed", on=True, ref=(p_ivc, V_ivc, st.T, cpu / (cpu - Ru)))

    # ------------------------------------------------------------------------------------------
    def evaluate(self, st: Trapped, rpm, adv, hp=None):
        hp = self.hp_closed(st, rpm) if hp is None else hp
        cl = self.closed(st.m, st.xr, st.T, rpm, hp, st.ivc_deg - CYCLE_DEG, st.evo_deg, st.step_deg, adv)
        p_max = float(np.max(cl.col(tz.O_P)))
        sp = 2.0 * self.geom.stroke * rpm / SEC_PER_MIN
        W = float(cl.col(tz.O_W)[-1])
        obj = W - fmep_total(p_max, sp, self.gspec["friction"]) * self.geom.swept_volume
        ca10, ca50, ca90 = burn_points(cl)
        return {"cl": cl, "W_closed": W, "p_max": p_max, "objective": obj, "ca10": ca10, "ca50": ca50,
                "ca90": ca90, "knock": knock_analysis(cl, rpm, self.kspec), "adv": adv}

    def optimise(self, st: Trapped, rpm: float) -> dict:
        """MBT, KLSA and the chosen spark for a trapped state at rpm."""
        s = self.sspec
        lo, hi = float(s["advance_min_deg"]), float(s["advance_max_deg"])
        hp = self.hp_closed(st, rpm)
        cache = {}

        def ev(a):
            a = float(a)
            if a not in cache:
                cache[a] = self.evaluate(st, rpm, a, hp)
            return cache[a]

        # MBT: coarse scan, then golden section around the best coarse point
        n_coarse = max(2, int(round((hi - lo) / float(s["mbt_coarse_step_deg"]))))
        grid = np.linspace(lo, hi, n_coarse + 1)
        vals = [ev(a)["objective"] for a in grid]
        k = int(np.argmax(vals))
        a, b = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
        c = b - GOLDEN * (b - a)
        d = a + GOLDEN * (b - a)
        while b - a > float(s["mbt_tol_deg"]):
            if ev(c)["objective"] > ev(d)["objective"]:
                b, d = d, c
                c = b - GOLDEN * (b - a)
            else:
                a, c = c, d
                d = a + GOLDEN * (b - a)
        mbt = 0.5 * (a + b)
        thr = float(self.kspec["integral_threshold"])
        klsa, status = klsa_search(lambda x: ev(x)["knock"]["margin"], lo, hi, float(s["klsa_tol_deg"]), thr,
                                   float(s["klsa_scan_step_deg"]))
        margin = float(self.kspec["klsa_margin_deg"])
        if status == "unlimited" or mbt <= klsa - margin:
            adv, limited = mbt, False
        else:
            adv, limited = max(klsa - margin, lo), True
        r = ev(adv)
        return {"mbt_deg": mbt, "klsa_deg": klsa, "klsa_status": status, "spark_deg": adv,
                "knock_limited": limited, "eval": r, "n_evals": len(cache), "ca50_mbt": ev(mbt)["ca50"]}
