"""Engine maps for the MVEM, generated from thermo/ crank-angle cycle runs (M2 model).

`python -m mvem.maps` sweeps the coarse speed x MAP grid in spec/mvem.yaml and, at every point, the
Wiebe start of combustion (SOC). Each point is a converged two-zone cycle (thermo.cycle.run_cycle)
with the manifold pressure as the intake boundary and an exhaust backpressure that grows with flow.
Stored (committed fixture mvem/data/thermo_maps.json, reproducible bit-for-bit on one machine):
  eta_gross[rpm][map][soc]  gross indicated efficiency (compression + expansion strokes)
  ve[rpm][map]              volumetric efficiency referenced to manifold density (reference SOC)
  pmep[rpm][map]            pumping mean effective pressure [Pa] (reference SOC)
  soc_mbt[rpm][map]         SOC of maximum gross efficiency (parabolic refinement)
All values are MODELLED.
"""
from __future__ import annotations

import argparse
import copy
import json
from functools import lru_cache
from multiprocessing import Pool

import numpy as np

from mvem import params
from redline.spec import ROOT

FIXTURE_VERSION = "m8-maps-1"
LOAD_KEY = "mvem_grid"


def soc_grid(ms: dict | None = None) -> np.ndarray:
    s = (ms or params.mvem_spec())["maps"]["soc_deg"]
    n = int(round((float(s["stop"]) - float(s["start"])) / float(s["step"]))) + 1
    return float(s["start"]) + float(s["step"]) * np.arange(n)


def exhaust_pressure(rpm: float, p_map: float, ms: dict) -> float:
    ex, amb = ms["exhaust"], ms["ambient"]
    flow_ratio = (rpm / float(ex["rated_rpm"])) * (p_map / float(amb["p"]))
    return float(amb["p"]) + float(ex["dp_rated"]) * flow_ratio ** float(ex["dp_exponent"])


@lru_cache(maxsize=1)
def _tables():
    from thermo.props import build_tables, thermo_spec
    return build_tables(thermo_spec())


def run_point(args) -> list[dict]:
    """All SOC values at one (rpm, MAP) grid point; returns thermo metrics per SOC."""
    from thermo.cycle import run_cycle
    from thermo.props import thermo_spec
    rpm, p_map, socs = args
    ms = params.mvem_spec()
    spec = copy.deepcopy(params.engine_spec())
    ts = copy.deepcopy(thermo_spec())
    ts["loads"][LOAD_KEY] = {"p_intake": float(p_map), "T_intake": float(ms["manifold"]["temperature"]),
                             "p_exhaust": exhaust_pressure(rpm, p_map, ms),
                             "T_exhaust": float(ms["exhaust"]["temperature"])}
    out = []
    for soc in socs:
        spec["combustion"]["heat_release_start_deg"] = float(soc)
        r = run_cycle(float(rpm), LOAD_KEY, step_deg=float(ms["maps"]["step_deg"]), tables=_tables(),
                      spec=spec, tspec=ts)
        m = r.metrics
        out.append({"eta_gross": float(m["eta_ind_gross"]), "ve": float(m["volumetric_efficiency"]),
                    "pmep": float(m["pmep_bar"]) * 1.0e5, "converged": bool(m["converged"])})
    return out


def mbt_refine(soc: np.ndarray, eta: np.ndarray) -> float:
    i = int(np.argmax(eta))
    if i == 0 or i == len(eta) - 1:
        return float(soc[i])
    y0, y1, y2 = eta[i - 1], eta[i], eta[i + 1]
    den = y0 - 2.0 * y1 + y2
    h = soc[1] - soc[0]
    return float(soc[i] + (0.5 * h * (y0 - y2) / den if den < 0.0 else 0.0))


def generate(workers: int = 4) -> dict:
    ms = params.mvem_spec()
    rpms = [float(x) for x in ms["maps"]["rpm"]]
    maps = [float(x) for x in ms["maps"]["map"]]
    socs = soc_grid(ms)
    ref = float(params.engine_spec()["combustion"]["heat_release_start_deg"])
    iref = int(np.argmin(np.abs(socs - ref)))
    jobs = [(r, p, socs.tolist()) for r in rpms for p in maps]
    with Pool(workers) as pool:
        res = pool.map(run_point, jobs)
    nr, nm, ns = len(rpms), len(maps), len(socs)
    eta = np.zeros((nr, nm, ns))
    ve = np.zeros((nr, nm))
    pmep = np.zeros((nr, nm))
    conv = True
    for k, pts in enumerate(res):
        i, j = divmod(k, nm)
        eta[i, j] = [p["eta_gross"] for p in pts]
        ve[i, j] = pts[iref]["ve"]
        pmep[i, j] = pts[iref]["pmep"]
        conv &= all(p["converged"] for p in pts)
    mbt = np.array([[mbt_refine(socs, eta[i, j]) for j in range(nm)] for i in range(nr)])
    return {"version": FIXTURE_VERSION, "note": "MODELLED: thermo/ two-zone cycle runs, see mvem/maps.py",
            "rpm": rpms, "map": maps, "soc_deg": socs.tolist(), "soc_ref_deg": float(socs[iref]),
            "ignition_delay_deg": params.ignition_delay_deg(), "all_converged": bool(conv),
            "eta_gross": eta.tolist(), "ve": ve.tolist(), "pmep": pmep.tolist(), "soc_mbt": mbt.tolist()}


def fixture_path():
    return ROOT / params.mvem_spec()["maps"]["fixture"]


@lru_cache(maxsize=1)
def load_maps() -> dict:
    with open(fixture_path()) as f:
        d = json.load(f)
    for k in ("rpm", "map", "soc_deg", "eta_gross", "ve", "pmep", "soc_mbt"):
        d[k] = np.asarray(d[k], dtype=float)
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    d = generate(a.workers)
    p = fixture_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(d, f, indent=0, sort_keys=True)
        f.write("\n")
    print(f"wrote {p} (all_converged={d['all_converged']}); MODELLED")


if __name__ == "__main__":
    main()


def mbt_advance_table(d: dict | None = None) -> np.ndarray:
    """MBT spark advance [deg BTDC] on the (rpm, MAP) grid = ignition delay - SOC_MBT.

    Points whose efficiency maximum sits on the edge of the SOC sweep (residual-dominated
    corner at the lowest MAP, where the cycle barely burns) take the value of the next MAP column.
    """
    d = d or load_maps()
    soc = d["soc_deg"]
    mbt = np.array(d["soc_mbt"], dtype=float)
    bad = (mbt <= soc[0]) | (mbt >= soc[-1])
    for i in range(mbt.shape[0]):
        for j in range(mbt.shape[1] - 2, -1, -1):
            if bad[i, j]:
                mbt[i, j] = mbt[i, j + 1]
    return float(d["ignition_delay_deg"]) - mbt
