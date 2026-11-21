"""WLTC Class 3b speed trace (UN GTR 15), committed as vehicle/data/wltc_class3b.csv (source in its header)."""
from __future__ import annotations

import numpy as np

from redline.spec import ROOT
from vehicle.params import KMH_PER_MS, vehicle_spec


def load_trace(vs: dict | None = None):
    """(t [s], v [m/s]) at 1 Hz, t = 0..1800 s."""
    path = ROOT / (vs or vehicle_spec())["simulation"]["wltc_table"]
    data = np.loadtxt(path, delimiter=",", comments="#", skiprows=0, dtype=str)
    rows = np.array([r for r in data if not r[0].startswith("t_")], dtype=float)
    return rows[:, 0], rows[:, 1] / KMH_PER_MS


def distance(t, v) -> float:
    """Trapezoidal distance [m] (the GTR 15 convention for a 1 Hz trace)."""
    return float(np.trapezoid(v, t))


def phase_distances(t, v, vs: dict | None = None) -> dict:
    """{phase: distance [m]} over the GTR 15 phase spans."""
    out = {}
    for name, ph in (vs or vehicle_spec())["simulation"]["wltc_phases"].items():
        a, b = ph["span_s"]
        m = (t >= a) & (t <= b)
        out[name] = distance(t[m], v[m])
    return out
