"""Parameter access for the MVEM: spec/engine.yaml (single source of truth) + spec/mvem.yaml."""
from __future__ import annotations

import math

from redline import spec as rspec

MVEM_SPEC_PATH = rspec.ROOT / "spec" / "mvem.yaml"
ECU_SPEC_PATH = rspec.ROOT / "spec" / "ecu.yaml"

CYCLE_DEG = 720.0
REV_DEG = 360.0
SEC_PER_MIN = 60.0
REV_PER_CYCLE = 2.0
RPM_TO_RAD_S = 2.0 * math.pi / SEC_PER_MIN


def engine_spec() -> dict:
    return rspec.load()


def mvem_spec(path: str | None = None) -> dict:
    return rspec.load(str(path or MVEM_SPEC_PATH))


def ecu_spec(path: str | None = None) -> dict:
    return rspec.load(str(path or ECU_SPEC_PATH))


def geometry():
    return rspec.geometry(engine_spec())


def ignition_delay_deg(spec: dict | None = None) -> float:
    """Spark-to-SOC delay [deg] implied by the M2 reference point in spec/engine.yaml."""
    c = (spec or engine_spec())["combustion"]
    return float(c["spark_advance_deg"]) + float(c["heat_release_start_deg"])


def firing_tdc_deg(spec: dict | None = None) -> list[float]:
    """Firing-TDC cycle angle of cylinders 1..N (ADR 0001: cyl k fires 720/N * position in order)."""
    e = (spec or engine_spec())["engine"]
    n = int(e["cylinders"])
    out = [0.0] * n
    for pos, cyl in enumerate(e["firing_order"]):
        out[cyl - 1] = CYCLE_DEG / n * pos
    return out
