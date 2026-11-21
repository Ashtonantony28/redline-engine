"""spec/vehicle.yaml access plus derived quantities (test mass, rotating inertia). MODELLED values."""
from __future__ import annotations

import math
from functools import lru_cache

from redline.spec import ROOT, load

VEHICLE_SPEC = ROOT / "spec" / "vehicle.yaml"
SEC_PER_MIN = 60.0
KMH_PER_MS = 3.6
RAD_S_PER_RPM = 2.0 * math.pi / SEC_PER_MIN


@lru_cache(maxsize=None)
def vehicle_spec(path: str | None = None) -> dict:
    return load(str(path or VEHICLE_SPEC))


def engine_spec() -> dict:
    return load()


def wltp_test_mass(vs: dict | None = None) -> float:
    """GTR 15 test mass: running order + 25 kg + 15 % of max payload (optional equipment = 0)."""
    v = (vs or vehicle_spec())["vehicle"]
    return v["mass_running_order"] + v["wltp_extra_mass"] + v["wltp_payload_share"] * v["max_payload"]


def engine_inertia(es: dict | None = None, vs: dict | None = None) -> float:
    """Engine-side rotating inertia [kg m^2]: crank + flywheel (spec/engine.yaml) + clutch/input shaft."""
    m = (es or engine_spec())["masses"]
    return m["crank_inertia"] + m["flywheel_inertia"] + (vs or vehicle_spec())["driveline"]["gearbox_inertia_input"]


def overall_ratio(gear: int, vs: dict | None = None) -> float:
    """Engine speed / wheel speed in 1-based gear."""
    d = (vs or vehicle_spec())["driveline"]
    return d["gear_ratios"][gear - 1] * d["final_drive"]


def equivalent_mass(mass: float, gear: int | None, vs: dict | None = None, es: dict | None = None) -> float:
    """Translational mass + rotating inertias referred to the wheel rim [kg].

    gear None = clutch open / neutral (engine decoupled), else the engine inertia is reflected
    through the overall ratio squared (losses ignored for inertia, as in [G&S] eq. 2.6)."""
    vs = vs or vehicle_spec()
    v = vs["vehicle"]
    r = v["wheel_radius"]
    m = mass + v["wheels"] * v["wheel_inertia_each"] / r**2
    if gear is not None:
        m += engine_inertia(es, vs) * overall_ratio(gear, vs) ** 2 / r**2
    return m


def engine_rpm(speed: float, gear: int, vs: dict | None = None) -> float:
    """Engine speed [rpm] locked to vehicle speed [m/s] in a gear (no clutch slip)."""
    vs = vs or vehicle_spec()
    return speed / vs["vehicle"]["wheel_radius"] * overall_ratio(gear, vs) / RAD_S_PER_RPM
