"""Friction mean effective pressure: Chen-Flynn correlation.

    FMEP = A + B p_max + C S_p + D S_p^2

S. K. Chen and P. F. Flynn, "Development of a Single Cylinder Compression Ignition Research Engine",
SAE Technical Paper 650733 (1965), give the linear form FMEP [bar] = 0.137 + 0.005 p_max [bar]
+ 0.162 S_p [m/s]; the quadratic mean-piston-speed term is the common extension used by 1D engine codes.
Coefficients live in spec/gasdyn.yaml `friction` (SI units). Brake = net indicated (720 deg) - FMEP.

`fmep_total` adds, when `friction.piston_cooling_jets` is true, the oil-pump power of the piston cooling
jets (M6, ADR 0006; thermo/oiljet.py: jet flow x supply pressure / pump efficiency, zero below the jet
check-valve speed). The Chen-Flynn constant already stands for the base oil pump; the jets are extra delivery.
"""
from __future__ import annotations


def fmep_chen_flynn(p_max: float, mean_piston_speed: float, fr: dict) -> float:
    """FMEP [Pa] for peak cylinder pressure [Pa] and mean piston speed [m/s]."""
    if fr["model"] != "chen-flynn":
        raise ValueError(fr["model"])
    return (fr["const"] + fr["peak_pressure_factor"] * p_max + fr["speed_factor"] * mean_piston_speed
            + fr["speed_sq_factor"] * mean_piston_speed**2)


def fmep_total(p_max: float, mean_piston_speed: float, fr: dict) -> float:
    """Chen-Flynn FMEP plus the piston-cooling-jet pump FMEP [Pa] (if enabled in spec/gasdyn.yaml)."""
    return fmep_chen_flynn(p_max, mean_piston_speed, fr) + fmep_jets(mean_piston_speed, fr)


def fmep_jets(mean_piston_speed: float, fr: dict) -> float:
    """Piston-cooling-jet pump FMEP [Pa] (thermo/oiljet.py), zero when friction.piston_cooling_jets is off."""
    if not fr.get("piston_cooling_jets", False):
        return 0.0
    from thermo.oiljet import fmep_jets_from_piston_speed
    return fmep_jets_from_piston_speed(mean_piston_speed)
