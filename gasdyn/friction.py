"""Friction mean effective pressure: Chen-Flynn correlation.

    FMEP = A + B p_max + C S_p + D S_p^2

S. K. Chen and P. F. Flynn, "Development of a Single Cylinder Compression Ignition Research Engine",
SAE Technical Paper 650733 (1965), give the linear form FMEP [bar] = 0.137 + 0.005 p_max [bar]
+ 0.162 S_p [m/s]; the quadratic mean-piston-speed term is the common extension used by 1D engine codes.
Coefficients live in spec/gasdyn.yaml `friction` (SI units). Brake = net indicated (720 deg) - FMEP.
"""
from __future__ import annotations


def fmep_chen_flynn(p_max: float, mean_piston_speed: float, fr: dict) -> float:
    """FMEP [Pa] for peak cylinder pressure [Pa] and mean piston speed [m/s]."""
    if fr["model"] != "chen-flynn":
        raise ValueError(fr["model"])
    return (fr["const"] + fr["peak_pressure_factor"] * p_max + fr["speed_factor"] * mean_piston_speed
            + fr["speed_sq_factor"] * mean_piston_speed**2)
