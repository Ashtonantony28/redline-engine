"""Mass properties of CAD solids in SI units (ADR 0001).

build123d returns volume-based properties in CAD units (mm); this module multiplies by the
spec density and converts to kg, m and kg m^2. Inertia about a parallel axis uses the
parallel-axis theorem from the CoM tensor.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class MassProps:
    mass: float                       # kg
    com: np.ndarray                   # m
    inertia_com: np.ndarray           # kg m^2, about the CoM, body axes = global axes
    parts: list = field(default_factory=list)

    def inertia_about(self, point) -> np.ndarray:
        """Inertia tensor about `point` (m) with global axes (parallel-axis theorem)."""
        d = self.com - np.asarray(point, dtype=float)
        return self.inertia_com + self.mass * (np.dot(d, d) * np.eye(3) - np.outer(d, d))

    def first_moment(self) -> np.ndarray:
        """m * r_CoM (kg m)."""
        return self.mass * self.com

    def __add__(self, other: "MassProps") -> "MassProps":
        m = self.mass + other.mass
        com = (self.mass * self.com + other.mass * other.com) / m
        inertia = self.inertia_about(com) + other.inertia_about(com)
        return MassProps(m, com, inertia, self.parts + other.parts)

    def to_dict(self, axis_point=None) -> dict:
        out = {
            "mass_kg": self.mass,
            "com_m": self.com.tolist(),
            "inertia_about_com_kgm2": self.inertia_com.tolist(),
        }
        if axis_point is not None:
            out["inertia_about_crank_axis_point_kgm2"] = self.inertia_about(axis_point).tolist()
        return out


def mass_props(shape, density: float, cad_length_per_m: float, name: str = "") -> MassProps:
    k = cad_length_per_m
    vol = shape.volume / k**3
    mass = density * vol
    com = np.array(tuple(shape.center()), dtype=float) / k
    inertia = np.array(shape.matrix_of_inertia, dtype=float) * density / k**5
    return MassProps(mass, com, inertia, [name] if name else [])
