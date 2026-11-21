"""Parameter access for the M5 dynamics package.

Geometry, masses, firing order and combustion come from spec/engine.yaml (single source of
truth). Dynamics-only parameters (bearings, torsional stations, damper, flywheel targets, the
placeholder pressure generator) come from spec/dynamics.yaml. ADR 0001 conventions throughout.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from redline import spec as rspec

DYN_SPEC_PATH = rspec.ROOT / "spec" / "dynamics.yaml"

FULL_CYCLE = 4 * math.pi      # one four-stroke cycle [rad]
REV = 2 * math.pi             # one crank revolution [rad]
RPM_TO_RAD_S = 2 * math.pi / 60.0


def engine_spec(path: str | None = None) -> dict:
    return rspec.load(path)


def dyn_spec(path: str | None = None) -> dict:
    return rspec.load(str(path or DYN_SPEC_PATH))


def geometry(spec: dict | None = None) -> rspec.Geometry:
    return rspec.geometry(spec or engine_spec())


def omega_from_rpm(rpm: float) -> float:
    return rpm * RPM_TO_RAD_S


@dataclass(frozen=True)
class Masses:
    """Two-mass conrod lumping from the spec `masses` budget."""

    piston_assembly: float
    conrod: float
    small_end_fraction: float

    @property
    def reciprocating(self) -> float:
        """m_rec: piston assembly + small-end share of the conrod [kg]."""
        return self.piston_assembly + self.small_end_fraction * self.conrod

    @property
    def rod_big_end(self) -> float:
        """Rotating share of the conrod, lumped at the crankpin [kg]."""
        return (1.0 - self.small_end_fraction) * self.conrod


def masses(spec: dict | None = None) -> Masses:
    m = (spec or engine_spec())["masses"]
    return Masses(m["piston_assembly"], m["conrod"], m["conrod_small_end_fraction"])


def firing_offsets(spec: dict | None = None) -> np.ndarray:
    """Firing phase offset [rad] of cylinders 1..n (index 0..n-1), ADR 0001.

    Cylinder k fires at 180 deg * (its position in the firing order) for an I4,
    generally 720/n deg per position.
    """
    e = (spec or engine_spec())["engine"]
    order = list(e["firing_order"])
    n = e["cylinders"]
    step = FULL_CYCLE / n
    return np.array([step * order.index(k + 1) for k in range(n)])


def throw_angles(spec: dict | None = None) -> np.ndarray:
    """Geometric crank-throw angle of each cylinder [rad] (firing offset modulo one revolution)."""
    return np.mod(firing_offsets(spec), REV)


def cylinder_positions(spec: dict | None = None, dspec: dict | None = None) -> np.ndarray:
    """Axial (x, front to rear) position of each bore centre relative to the engine centre [m]."""
    n = (spec or engine_spec())["engine"]["cylinders"]
    pitch = (dspec or dyn_spec())["layout"]["cylinder_pitch"]
    return (np.arange(n) - (n - 1) / 2.0) * pitch


def shear_modulus(spec: dict | None = None, dspec: dict | None = None) -> float:
    e_mod = float((spec or engine_spec())["materials"]["steel_crank"]["E"])  # PyYAML reads "210.0e9" as str
    nu = (dspec or dyn_spec())["crank"]["poisson"]
    return e_mod / (2.0 * (1.0 + nu))
