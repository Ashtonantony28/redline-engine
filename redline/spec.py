"""Load spec/engine.yaml, the single source of truth (ADR 0001)."""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import math

import yaml

ROOT = Path(__file__).resolve().parent.parent
SPEC_PATH = ROOT / "spec" / "engine.yaml"


@lru_cache(maxsize=None)
def load(path: str | None = None) -> dict:
    with open(path or SPEC_PATH) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class Geometry:
    bore: float
    stroke: float
    conrod_length: float
    compression_ratio: float
    cylinders: int

    @property
    def crank_radius(self) -> float:
        return self.stroke / 2

    @property
    def piston_area(self) -> float:
        return math.pi * self.bore**2 / 4

    @property
    def swept_volume(self) -> float:
        return self.piston_area * self.stroke

    @property
    def clearance_volume(self) -> float:
        return self.swept_volume / (self.compression_ratio - 1)

    @property
    def displacement(self) -> float:
        return self.swept_volume * self.cylinders

    @property
    def lam(self) -> float:
        """Crank radius / conrod length."""
        return self.crank_radius / self.conrod_length


def geometry(spec: dict | None = None) -> Geometry:
    e = (spec or load())["engine"]
    return Geometry(e["bore"], e["stroke"], e["conrod_length"], e["compression_ratio"], e["cylinders"])
