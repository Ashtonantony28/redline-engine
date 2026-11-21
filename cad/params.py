"""CAD parameters: spec/engine.yaml (single source of truth) + spec/cad.yaml (CAD detail).

`P` exposes every length already converted to CAD units (mm). Angles stay in degrees at this
boundary (build123d rotations take degrees) and are suffixed `_deg`.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

from redline.spec import ROOT, geometry, load

CAD_SPEC_PATH = ROOT / "spec" / "cad.yaml"
AUTO = "auto"
DEG_PER_REV = 360.0


@dataclass(frozen=True, eq=False)  # hashed by identity so parts can be cached per Params
class Params:
    engine: dict
    cad: dict
    k: float  # CAD length units per metre

    # ---- helpers -------------------------------------------------------------------------
    def L(self, metres: float) -> float:
        """Spec length (m) -> CAD length (mm)."""
        return metres * self.k

    def to_si_length(self, cad_len: float) -> float:
        return cad_len / self.k

    def section(self, name: str) -> dict:
        return self.cad[name]

    def density(self, part: str) -> float:
        """kg/m^3 for a part, via the material mapping in cad.yaml."""
        return self.engine["materials"][self.cad["materials"][part]]["rho"]

    # ---- engine geometry -----------------------------------------------------------------
    @property
    def geom(self):
        return geometry(self.engine)

    @property
    def bore(self) -> float:
        return self.L(self.geom.bore)

    @property
    def bore_r(self) -> float:
        return self.bore / 2

    @property
    def crank_r(self) -> float:
        return self.L(self.geom.crank_radius)

    @property
    def rod_len(self) -> float:
        return self.L(self.geom.conrod_length)

    @property
    def cylinders(self) -> int:
        return self.geom.cylinders

    @property
    def spacing(self) -> float:
        return self.L(self.cad["layout"]["bore_spacing"])

    def cyl_x(self, idx: int) -> float:
        """x of cylinder idx (0-based), engine centred on x = 0, cylinder 1 at the front (-x)."""
        return (idx - (self.cylinders - 1) / 2) * self.spacing

    def main_x(self, j: int) -> float:
        """x of main journal j (0-based); one main between and outside each cylinder."""
        return (j - self.cylinders / 2) * self.spacing

    @property
    def n_mains(self) -> int:
        return self.cylinders + 1

    def firing_phase_deg(self, idx: int) -> float:
        """Firing phase of cylinder idx (0-based) per ADR 0001: 720/n * position in firing order."""
        order = self.engine["engine"]["firing_order"]
        cycle_deg = 2 * DEG_PER_REV
        return order.index(idx + 1) * cycle_deg / self.cylinders

    def pin_angle_deg(self, idx: int, crank_deg: float = 0.0) -> float:
        """Mechanical crankpin angle from top (bore axis) of cylinder idx at a crank angle."""
        return (crank_deg - self.firing_phase_deg(idx)) % DEG_PER_REV

    # ---- stack-up along the bore axis (z up, crank axis z = 0) ---------------------------
    @property
    def compression_height(self) -> float:
        return self.L(self.cad["piston"]["compression_height"])

    @property
    def deck_z(self) -> float:
        h = self.cad["head"]
        return self.crank_r + self.rod_len + self.compression_height + self.L(h["piston_deck_clearance"])

    @property
    def head_face_z(self) -> float:
        return self.deck_z + self.L(self.cad["head"]["gasket_thickness"])

    @property
    def crown_tdc_z(self) -> float:
        return self.crank_r + self.rod_len + self.compression_height

    # ---- combustion chamber --------------------------------------------------------------
    @property
    def pent_roof_shape_factor(self) -> float:
        """V_roof = tan(theta) * R^3 * (pi - 4/3) for a ridge through the bore axis."""
        return math.pi - 4.0 / 3.0

    def valve_relief_volume(self) -> float:
        """Analytic volume of the optional crown valve reliefs (CAD units^3)."""
        p = self.cad["piston"]
        if not p["valve_reliefs"]:
            return 0.0
        v = self.engine["valves"]
        depth = self.L(p["valve_relief_depth"])
        total = 0.0
        for kind in ("intake", "exhaust"):
            d = self.L(v[f"{kind}_head_diameter"] + p["valve_relief_margin"])
            total += v["per_cylinder"][kind] * math.pi * d**2 / 4 * depth
        return total

    @property
    def clearance_volume_target(self) -> float:
        """Geometric clearance volume from the spec CR (CAD units^3)."""
        return self.geom.clearance_volume * self.k**3

    @property
    def roof_angle_deg(self) -> float:
        h = self.cad["head"]
        if h["roof_angle_deg"] != AUTO:
            return float(h["roof_angle_deg"])
        cyl_gap = self.head_face_z - self.crown_tdc_z
        v_cyl = math.pi * self.bore_r**2 * cyl_gap
        v_roof = self.clearance_volume_target - v_cyl - self.valve_relief_volume()
        tan_t = v_roof / (self.bore_r**3 * self.pent_roof_shape_factor)
        return math.degrees(math.atan(tan_t))

    @property
    def roof_height(self) -> float:
        return self.bore_r * math.tan(math.radians(self.roof_angle_deg))

    def valve_layout(self):
        """[(kind, side_sign, x_local, y_local, head_diameter)] per cylinder, in the head plane.

        Intake on -y, exhaust on +y; each valve's axis is normal to its roof face. x/y are the
        horizontal positions where the valve axis pierces its roof face.
        """
        h = self.cad["head"]
        v = self.engine["valves"]
        cos_t = math.cos(math.radians(self.roof_angle_deg))
        out = []
        for kind, sign in (("intake", -1.0), ("exhaust", 1.0)):
            d = self.L(v[f"{kind}_head_diameter"])
            n = v["per_cylinder"][kind]
            pitch = d + self.L(h["valve_gap"])
            y = sign * (d / 2 * cos_t + self.L(h["valve_ridge_margin"]))
            for i in range(n):
                x = (i - (n - 1) / 2) * pitch
                out.append((kind, sign, x, y, d))
        return out


@lru_cache(maxsize=None)
def params(engine_path: str | None = None, cad_path: str | None = None) -> Params:
    cad = load(str(cad_path or CAD_SPEC_PATH))
    return Params(load(engine_path), cad, cad["units"]["cad_length_per_m"])
