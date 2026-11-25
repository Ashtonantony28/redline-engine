"""CAD parameters: spec/engine.yaml (single source of truth) + spec/cad.yaml (CAD detail), with the valvetrain,
dynamics, gas-dynamics and FEA specs read where the CAD must agree with them (valve angles and piston-to-valve
gap from spec/valvetrain.yaml, bearing clearances from spec/dynamics.yaml, runner/primary sizes from
spec/gasdyn.yaml, main-cap material from spec/fea.yaml).

`P` exposes every length already converted to CAD units (mm). Angles stay in degrees at this
boundary (build123d rotations take degrees) and are suffixed `_deg`.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from redline.spec import ROOT, geometry, load

CAD_SPEC_PATH = ROOT / "spec" / "cad.yaml"
VT_SPEC_PATH = ROOT / "spec" / "valvetrain.yaml"
DYN_SPEC_PATH = ROOT / "spec" / "dynamics.yaml"
GAS_SPEC_PATH = ROOT / "spec" / "gasdyn.yaml"
FEA_SPEC_PATH = ROOT / "spec" / "fea.yaml"
AUTO = "auto"
DEG_PER_REV = 360.0
CHAMBER_GRID_MM = 0.1         # numeric integration step of the chamber-volume solve (CAD units)
CHAMBER_SOLVE_TOL_MM = 1e-3   # bisection tolerance on the chamber wall margin
FOOTPRINT_SAMPLES = 96        # points per ellipse when building the chamber footprint hulls


def _spec(path) -> dict:
    return load(str(path))


@dataclass(frozen=True, eq=False)  # hashed by identity so parts can be cached per Params
class Params:
    engine: dict
    cad: dict
    k: float  # CAD length units per metre

    # ---- helpers -------------------------------------------------------------------------
    def L(self, metres: float) -> float:
        """Spec length (m) -> CAD length (mm)."""
        return float(metres) * self.k

    def to_si_length(self, cad_len: float) -> float:
        return cad_len / self.k

    def section(self, name: str) -> dict:
        return self.cad[name]

    @property
    def vt(self) -> dict:
        return _spec(VT_SPEC_PATH)

    @property
    def dyn(self) -> dict:
        return _spec(DYN_SPEC_PATH)

    @property
    def gas(self) -> dict:
        return _spec(GAS_SPEC_PATH)

    @property
    def fea(self) -> dict:
        return _spec(FEA_SPEC_PATH)

    def material(self, part: str) -> dict:
        """Material record of a part: engine.yaml materials first, then CAD-only material_data in cad.yaml
        (an entry `{from: fea.materials.main_cap}` points into spec/fea.yaml instead of repeating it)."""
        key = self.cad["materials"][part]
        if key in self.engine["materials"]:
            return self.engine["materials"][key]
        rec = self.cad["material_data"][key]
        if "from" in rec:
            node = {"fea": self.fea, "engine": self.engine}[rec["from"].split(".")[0]]
            for k in rec["from"].split(".")[1:]:
                node = node[k]
            return node
        return rec

    def density(self, part: str) -> float:
        """kg/m^3 for a part, via the material mapping in cad.yaml."""
        return float(self.material(part)["rho"])

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

    @property
    def tdc_pin_z(self) -> float:
        return self.crank_r + self.rod_len

    # ---- combustion chamber: asymmetric pent roof from the valve angles ---------------------
    def valve_angle_deg(self, kind: str) -> float:
        """Valve axis angle from the bore axis (spec/valvetrain.yaml); each roof face is normal to its valves."""
        return float(self.vt["combustion_chamber"][f"{kind}_valve_angle_deg"])

    def side_sign(self, kind: str) -> float:
        return -1.0 if kind == "intake" else 1.0

    def roof_z(self, y: float) -> float:
        """Height of the roof plane over the head face at y: each face meets the head face at the bore edge."""
        t_in = math.tan(math.radians(self.valve_angle_deg("intake")))
        t_ex = math.tan(math.radians(self.valve_angle_deg("exhaust")))
        R = self.bore_r
        return self.head_face_z + min((R + y) * t_in, (R - y) * t_ex)

    @property
    def ridge_y(self) -> float:
        """y of the roof ridge (intersection of the two roof faces)."""
        t_in = math.tan(math.radians(self.valve_angle_deg("intake")))
        t_ex = math.tan(math.radians(self.valve_angle_deg("exhaust")))
        return self.bore_r * (t_ex - t_in) / (t_ex + t_in)

    @property
    def roof_height(self) -> float:
        return self.roof_z(self.ridge_y) - self.head_face_z

    def valve_y(self, kind: str) -> float:
        """y of the valve-face centre, solved so the closed valve face's outer edge sits valve_to_crown_at_tdc
        from the flat crown along the valve axis (the gap spec/valvetrain.yaml and the M3 clearance model use)."""
        a = math.radians(self.valve_angle_deg(kind))
        d = self.L(self.engine["valves"][f"{kind}_head_diameter"])
        g0 = self.L(self.vt["combustion_chamber"]["valve_to_crown_at_tdc"])
        gap = self.head_face_z - self.crown_tdc_z
        y_edge = self.bore_r - (g0 * math.cos(a) - gap) / math.tan(a)
        y = y_edge - d / 2 * math.cos(a)
        inner = y - d / 2 * math.cos(a)
        margin = self.L(self.cad["head"]["valve_ridge_margin"])
        ridge = self.side_sign(kind) * self.ridge_y
        if inner < max(ridge, 0.0) + margin - 1e-9:
            raise ValueError(f"{kind} valve crosses the roof ridge (inner edge {inner:.2f} mm)")
        return self.side_sign(kind) * y

    def valve_layout(self):
        """[(kind, side_sign, x_local, y_local, head_diameter)] per cylinder, in the head plane.

        Intake on -y, exhaust on +y; each valve's axis is normal to its roof face. x/y are the
        horizontal positions where the valve axis pierces its roof face.
        """
        h = self.cad["head"]
        v = self.engine["valves"]
        out = []
        for kind in ("intake", "exhaust"):
            d = self.L(v[f"{kind}_head_diameter"])
            n = v["per_cylinder"][kind]
            pitch = d + self.L(h["valve_gap"][kind])
            y = self.valve_y(kind)
            for i in range(n):
                x = (i - (n - 1) / 2) * pitch
                out.append((kind, self.side_sign(kind), x, y, d))
        return out

    def valve_axis(self, kind: str) -> np.ndarray:
        """Unit vector along a valve axis, from the face up the stem."""
        a = math.radians(self.valve_angle_deg(kind))
        return np.array([0.0, self.side_sign(kind) * math.sin(a), math.cos(a)])

    def valve_origin(self, idx: int, kind: str, x_local: float) -> np.ndarray:
        y = self.valve_y(kind)
        return np.array([self.cyl_x(idx) + x_local, y, self.roof_z(y)])

    # ---- chamber footprint (squish lands outside it) ------------------------------------------
    def max_lift(self, kind: str) -> float:
        return self.L(self.engine["valves"][f"{kind}_max_lift"])

    def footprint_hulls(self, margin: float):
        """Convex outlines (x, y) around each valve of one cylinder (local x): the valve-head ellipse in plan swept
        by the opening stroke (the head moves along its axis, so in plan toward the ridge) and offset by `margin`."""
        hulls = []
        for kind, sign, x, y, d in self.valve_layout():
            a = math.radians(self.valve_angle_deg(kind))
            ax, ay = d / 2 + margin, d / 2 * math.cos(a) + margin
            shift = -sign * self.max_lift(kind) * math.sin(a)       # toward the bore axis
            t = np.linspace(0.0, 2 * math.pi, FOOTPRINT_SAMPLES, endpoint=False)
            pts = np.vstack([np.c_[x + ax * np.cos(t), y + ay * np.sin(t)],
                             np.c_[x + ax * np.cos(t), y + shift + ay * np.sin(t)]])
            hulls.append(_convex_hull(pts))
        return hulls

    def core_half(self):
        """Half extents (x, y) of the central chamber strip that joins the valve recesses around the plug."""
        h = self.cad["head"]
        return self.L(h["chamber_core_half_length"]), self.L(h["chamber_core_half_width"])

    def chamber_volume_numeric(self, margin: float) -> float:
        """Roof volume over the head face inside the footprint (mm^3), by grid integration."""
        R, step = self.bore_r, CHAMBER_GRID_MM
        g = np.arange(-R, R, step) + step / 2
        X, Y = np.meshgrid(g, g)
        inside = X**2 + Y**2 < R**2
        foot = np.zeros_like(inside)
        for hull in self.footprint_hulls(margin):
            foot |= _in_convex(hull, X, Y)
        cx, cy = self.core_half()
        foot |= (np.abs(X) < cx) & (np.abs(Y - self.ridge_y) < cy)
        t_in = math.tan(math.radians(self.valve_angle_deg("intake")))
        t_ex = math.tan(math.radians(self.valve_angle_deg("exhaust")))
        H = np.minimum((R + Y) * t_in, (R - Y) * t_ex)
        return float((H * (inside & foot)).sum() * step * step)

    @property
    def clearance_volume_target(self) -> float:
        """Geometric clearance volume from the spec CR (CAD units^3)."""
        return self.geom.clearance_volume * self.k**3

    def pocket_volume(self) -> float:
        """Volume of the crown valve pockets of one piston (mm^3), from the pocket solids."""
        from cad.piston import valve_pockets
        pk = valve_pockets(self)
        return pk.volume if pk is not None else 0.0

    @property
    def chamber_margin(self) -> float:
        """Chamber wall offset around the valve sweep, solved so the CAD compression ratio equals the spec one."""
        h = self.cad["head"]
        if h["chamber_margin"] != AUTO:
            return self.L(h["chamber_margin"])
        return _solve_margin(self)


@lru_cache(maxsize=None)
def _solve_margin(p: Params) -> float:
    cyl_gap = math.pi * p.bore_r**2 * (p.head_face_z - p.crown_tdc_z)
    target = p.clearance_volume_target - cyl_gap - p.pocket_volume()
    lo, hi = 0.0, p.bore_r
    if p.chamber_volume_numeric(lo) > target:
        raise ValueError("chamber with zero margin already exceeds the clearance volume: CR not reachable")
    while hi - lo > CHAMBER_SOLVE_TOL_MM:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if p.chamber_volume_numeric(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def _convex_hull(pts: np.ndarray) -> np.ndarray:
    """Andrew's monotone chain; returns hull vertices counter-clockwise."""
    pts = np.unique(np.round(pts, 9), axis=0)
    pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower, upper = [], []
    for q in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in pts[::-1]:
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    return np.array(lower[:-1] + upper[:-1])


def _in_convex(hull: np.ndarray, X, Y):
    inside = np.ones_like(X, dtype=bool)
    n = len(hull)
    for i in range(n):
        x0, y0 = hull[i]
        x1, y1 = hull[(i + 1) % n]
        inside &= (x1 - x0) * (Y - y0) - (y1 - y0) * (X - x0) >= 0
    return inside


@lru_cache(maxsize=None)
def params(engine_path: str | None = None, cad_path: str | None = None) -> Params:
    cad = load(str(cad_path or CAD_SPEC_PATH))
    return Params(load(engine_path), cad, cad["units"]["cad_length_per_m"])
