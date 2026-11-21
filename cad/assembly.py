"""Engine assembly at a given crank angle (ADR 0001: 0 deg = firing TDC of cylinder 1).

The crankshaft rotates about +x by -theta so that the crankpin of cylinder k sits at
(y, z) = r (sin(theta - phase_k), cos(theta - phase_k)); each conrod is rotated about its
big-end axis by the exact slider-crank conrod angle, and each piston sits at the exact
slider-crank pin height (kin/slider_crank.py, the M0-verified reference).
"""
from __future__ import annotations

from functools import lru_cache
import math

from build123d import Compound, Pos, Rot

from cad.block import block
from cad.conrod import conrod
from cad.crankshaft import crankshaft
from cad.head import head
from cad.params import Params, params
from cad.piston import piston_parts
from kin.slider_crank import piston_position


@lru_cache(maxsize=None)
def parts(p: Params | None = None) -> dict:
    """Part solids in their local frames: {name: (shape, material key)}. Cached."""
    p = p or params()
    out = {name: v for name, v in piston_parts(p).items()}
    out["conrod"] = (conrod(p), "conrod")
    out["crankshaft"] = (crankshaft(p), "crankshaft")
    out["block"] = (block(p), "block")
    out["head"] = (head(p), "head")
    return out


def cylinder_pose(p: Params, idx: int, crank_deg: float):
    """(pin_y, pin_z, rod_angle_deg, wrist_pin_z) for cylinder idx at a crank angle."""
    a = math.radians(p.pin_angle_deg(idx, crank_deg))
    r, l = p.crank_r, p.rod_len
    py, pz = r * math.sin(a), r * math.cos(a)
    wz = float(piston_position(a, r, l))
    rod_deg = math.degrees(math.asin(py / l))  # Rot about x by rod_deg maps +z to (-sin, cos)
    return py, pz, rod_deg, wz


def placed_parts(crank_deg: float = 0.0, p: Params | None = None) -> dict:
    """{instance name: (located shape, material key)} at the given crank angle."""
    p = p or params()
    src = parts(p)
    out = {
        "block": src["block"],
        "head": src["head"],
        "crankshaft": (Rot(-crank_deg, 0, 0) * src["crankshaft"][0], "crankshaft"),
    }
    piston_names = [n for n in src if n not in ("conrod", "crankshaft", "block", "head")]
    for i in range(p.cylinders):
        py, pz, rod_deg, wz = cylinder_pose(p, i, crank_deg)
        x = p.cyl_x(i)
        out[f"conrod_{i + 1}"] = (Pos(x, py, pz) * Rot(rod_deg, 0, 0) * src["conrod"][0], "conrod")
        for n in piston_names:
            shape, mat = src[n]
            out[f"{n}_{i + 1}"] = (Pos(x, 0.0, wz) * shape, mat)
    return out


def assembly(crank_deg: float = 0.0, p: Params | None = None) -> Compound:
    children = []
    for name, (shape, _mat) in placed_parts(crank_deg, p).items():
        shape.label = name
        children.append(shape)
    return Compound(children=children, label=f"redline_engine_{crank_deg:g}deg")
