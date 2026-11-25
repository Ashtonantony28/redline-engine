"""Plain-bearing half shells for the main journals and the big ends.

Each shell is a half tube: outer radius = housing bore (block saddle / main cap, rod / rod cap),
inner radius = journal radius + the running radial clearance of spec/dynamics.yaml `bearings`
(main 25 um, big end 22 um), width = the bearing land length of the same spec. The upper/lower
halves split on the housing split plane (z = 0 in the housing frame). Main shells are built at
the main-journal stations in engine frame; big-end shells in conrod-local frame.

The journal diameters (52 / 45 mm) come from spec/cad.yaml `crankshaft`; spec/dynamics.yaml uses
the same diameters for its oil-film model (checked by tests/test_m5_dynamics.py).
"""
from __future__ import annotations

from cad.conrod import big_end_bore_r
from cad.maincap import main_bore_r
from cad.params import Params
from cad.primitives import box, tube_x

BIG = 1.0e3


def _half(r_out: float, r_in: float, x: float, w: float, upper: bool):
    t = tube_x(r_out, r_in, x - w / 2, x + w / 2)
    return t & (box(-BIG, BIG, -BIG, BIG, 0.0, BIG) if upper else box(-BIG, BIG, -BIG, BIG, -BIG, 0.0))


def main_shell(p: Params, j: int, upper: bool):
    b = p.dyn["bearings"]["main"]
    r_in = p.L(p.section("crankshaft")["main_diameter"]) / 2 + p.L(b["radial_clearance"])
    return _half(main_bore_r(p), r_in, p.main_x(j), p.L(b["length"]), upper)


def big_end_shell(p: Params, upper: bool):
    b = p.dyn["bearings"]["big_end"]
    r_in = p.L(p.section("crankshaft")["pin_diameter"]) / 2 + p.L(b["radial_clearance"])
    return _half(big_end_bore_r(p), r_in, 0.0, p.L(b["length"]), upper)
