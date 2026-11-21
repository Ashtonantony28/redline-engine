#!/usr/bin/env python3
"""Apply the placement to redline_ecu: outline, parts, In1 ground plane, hole keep-outs, net classes.

Runs under KiCad 10's Python (pcbnew).  Re-runnable: it removes every track, via, zone and
board drawing and re-places every footprint from the placement JSON (layout/autoplace.py).

    layout/run.sh place      (wraps: python3 layout/place.py <board.kicad_pcb> <placement.json>)

Adapted from drone-design ecad/layout/place.py (outline, zones, keep-outs, net class patterns).
"""
import json
import math
import os
import sys

import pcbnew

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from board import CORNER_R, EDGE_INSET, H, HOLE_KEEPOUT_R, ORIGIN, POWER_NETS, W, WIDE_NETS, corner_holes  # noqa: E402

MM = pcbnew.FromMM
LIBRARY_CHECKS = ("lib_footprint_issues", "lib_footprint_mismatch")
SILK_REFS = ("J", "TP", "SW", "U", "Q", "D", "H", "Y")


def pt(x, y):
    return pcbnew.VECTOR2I(MM(ORIGIN[0] + x), MM(ORIGIN[1] + y))


def footprints_by_path(footprints):
    out = {}
    for fp in footprints:
        path = next((f.GetText() for f in fp.GetFields() if f.GetName() == "Path"), fp.GetReference())
        out[path] = fp
    return out


def put(fp, x, y, rot, side):
    if fp.IsFlipped() != (side == "B"):
        fp.Flip(fp.GetPosition(), pcbnew.FLIP_DIRECTION_LEFT_RIGHT)
    fp.SetOrientationDegrees(0)
    fp.SetPosition(pt(x, y))
    fp.SetOrientationDegrees(rot)


def set_outline(board):
    hx, hy, r = W / 2, H / 2, CORNER_R
    k = r * (1 - math.sqrt(0.5))

    def line(a, b):
        s = pcbnew.PCB_SHAPE(board, pcbnew.SHAPE_T_SEGMENT)
        s.SetStart(pt(*a))
        s.SetEnd(pt(*b))
        return s

    def arc(start, mid, end):
        s = pcbnew.PCB_SHAPE(board, pcbnew.SHAPE_T_ARC)
        s.SetArcGeometry(pt(*start), pt(*mid), pt(*end))
        return s

    shapes = [
        line((-hx + r, -hy), (hx - r, -hy)), arc((hx - r, -hy), (hx - k, -hy + k), (hx, -hy + r)),
        line((hx, -hy + r), (hx, hy - r)), arc((hx, hy - r), (hx - k, hy - k), (hx - r, hy)),
        line((hx - r, hy), (-hx + r, hy)), arc((-hx + r, hy), (-hx + k, hy - k), (-hx, hy - r)),
        line((-hx, hy - r), (-hx, -hy + r)), arc((-hx, -hy + r), (-hx + k, -hy + k), (-hx + r, -hy)),
    ]
    for s in shapes:
        s.SetLayer(pcbnew.Edge_Cuts)
        s.SetWidth(MM(0.1))
        board.Add(s)


def add_text(board):
    for text, x, y, layer, size in (("REDLINE ECU r1", 0.0, -42.5, pcbnew.F_SilkS, 1.2),
                                    ("J2 SENSORS", -30.7, 29.0, pcbnew.F_SilkS, 1.0),
                                    ("J1 POWER / ACTUATORS", 30.7, 29.5, pcbnew.F_SilkS, 1.0)):
        t = pcbnew.PCB_TEXT(board)
        t.SetText(text)
        t.SetLayer(layer)
        t.SetPosition(pt(x, y))
        t.SetTextSize(pcbnew.VECTOR2I(MM(size), MM(size)))
        t.SetTextThickness(MM(0.15))
        board.Add(t)


def tidy_text(footprints):
    """Passives get their designator on the fab layer; ICs, connectors and test points keep a silk label."""
    for fp in footprints:
        ref = fp.Reference()
        keep = fp.GetReference().startswith(SILK_REFS) and not fp.GetReference().startswith("TP")
        ref.SetLayer(pcbnew.F_SilkS if keep else pcbnew.F_Fab)
        ref.SetTextSize(pcbnew.VECTOR2I(MM(0.8), MM(0.8)))
        ref.SetTextThickness(MM(0.12))
        ref.SetVisible(True)
        fp.Value().SetLayer(pcbnew.F_Fab)
        fp.Value().SetVisible(False)
        if not keep:
            ref.SetPosition(fp.GetPosition())


def add_zone(board, net, layer, outline, priority=0, name="", clearance=0.2, thermal_tht=True, min_width=0.2):
    z = pcbnew.ZONE(board)
    z.SetLayer(layer)
    z.SetNetCode(board.GetNetcodeFromNetname(net))
    z.SetZoneName(name or net)
    z.SetAssignedPriority(priority)
    z.SetLocalClearance(MM(clearance))
    z.SetMinThickness(MM(min_width))
    # SMD pads solid; through-hole pins (header, holes) on thermal spokes so that they can be soldered
    z.SetPadConnection(pcbnew.ZONE_CONNECTION_THT_THERMAL if thermal_tht else pcbnew.ZONE_CONNECTION_FULL)
    z.SetThermalReliefGap(MM(0.3))
    z.SetThermalReliefSpokeWidth(MM(0.5))
    z.SetIslandRemovalMode(pcbnew.ISLAND_REMOVAL_MODE_ALWAYS)
    poly = z.Outline()
    poly.NewOutline()
    for x, y in outline:
        poly.Append(MM(ORIGIN[0] + x), MM(ORIGIN[1] + y))
    board.Add(z)
    return z


def add_keepout(board, name, outline, layers):
    z = pcbnew.ZONE(board)
    z.SetIsRuleArea(True)
    z.SetZoneName(name)
    ls = pcbnew.LSET()
    for layer in layers:
        ls.AddLayer(layer)
    z.SetLayerSet(ls)
    z.SetDoNotAllowTracks(True)
    z.SetDoNotAllowVias(True)
    z.SetDoNotAllowZoneFills(True)
    z.SetDoNotAllowPads(False)
    z.SetDoNotAllowFootprints(False)
    poly = z.Outline()
    poly.NewOutline()
    for x, y in outline:
        poly.Append(MM(ORIGIN[0] + x), MM(ORIGIN[1] + y))
    board.Add(z)


def board_rect(inset):
    """The outline inset by `inset`, corners chamfered so the pour never comes near the arcs."""
    hx, hy, c = W / 2 - inset, H / 2 - inset, CORNER_R
    return [(-hx + c, -hy), (hx - c, -hy), (hx, -hy + c), (hx, hy - c), (hx - c, hy), (-hx + c, hy), (-hx, hy - c), (-hx, -hy + c)]


def assign_netclasses(pro_path):
    with open(pro_path) as f:
        pro = json.load(f)
    patterns = pro["net_settings"].setdefault("netclass_patterns", [])
    patterns[:] = [p for p in patterns if p.get("netclass") != "Power"]
    patterns.extend({"netclass": "Power", "pattern": n} for n in sorted(set(POWER_NETS + WIDE_NETS + ("GND",))))
    # The footprints come from the Zener package cache, not a KiCad library table, and silk.py
    # trims their silkscreen on purpose: the two library comparisons say nothing here.
    severities = pro["board"]["design_settings"].setdefault("rule_severities", {})
    for check in LIBRARY_CHECKS:
        severities[check] = "ignore"
    with open(pro_path, "w") as f:
        json.dump(pro, f, indent=2)
        f.write("\n")


def main(board_path, placement_path):
    assign_netclasses(os.path.splitext(board_path)[0] + ".kicad_pro")
    board = pcbnew.LoadBoard(board_path)
    # take every container before the first Remove(): the SWIG list accessors stop working after it
    footprints = list(board.GetFootprints())
    stale = list(board.GetTracks()) + list(board.Zones())
    stale += [d for d in board.GetDrawings() if d.GetLayer() in (pcbnew.Edge_Cuts, pcbnew.F_SilkS, pcbnew.B_SilkS)]
    for item in stale:
        board.Remove(item)
    set_outline(board)
    add_text(board)
    fps = footprints_by_path(footprints)
    placement = json.load(open(placement_path))
    missing = sorted(set(fps) - set(placement))
    assert not missing, "not placed: %s" % ", ".join(missing)
    for path, (x, y, rot, side) in placement.items():
        put(fps[path], x, y, rot, side)
    tidy_text(footprints)
    # L2: one solid, unbroken ground plane
    add_zone(board, "GND", pcbnew.In1_Cu, board_rect(EDGE_INSET), name="L2_GND", clearance=0.2)
    cu = [pcbnew.F_Cu, pcbnew.In1_Cu, pcbnew.In2_Cu, pcbnew.B_Cu]
    for i, (hx, hy) in enumerate(corner_holes()):
        circle = [(hx + HOLE_KEEPOUT_R * math.cos(a * math.pi / 18), hy + HOLE_KEEPOUT_R * math.sin(a * math.pi / 18)) for a in range(36)]
        add_keepout(board, "HOLE%d_KEEPOUT" % (i + 1), circle, cu)
    pcbnew.ZONE_FILLER(board).Fill(board.Zones())
    board.Save(board_path)
    print("placed %d footprints" % len(placement))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
