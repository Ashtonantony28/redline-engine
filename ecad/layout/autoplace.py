#!/usr/bin/env python3
"""Placement for redline_ecu: the fixed parts from a table, every other part clustered by its nets.

    python3 layout/autoplace.py <board.json> <placement.json>

Input is a geometry dump of the board in any state (layout/dump.py); only the footprints' pad
offsets and courtyards are used.  Output maps each footprint path to (x, y, rotation, side) in
board-centre millimetres; layout/place.py applies it under KiCad's Python.

Floor plan (docs/ecad.md, layout section).  Both Micro-Fit headers face out of the bottom
edge: J1 (battery, actuators, CAN) on the right, J2 (sensors) on the left, so that the power
stage and the actuator drivers form the right half and the analogue front ends the left half,
with the MCU on top between them:

    +--------------------------------------------------------------+
    | Tag-Connect   USB-C             | buck + 3V3 LDO               |
    |  level shifters    MCU          | bulk caps                    |
    |  CJ125   TPIC8101               | CAN    heater FET   ideal    |
    |  LM1815  trackers  low-side x6  | gate drivers  ETB   diode    |
    |  [ J2 sensors     ]    [ J1 battery / actuators ]   TVS        |
    +--------------------------------------------------------------+

Every part not in FIXED is placed greedily: parts with the most connections to placed parts
first, each at the free spot (courtyard + GAP, all four rotations) nearest the mean of the pads
it connects to, scored by the distance to those pads.  Supply-only parts (decoupling) are drawn
to the IC their name refers to (C_U8_VDD1 -> U8 of the same module) or, for the MCU, to the
supply pin in MCU_DECAPS.  Deterministic.
"""
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from board import CONN_Y, H, HOLE_KEEPOUT_R, W, corner_holes, inside_outline  # noqa: E402

GRID = 0.1
GAP = 0.3                  # courtyard to courtyard, on top of the courtyards own margin
EDGE_KEEP = 1.0            # courtyard to board edge (connectors and USB-C excepted)
SEARCH = 30.0              # mm: furthest a part may land from its target
GLOBAL = {"GND", "3V3", "5V", "VBAT", "VDDA"}

# path -> (x, y, rotation, side).  Connectors are aligned to the edge in align_edge().
FIXED = {
    "CONN.J1.X": (30.7, None, 180, "F"),
    "CONN.J2.X": (-30.7, None, 180, "F"),
    "COM.J3.X": (6.0, None, 180, "F"),
    "COM.J4.X": (-34.0, -38.0, 0, "F"),
    "MCU.U1.X": (-10.0, -20.0, 0, "F"),
    # --- power input, right column: TVS at the connector, then the ideal diode, bulk, buck, LDO ---
    "PWR.D1.D": (53.0, 22.0, 90, "F"),
    "PWR.Q1.X": (51.0, 9.0, 0, "F"),
    "PWR.U2.X": (43.5, 9.5, 0, "F"),
    "PWR.C1.X": (51.5, -3.5, 0, "F"),
    "PWR.C2.X": (51.5, -14.5, 0, "F"),
    "PWR.U3.X": (40.0, -24.0, 0, "F"),
    "PWR.U4.X": (44.0, -36.0, 0, "F"),
    # --- actuator drivers above J1: injector switches with the drain tab towards the header ---
    "OUT.Q2.X": (32.0, 23.5, 270, "F"),
    "OUT.Q3.X": (22.0, 23.5, 270, "F"),
    "OUT.Q4.X": (12.0, 23.5, 270, "F"),
    "OUT.Q5.X": (2.0, 23.5, 270, "F"),
    "OUT.Q6.X": (38.0, 9.0, 270, "F"),
    "OUT.Q7.X": (4.0, 10.0, 270, "F"),
    "OUT.U8.X": (26.0, 10.0, 0, "F"),
    "OUT.U9.X": (15.0, 10.0, 0, "F"),
    "ETB.U7.X": (41.0, 19.5, 270, "F"),   # motor and VM pins face J1
    "KL.Q8.X": (28.0, -3.0, 270, "F"),
    "KL.U10.X": (20.0, -3.0, 0, "F"),
    "COM.U15.X": (10.0, -4.0, 0, "F"),
    # --- analogue front ends above J2 ---
    "SENS.U11.X": (-18.0, 21.0, 0, "F"),
    "PWR.U5.X": (-8.0, 14.0, 0, "F"),
    "PWR.U6.X": (-8.0, 6.0, 0, "F"),
    "KL.U12.X": (-31.0, 9.0, 0, "F"),
    "KL.U13.X": (-46.0, 15.0, 0, "F"),
    "KL.U14.X": (-30.0, -6.0, 0, "F"),
    "KL.U22.X": (-40.0, -4.0, 0, "F"),
}
for _i, (_x, _y) in enumerate(corner_holes()):
    FIXED["H%d" % (_i + 1)] = (_x, _y, 0, "F")

# Escape corridors kept free of passives: (footprint, pins, length in mm away from the part).
# The DRV8873 motor outputs need a straight 1 mm run to J1.13-16.
KEEP_CLEAR = [("ETB.U7.X", ("14", "15"), 6.0), ("ETB.U7.X", ("18", "19"), 6.0)]

# Decoupling at the MCU pins: capacitor path -> MCU pad (placed like any part, drawn to that pad)
MCU_DECAPS = {
    "MCU.C_VDD_6.C": "6", "MCU.C_VDD_11.C": "11", "MCU.C_VDDA_100N.C": "21", "MCU.C_VDDA_1U.C": "20",
    "MCU.C_VDD_27.C": "27", "MCU.C_VCAP_48.C": "48", "MCU.C_VDD_50.C": "50", "MCU.C_VCAP_73.C": "73",
    "MCU.C_VDD_75.C": "75", "MCU.C_VDD_100.C": "100", "MCU.C_VDD_BULK.C": "75",
}
# parts drawn to a fixed spot (mm) instead of their nets
SPOT = {
    "TP1.TP": (46.0, 2.0),     # VBAT
    "TP2.TP": (34.0, -32.0),   # 5V
    "TP3.TP": (38.0, -40.0),   # 3V3
    "TP4.TP": (-22.0, -10.0),  # VDDA
    "TP7.TP": (36.0, -40.0),   # GND by the supply test points
    "TP8.TP": (-12.0, 2.0),    # GND by the sensor supplies
    "MCU.LED_STATUS.LED": (-24.0, -42.0),
    "MCU.LED_FAULT.LED": (-20.0, -42.0),
    "PWR.LED_PWR.LED": (30.0, -42.0),
}


def rot_pt(x, y, rot):
    """KiCad orientation (degrees, counter-clockwise on screen with y down)."""
    a = math.radians(rot)
    c, s = round(math.cos(a)), round(math.sin(a))
    return (x * c + y * s, -x * s + y * c)


class Part:
    def __init__(self, f):
        self.ref, self.path = f["ref"], f["path"]
        x0, y0 = f["at"]
        r = -f["rot"]
        # pad offsets and courtyard at orientation 0
        self.pads = []
        for p in f["pads"]:
            dx, dy = rot_pt(p["at"][0] - x0, p["at"][1] - y0, r)
            self.pads.append((p["number"], p["net"], dx, dy))
        cy = f["courtyard"] or [min(p["box"][0] for p in f["pads"]), min(p["box"][1] for p in f["pads"]),
                                max(p["box"][2] for p in f["pads"]), max(p["box"][3] for p in f["pads"])]
        pts = [rot_pt(cx - x0, cyy - y0, r) for cx in (cy[0], cy[2]) for cyy in (cy[1], cy[3])]
        self.box0 = (min(q[0] for q in pts), min(q[1] for q in pts), max(q[0] for q in pts), max(q[1] for q in pts))
        self.nets = {p[1] for p in self.pads if p[1] and not p[1].startswith("unconnected")}
        self.at = None

    def box(self, x, y, rot):
        pts = [rot_pt(a, b, rot) for a in (self.box0[0], self.box0[2]) for b in (self.box0[1], self.box0[3])]
        return (x + min(q[0] for q in pts), y + min(q[1] for q in pts), x + max(q[0] for q in pts), y + max(q[1] for q in pts))

    def pad_xy(self, x, y, rot):
        out = []
        for num, net, dx, dy in self.pads:
            px, py = rot_pt(dx, dy, rot)
            out.append((num, net, x + px, y + py))
        return out


OCC = None


class Occupancy:
    def __init__(self):
        self.nx, self.ny = int(W / GRID) + 1, int(H / GRID) + 1
        self.rows = [bytearray(self.nx) for _ in range(self.ny)]

    def _span(self, box):
        x0 = max(0, int(math.floor((box[0] + W / 2) / GRID)))
        x1 = min(self.nx, int(math.ceil((box[2] + W / 2) / GRID)) + 1)
        y0 = max(0, int(math.floor((box[1] + H / 2) / GRID)))
        y1 = min(self.ny, int(math.ceil((box[3] + H / 2) / GRID)) + 1)
        return x0, x1, y0, y1

    def free(self, box):
        x0, x1, y0, y1 = self._span(box)
        zero = bytes(x1 - x0)
        return all(self.rows[y][x0:x1] == zero for y in range(y0, y1))

    def take(self, box):
        x0, x1, y0, y1 = self._span(box)
        one = b"\x01" * (x1 - x0)
        for y in range(y0, y1):
            self.rows[y][x0:x1] = one


def legal(box, edge_ok=False):
    m = 0.0 if edge_ok else EDGE_KEEP
    for x in (box[0], box[2]):
        for y in (box[1], box[3]):
            if not inside_outline(x, y, m):
                return False
    for hx, hy in corner_holes():
        dx = max(box[0] - hx, 0, hx - box[2])
        dy = max(box[1] - hy, 0, hy - box[3])
        if math.hypot(dx, dy) < HOLE_KEEPOUT_R + 0.5:
            return False
    return True


def search(p, tx, ty, rots, score):
    """Rings of growing radius around (tx, ty); the best legal spot of the first ring that has one."""
    steps = int(SEARCH / GRID)
    for r in range(0, steps, 2):
        if r == 0:
            ring = [(0, 0)]
        else:
            ring = [(i, -r) for i in range(-r, r + 1, 2)] + [(i, r) for i in range(-r, r + 1, 2)] + \
                   [(-r, j) for j in range(-r + 2, r, 2)] + [(r, j) for j in range(-r + 2, r, 2)]
        cands = []
        for i, j in ring:
            x, y = round((tx + i * GRID) / GRID) * GRID, round((ty + j * GRID) / GRID) * GRID
            for rot in rots:
                b = p.box(x, y, rot)
                if legal(b) and OCC.free(b):
                    cands.append((score(x, y, rot) + 0.3 * math.hypot(x - tx, y - ty), rot, x, y))
        if cands:
            return min(cands)
    return None


def nearest(p, tx, ty, rots, score):
    best = search(p, tx, ty, rots, score)
    assert best, "no room for %s near (%.1f, %.1f)" % (p.path, tx, ty)
    if abs(best[2] - tx) > 0.05 or abs(best[3] - ty) > 0.05:
        print("  %s moved from (%.1f, %.1f) to (%.1f, %.1f)" % (p.path, tx, ty, best[2], best[3]))
    return best[2], best[3]


def inflate(b, g):
    return (b[0] - g, b[1] - g, b[2] + g, b[3] + g)


def main(src, dst):
    global OCC
    OCC = Occupancy()
    data = json.load(open(src))
    parts = {f["path"]: Part(f) for f in data["footprints"]}
    occ = OCC
    placed = {}

    def put(path, x, y, rot, edge_ok=False):
        p = parts[path]
        x, y = round(x / GRID) * GRID, round(y / GRID) * GRID
        b = p.box(x, y, rot)
        assert path.startswith("H") or legal(b, edge_ok), "%s does not fit at (%.1f, %.1f)" % (path, x, y)
        assert occ.free(b), "%s overlaps a placed part at (%.1f, %.1f)" % (path, x, y)
        occ.take(inflate(b, GAP / 2))
        p.at = (round(x, 3), round(y, 3), rot)
        placed[path] = p

    for path, (x, y, rot, side) in FIXED.items():
        p = parts[path]
        if y is None:   # connector: courtyard front on the edge
            b = p.box(0, 0, rot)   # x in the table is the courtyard centre
            y = (CONN_Y - b[3]) if path.startswith("CONN.") else (-H / 2 - b[1])
            put(path, x - (b[0] + b[2]) / 2, y, rot, edge_ok=True)
        elif path.startswith("H"):
            put(path, x, y, rot, edge_ok=True)
        else:   # the nearest free spot to the table's, in the table's rotation
            x, y = nearest(parts[path], x, y, (rot,), lambda *a: 0.0)
            put(path, x, y, rot)

    for path, pins, length in KEEP_CLEAR:
        p = parts[path]
        pts = [(x, y) for num, net, x, y in p.pad_xy(*p.at) if num in pins]
        mx, my = sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts)
        dx, dy = mx - p.at[0], my - p.at[1]
        x0, y0 = min(q[0] for q in pts) - 0.5, min(q[1] for q in pts) - 0.5
        x1, y1 = max(q[0] for q in pts) + 0.5, max(q[1] for q in pts) + 0.5
        if abs(dx) > abs(dy):
            x0, x1 = (x0, x1 + length) if dx > 0 else (x0 - length, x1)
        else:
            y0, y1 = (y0, y1 + length) if dy > 0 else (y0 - length, y1)
        occ.take((x0, y0, x1, y1))

    pad_at = {}

    def refresh(path):
        p = parts[path]
        for num, net, x, y in p.pad_xy(*p.at):
            pad_at.setdefault(net, []).append((path, num, x, y))

    for path in placed:
        refresh(path)

    def anchors(path):
        """Points this part wants to be near: (x, y, weight)."""
        p = parts[path]
        if path in SPOT:
            return [(SPOT[path][0], SPOT[path][1], 1.0, None)]
        if path in MCU_DECAPS:
            mcu = parts["MCU.U1.X"]
            pin = next((x, y) for num, net, x, y in mcu.pad_xy(*mcu.at) if num == MCU_DECAPS[path])
            return [(pin[0], pin[1], 1.0, next(net for num, net, x, y in mcu.pad_xy(*mcu.at) if num == MCU_DECAPS[path]))]
        pts = []
        # a capacitor named after an IC (C_U5_OUT, C_U8_VDD1) sits at that IC's pins on its nets
        m = re.match(r"\w+\.C_(U\d+)_", path)
        if m:
            target = "%s.%s.X" % (path.split(".")[0], m.group(1))
            if target in placed:
                q = parts[target]
                pts = [(x, y, 1.0, net) for num, net, x, y in q.pad_xy(*q.at) if net in p.nets - {"GND"}]
                if pts:
                    return pts
        for net in p.nets - GLOBAL:
            for other, num, x, y in pad_at.get(net, []):
                if other != path:
                    pts.append((x, y, 1.0, net))
        if pts:
            return pts
        # supply-only part: the IC its name points at, else the nearest supply pin of the module
        module = path.split(".")[0]
        m = re.search(r"_(U\d+|Q\d+|D\d+)", path)
        if m:
            target = "%s.%s.X" % (module, m.group(1))
            if target in placed:
                q = parts[target]
                return [(x, y, 1.0, net) for num, net, x, y in q.pad_xy(*q.at) if net in p.nets - {"GND"}] or [(q.at[0], q.at[1], 1.0, None)]
        for net in sorted(p.nets - {"GND"}):
            for other, num, x, y in pad_at.get(net, []):
                if other.startswith(module + ".") and other != path:
                    pts.append((x, y, 1.0, net))
        return pts

    def links(path):
        p = parts[path]
        if path in SPOT or path in MCU_DECAPS:
            return 100
        m = re.match(r"\w+\.C_(U\d+)_", path)
        if m and "%s.%s.X" % (path.split(".")[0], m.group(1)) in placed:
            return 50      # decoupling first, while the IC's pins still have room round them
        n = sum(1 for net in p.nets - GLOBAL for other, _, _, _ in pad_at.get(net, []) if other != path)
        if n == 0 and anchors(path):
            return 0.5
        return n

    todo = sorted(p for p in parts if p not in placed)
    while todo:
        todo.sort(key=lambda q: (-links(q), q))
        path = todo.pop(0)
        p = parts[path]
        pts = anchors(path)
        if not pts:
            pts = [(0.0, 0.0, 1.0, None)]
        tx = sum(a[0] * a[2] for a in pts) / sum(a[2] for a in pts)
        ty = sum(a[1] * a[2] for a in pts) / sum(a[2] for a in pts)

        def score(x, y, rot):
            total = 0.0
            mine = p.pad_xy(x, y, rot)
            for ax, ay, w, anet in pts:
                same = [(px, py) for _, net, px, py in mine if net == anet] or [(px, py) for _, net, px, py in mine]
                total += w * min(math.hypot(px - ax, py - ay) for px, py in same)
            return total

        best = search(p, tx, ty, (0, 90, 180, 270), score)
        assert best, "no room for %s near (%.1f, %.1f)" % (path, tx, ty)
        _, rot, x, y = best
        put(path, x, y, rot)
        refresh(path)
    out = {path: [p.at[0], p.at[1], p.at[2], "F"] for path, p in sorted(parts.items())}
    json.dump(out, open(dst, "w"), indent=0)
    print("autoplace: %d parts (%d fixed)" % (len(out), len(FIXED)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
