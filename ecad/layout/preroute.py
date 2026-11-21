#!/usr/bin/env python3
"""Scripted copper for redline_ecu: ground vias at the SMD ground pads (host Python).

    python3 layout/preroute.py <placed.json> <preroute.json>

L2 (In1.Cu) is one solid ground plane, so every SMD ground pad gets its own via to it:

  * exposed pads and other large ground pads (both sides >= EP_MIN) get a grid of thermal
    vias inside the pad (DRV8873 PowerPAD, buck EP, TO-252 / SOT-223 tabs on ground);
  * every other ground pad gets one via just outside the pad, on a short stub, on the side
    away from its own footprint's centre (outwards from an IC, so the via does not sit in the
    pin field), at the first spot that keeps clearance to all foreign copper;
  * the ground pin of a low-side switch (injector, VVT, pump returns) gets two.

The maze router (autoroute.py) treats these vias as fixed and joins any ground pad that did
not get one through the plane.  Output: tracks and vias (mm, board-centre coordinates) for
route.py apply.  Deterministic.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from board import CLR, HOLE_KEEPOUT_R, corner_holes, inside_outline  # noqa: E402

VIA_D, VIA_DRILL = 0.6, 0.3
VIA_PITCH = 0.75            # via centre to via centre (holes 0.45 mm apart, rule 0.25)
EP_MIN = 1.2                # a ground pad this wide both ways takes vias in the pad
EP_PITCH = 1.2
SPACE = CLR + 0.05          # via / stub to foreign copper, with margin for the router's grid
EDGE = 0.5
STUB_W = 0.3
DOUBLE = ("OUT.Q2.X", "OUT.Q3.X", "OUT.Q4.X", "OUT.Q5.X", "OUT.Q6.X", "OUT.Q7.X", "KL.Q8.X")


def box_dist(b, x, y):
    dx = max(b[0] - x, 0, x - b[2])
    dy = max(b[1] - y, 0, y - b[3])
    return math.hypot(dx, dy)


def seg_box_dist(a, b, box, steps=12):
    return min(box_dist(box, a[0] + (b[0] - a[0]) * t / steps, a[1] + (b[1] - a[1]) * t / steps) for t in range(steps + 1))


def main(src, dst):
    d = json.load(open(src))
    foreign = []        # (box, net, is_hole) of all SMD/TH pads
    for f in d["footprints"]:
        for p in f["pads"]:
            if not p["layers"] and not p["drill"]:
                continue    # paste-only aperture
            foreign.append((p["box"], p["net"], p["drill"] > 0 and not p["layers"] or p["npth"]))
    vias, tracks = [], []

    def via_ok(x, y, own_box=None):
        if not inside_outline(x, y, EDGE + VIA_D / 2):
            return False
        for hx, hy in corner_holes():
            if math.hypot(x - hx, y - hy) < HOLE_KEEPOUT_R + VIA_D / 2 + 0.1:
                return False
        for box, net, hole in foreign:
            dist = box_dist(box, x, y) - VIA_D / 2
            if net == "GND" and not hole:
                # no via in or against another ground pad (solder wicking), except its own pad
                if box is not own_box and dist < 0.1:
                    return False
            elif dist < (0.3 if hole else SPACE):
                return False
        for v in vias:
            if math.hypot(v["at"][0] - x, v["at"][1] - y) < VIA_PITCH:
                return False
        return True

    def stub_ok(a, b, own_box):
        for box, net, hole in foreign:
            if net != "GND" and seg_box_dist(a, b, box) < SPACE + STUB_W / 2:
                return False
        return True

    count_ep = count_out = 0
    missed = []
    for f in sorted(d["footprints"], key=lambda f: f["ref"]):
        cx, cy = f["at"]
        for p in f["pads"]:
            if p["net"] != "GND" or p["drill"] > 0 or "F" not in p["layers"]:
                continue
            b = p["box"]
            w, h = b[2] - b[0], b[3] - b[1]
            if w >= EP_MIN and h >= EP_MIN:
                nx, ny = max(1, int((w - 0.6) / EP_PITCH) + 1), max(1, int((h - 0.6) / EP_PITCH) + 1)
                x0, y0 = p["at"][0] - (nx - 1) * EP_PITCH / 2, p["at"][1] - (ny - 1) * EP_PITCH / 2
                placed = 0
                for i in range(nx):
                    for j in range(ny):
                        x, y = round(x0 + i * EP_PITCH, 3), round(y0 + j * EP_PITCH, 3)
                        if box_dist((b[0] + 0.3, b[1] + 0.3, b[2] - 0.3, b[3] - 0.3), x, y) < 1e-9 and via_ok(x, y, b):
                            vias.append({"net": "GND", "at": [x, y], "dia": VIA_D, "drill": VIA_DRILL})
                            placed += 1
                if placed:
                    count_ep += placed
                    continue
            want = 2 if f["path"] in DOUBLE else 1
            px, py = p["at"]
            # directions, best first: away from the footprint centre, then the others
            dirs = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            dirs.sort(key=lambda u: -((px - cx) * u[0] + (py - cy) * u[1]) if (abs(px - cx) + abs(py - cy)) > 0.3 else 0)
            got = 0
            for ux, uy in dirs:
                half = (w if ux else h) / 2
                for extra in (0.0, 0.2, 0.4, 0.7, 1.0, 1.5):
                    for side in (0.0, 0.5, -0.5):
                        dist = half + 0.15 + VIA_D / 2 + extra
                        x = round(px + ux * dist + (side * (w / 2) if uy else 0), 2)
                        y = round(py + uy * dist + (side * (h / 2) if ux else 0), 2)
                        a = (px + (side * (w / 2) if uy else 0), py + (side * (h / 2) if ux else 0))
                        if via_ok(x, y, b) and stub_ok(a, (x, y), b):
                            vias.append({"net": "GND", "at": [x, y], "dia": VIA_D, "drill": VIA_DRILL})
                            tracks.append({"net": "GND", "layer": "F.Cu", "start": [round(a[0], 3), round(a[1], 3)],
                                           "end": [x, y], "width": STUB_W})
                            got += 1
                            break
                    if got and (got >= want or True):
                        break
                if got >= want:
                    break
            count_out += got
            if not got:
                missed.append("%s.%s" % (f["ref"], p["number"]))
    json.dump({"tracks": tracks, "vias": vias}, open(dst, "w"), indent=0)
    print("preroute: %d ground vias in exposed pads, %d beside pads; %d ground pads left to the router%s" % (
        count_ep, count_out, len(missed), (": " + " ".join(missed)) if missed else ""))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
