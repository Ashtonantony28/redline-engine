#!/usr/bin/env python3
"""Grid maze router for redline_ecu: routes everything preroute.py left open.

    python3 layout/autoroute.py layout/build/pre.json layout/build/autoroute.json

Adapted from drone-design ecad/layout/autoroute.py (same method, generalised to a W x H board,
with the drone's IMU island, DShot, QSPI-matching and L3-plane code taken out and a search
window per net so that a 120 x 90 mm board stays fast).

Input is the geometry dump of the placed and pre-routed board; output is a list of tracks and
vias (mm, board-centre coordinates) for route.py to add.

Method: A* on a 0.1 mm grid, 8-direction moves plus through vias, over L1, L3 (In2.Cu, at a
premium) and L4; L2 stays an unbroken ground plane.  Nets are routed by negotiated congestion
(PathFinder): each net is routed as if alone and pays for cells inside another net's keep-off
halo; the price rises until nets stop overlapping.  What still overlaps is taken out and put
back strictly, with a targeted rip-up for the last connections.  Clearances are baked into the
occupancy maps (copper inflated by clearance + half the track width), so the output is clean by
construction; KiCad's DRC stays the judge (layout/run.sh drc).

Track classes (layout/board.py): signals 0.2 mm (T map), supplies 0.3 mm (Q map), high-current
nets 1.0 mm (W map) wherever the 1.0 mm track fits, necking to 0.3 mm elsewhere at a cost.

Pure Python + numpy, deterministic.
"""
import heapq
import json
import math
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from board import CLR, CORNER_R, H, HOLE_KEEPOUT_R, NARROW_POWER, PLANE_NETS, POWER_NETS, TRACK, W, WIDE, WIDE_NETS, corner_holes  # noqa: E402

GRID = 0.1
NX = int(round(W / GRID)) + 1
NY = int(round(H / GRID)) + 1
NN = NX * NY
EDGE = 0.45                 # copper to edge (rule 0.4)
HOLE_CLR = 0.27             # hole to foreign copper (rule 0.25)
VIA_D, VIA_DRILL = 0.6, 0.3
VIA_PITCH = 0.75
VIA_EXTRA = 0.02
VIA_HOLE_CLR = 0.15
FIXED_EXTRA = 0.02          # pads and scripted copper are off the grid: a diagonal step may cut a corner
VIA_PAD_GAP = 0.05          # via edge to any SMD pad: no via-in-pad
POWER_CLASS = set(POWER_NETS) | set(PLANE_NETS)
WIDE_SET = set(WIDE_NETS)
WIDTHS = {}
PLANES = {n: None for n in PLANE_NETS}
NEAR_PAD = 1.5              # a high-current net may neck down to 0.3 mm this close to its own pads
NARROW_COST = 8             # per step of a high-current net on a 0.3 mm neck away from its pads
FIRST = list(WIDE_NETS) + ["PWR.BUCK_SW", "PWR.BUCK_BOOT", "PWR.BUCK_FB", "MCU.HSE_IN", "MCU.HSE_OUT", "USB_P", "USB_N"]
ORTH, DIAG, VIA_COST = 10, 14, 300
WIDE_VIA_COST = 1500        # high-current nets stay on one layer where they can
PROX_COST = 6
PROX_DIST = 0.45
LANE_LEN = 1.5
LANE_HALF = 0.2
LANE_COST = 40
LANE_VIA_LEN = 1.0
LANE_VIA_HALF = 0.43
PRES_START = 10
PRES_GROWTH = 1.25
PRES_MAX = 3000
HIST_FRACTION = 0.25
MAX_ITERATIONS = int(os.environ.get("AUTOROUTE_ITERATIONS", 60))
RIPUP_ROUNDS = 6
RENEG_ITERATIONS = 30
RENEG_SCALE = (1, 0.3, 3, 0.1)
RENEG_PRES = 150
PATIENCE = 15
LAYERS = ("F", "I", "B")    # I = In2.Cu; In1.Cu is the ground plane
INNER_COST = 12
OUTER_ONLY = ("MCU.HSE", "PWR.BUCK")
MARGIN = 80                 # search window: cells round the bounding box of a net's pads
ASTAR_LIMIT = 1500000


def cellx(v):
    return int(round((v + W / 2) / GRID))


def celly(v):
    return int(round((v + H / 2) / GRID))


def coordx(i):
    return round(i * GRID - W / 2, 4)


def coordy(i):
    return round(i * GRID - H / 2, 4)


class Maps:
    """Occupancy maps.  Each holds, per cell, 0 (free), a net id (inside the keep-off halo of
    that net's copper only) or -1 (blocked for everyone)."""

    def __init__(self):
        self.X, self.Y = np.meshgrid(np.arange(NX) * GRID - W / 2, np.arange(NY) * GRID - H / 2)
        # keep-off maps for: a signal track centre (T), a 0.3 mm supply track (Q), a wide track or a
        # power-class via (W), a via (V)
        self.T = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.Q = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.W = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.V = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.deep = {L: np.zeros((NY, NX), bool) for L in LAYERS}
        self.viaban = np.zeros((NY, NX), bool)
        self.comp = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}

    def window(self, box, r):
        x0, y0 = max(0, cellx(box[0] - r) - 1), max(0, celly(box[1] - r) - 1)
        x1, y1 = min(NX, cellx(box[2] + r) + 2), min(NY, celly(box[3] + r) + 2)
        return slice(y0, y1), slice(x0, x1)

    def dist_rect(self, win, box):
        X, Y = self.X[win], self.Y[win]
        dx = np.maximum(np.maximum(box[0] - X, X - box[2]), 0)
        dy = np.maximum(np.maximum(box[1] - Y, Y - box[3]), 0)
        return np.hypot(dx, dy)

    def dist_seg(self, win, a, b):
        X, Y = self.X[win], self.Y[win]
        vx, vy = b[0] - a[0], b[1] - a[1]
        L2 = vx * vx + vy * vy
        if L2 == 0:
            return np.hypot(X - a[0], Y - a[1])
        t = np.clip(((X - a[0]) * vx + (Y - a[1]) * vy) / L2, 0, 1)
        return np.hypot(X - a[0] - t * vx, Y - a[1] - t * vy)

    @staticmethod
    def mark(arr, win, mask, nid):
        sub = arr[win]
        sub[mask & (sub != 0) & (sub != nid)] = -1
        sub[mask & (sub == 0)] = nid


class Router:
    def __init__(self, data):
        self.data = data
        self.m = Maps()
        self.net_id = {}
        self.net_name = {}
        self.parent = [0]
        self.comp_net = [0]
        self.anchor = {}
        self.plane_comp = {}
        self.fake = 0
        self.prox = {L: np.zeros((NY, NX), np.uint8) for L in LAYERS}
        self.lane = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.lane_via = np.zeros((NY, NX), np.int32)
        self.pad_boxes = {}
        self.build()

    def nid(self, net):
        if net not in self.net_id:
            self.net_id[net] = len(self.net_id) + 1
            self.net_name[self.net_id[net]] = net
        return self.net_id[net]

    def new_comp(self, net):
        self.parent.append(len(self.parent))
        self.comp_net.append(self.nid(net))
        return len(self.parent) - 1

    def find(self, c):
        while self.parent[c] != c:
            self.parent[c] = self.parent[self.parent[c]]
            c = self.parent[c]
        return c

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)

    def set_comp(self, layer, win, mask, c):
        sub = self.m.comp[layer][win]
        for other in np.unique(sub[mask]):
            if other and self.comp_net[other] == self.comp_net[c]:
                self.union(c, int(other))
        sub[mask] = c

    def halo(self, layer, kind, geom, net, extra=0.0):
        """Keep-off halo of a piece of copper as (map name, window, mask) for the four maps."""
        c = (HOLE_CLR if net is None else CLR) + extra
        m = self.m
        for name, r in (("T", c + TRACK / 2), ("Q", c + NARROW_POWER / 2), ("W", max(c, VIA_HOLE_CLR + extra) + max(WIDE, VIA_D) / 2),
                        ("V", max(c, VIA_HOLE_CLR + extra) + VIA_D / 2)):
            if kind == "rect":
                win = m.window(geom, r)
                d = m.dist_rect(win, geom)
            elif kind == "seg":
                a, b, w = geom
                win = m.window((min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])), r + w / 2)
                d = m.dist_seg(win, a, b) - w / 2
            else:
                p, dia = geom
                win = m.window((p[0], p[1], p[0], p[1]), r + dia / 2)
                d = m.dist_seg(win, p, p) - dia / 2
            yield name, win, d < r - 1e-9

    def stamp_halo(self, layer, kind, geom, net, extra=0.0):
        nid = -1 if net is None else self.nid(net)
        maps = {"T": self.m.T, "Q": self.m.Q, "W": self.m.W, "V": self.m.V}
        for name, win, mask in self.halo(layer, kind, geom, net, extra):
            self.m.mark(maps[name][layer], win, mask, nid)

    def add_pad(self, p):
        layers = LAYERS if p["drill"] > 0 else [L for L in ("F", "B") if L in p["layers"]]
        box = p["box"]
        if p["npth"] or not p["net"]:
            for L in layers:
                self.fake += 1
                self.stamp_halo(L, "rect", box, None if p["drill"] > 0 else "~nc%d" % self.fake, extra=FIXED_EXTRA)
        else:
            for L in layers:
                self.stamp_halo(L, "rect", box, p["net"], extra=FIXED_EXTRA)
            self.pad_boxes.setdefault(p["net"], []).append(box)
        r = VIA_D / 2 + (HOLE_CLR if p["drill"] > 0 else VIA_PAD_GAP)
        win = self.m.window(box, r)
        self.m.viaban[win] |= self.m.dist_rect(win, box) < r
        for L in layers:
            win = self.m.window(box, PROX_DIST)
            self.prox[L][win][self.m.dist_rect(win, box) < PROX_DIST] = PROX_COST
        if not p["net"] or p["npth"]:
            return None
        c = self.new_comp(p["net"])
        for L in layers:
            win = self.m.window(box, 0)
            if p.get("round"):
                rad = min(box[2] - box[0], box[3] - box[1]) / 2
                core = (box[0] + rad, box[1] + rad, box[2] - rad, box[3] - rad)
                dcore = self.m.dist_rect(win, core)
                mask = dcore <= rad - 0.05 + 1e-9
            else:
                rr = min(0.25, 0.25 * min(box[2] - box[0], box[3] - box[1]))
                core = (box[0] + rr, box[1] + rr, box[2] - rr, box[3] - rr)
                mask = self.m.dist_rect(win, core) <= rr - 0.05 + 1e-9
            if not mask.any():
                d = np.hypot(self.m.X[win] - p["at"][0], self.m.Y[win] - p["at"][1])
                mask = d <= d.min() + 1e-9
            self.set_comp(L, win, mask, c)
            half = TRACK / 2
            if p.get("round"):
                self.m.deep[L][win] |= dcore <= rad - half + 1e-9
            elif box[2] - box[0] >= 2 * half and box[3] - box[1] >= 2 * half:
                self.m.deep[L][win] |= self.m.dist_rect(win, (box[0] + half, box[1] + half, box[2] - half, box[3] - half)) < 1e-9
        if p["drill"] > 0 and p["net"] in PLANES:
            self.union(c, self.plane(p["net"]))
        return c

    def plane(self, net):
        if net not in self.plane_comp:
            self.plane_comp[net] = self.new_comp(net)
        return self.plane_comp[net]

    def add_track(self, t, exact=True):
        L = t["layer"][0]
        a, b, w = t["start"], t["end"], t["width"]
        self.stamp_halo(L, "seg", (a, b, w), t["net"], extra=FIXED_EXTRA)
        c = self.new_comp(t["net"])
        win = self.m.window((min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1])), 0.1)
        d = self.m.dist_seg(win, a, b)
        mask = d <= 0.0708
        self.set_comp(L, win, mask, c)
        if exact:
            ys, xs = np.nonzero(mask)
            vx, vy = b[0] - a[0], b[1] - a[1]
            L2 = vx * vx + vy * vy
            for yy, xx in zip(ys, xs):
                gy, gx = yy + win[0].start, xx + win[1].start
                px, py = coordx(gx), coordy(gy)
                tt = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - a[0]) * vx + (py - a[1]) * vy) / L2))
                self.anchor[(L, gy, gx)] = (round(a[0] + tt * vx, 4), round(a[1] + tt * vy, 4))
        return c

    def add_via(self, v, exact=True):
        p = v["at"]
        c = self.new_comp(v["net"])
        for L in LAYERS:
            self.stamp_halo(L, "via", (p, v["dia"]), v["net"], extra=VIA_EXTRA)
            win = self.m.window((p[0], p[1], p[0], p[1]), 0.1)
            d = np.hypot(self.m.X[win] - p[0], self.m.Y[win] - p[1])
            mask = d <= max(0.0708, d.min() + 1e-9)
            self.set_comp(L, win, mask, c)
            if exact:
                ys, xs = np.nonzero(mask)
                for yy, xx in zip(ys, xs):
                    self.anchor[(L, yy + win[0].start, xx + win[1].start)] = (p[0], p[1])
        win = self.m.window((p[0], p[1], p[0], p[1]), VIA_PITCH)
        self.m.viaban[win] |= np.hypot(self.m.X[win] - p[0], self.m.Y[win] - p[1]) < VIA_PITCH - 1e-9
        if v["net"] in PLANES:
            self.union(c, self.plane(v["net"]))
        return c

    def build(self):
        m = self.m
        kx, ky = W / 2 - CORNER_R, H / 2 - CORNER_R
        ax, ay = np.abs(m.X), np.abs(m.Y)
        corner = (ax > kx) & (ay > ky)
        hole = np.full((NY, NX), 1e9)
        for hx, hy in corner_holes():
            hole = np.minimum(hole, np.hypot(m.X - hx, m.Y - hy))
        for arrs, half in ((m.T, TRACK / 2), (m.Q, NARROW_POWER / 2), (m.W, WIDE / 2), (m.V, VIA_D / 2)):
            margin = EDGE + half
            bad = (ax > W / 2 - margin + 1e-9) | (ay > H / 2 - margin + 1e-9) | \
                  (corner & (np.hypot(ax - kx, ay - ky) > CORNER_R - margin)) | (hole < HOLE_KEEPOUT_R + half + 0.05)
            for L in LAYERS:
                arrs[L][bad] = -1
        bad = (ax > W / 2 - EDGE - VIA_D / 2) | (ay > H / 2 - EDGE - VIA_D / 2) | (hole < HOLE_KEEPOUT_R + VIA_D / 2 + 0.05)
        m.viaban |= bad
        for f in self.data["footprints"]:
            for p in f["pads"]:
                if p["layers"] or p["drill"]:
                    self.add_pad(p)
        count = {}
        for f in self.data["footprints"]:
            for p in f["pads"]:
                count[p["net"]] = count.get(p["net"], 0) + 1
        for f in self.data["footprints"]:
            for p in f["pads"]:
                b = p["box"]
                w, h = b[2] - b[0], b[3] - b[1]
                if p["drill"] or not p["layers"] or not p["net"] or count[p["net"]] < 2 or min(w, h) > 0.45 or max(w, h) < 2 * min(w, h):
                    continue
                nid = self.nid(p["net"])
                x, y = p["at"]
                for length, half, arrs in ((LANE_LEN, LANE_HALF, [self.lane[L] for L in ("F", "B") if L in p["layers"]]),
                                           (LANE_VIA_LEN, LANE_VIA_HALF, [self.lane_via])):
                    if w > h:
                        boxes = [(b[0] - length, y - half, b[0], y + half), (b[2], y - half, b[2] + length, y + half)]
                    else:
                        boxes = [(x - half, b[1] - length, x + half, b[1]), (x - half, b[3], x + half, b[3] + length)]
                    for box in boxes:
                        win = m.window(box, 0)
                        mask = m.dist_rect(win, box) < 1e-9
                        for arr in arrs:
                            m.mark(arr, win, mask, nid)
        for t in self.data["tracks"]:
            self.add_track(t)
        for v in self.data["vias"]:
            self.add_via(v)

    def free_maps(self, net, win):
        """Inside the window: cells a track centre of this net may use on each layer (narrow class),
        the cells where a wide track fits (high-current nets), and via cells."""
        nid = self.nid(net)
        m = self.m
        wide = net in WIDE_SET
        special = net in POWER_CLASS or wide
        narrow = m.Q if special else m.T
        free, wide_ok = {}, {}
        for L in LAYERS:
            sub = narrow[L][win]
            f = (sub == 0) | (sub == nid)
            if L == "I" and net.startswith(OUTER_ONLY):
                f = np.zeros_like(f)
            free[L] = f
            if wide:
                ws = m.W[L][win]
                wide_ok[L] = f & ((ws == 0) | (ws == nid))
        vmap = m.W if special else m.V
        via = ~m.viaban[win]
        for L in LAYERS:
            sub = vmap[L][win]
            via &= (sub == 0) | (sub == nid)
        lv = self.lane_via[win]
        via &= (lv == 0) | (lv == nid)
        return free, (wide_ok if wide else None), via

    def near_pads(self, net, win):
        near = np.zeros((win[0].stop - win[0].start, win[1].stop - win[1].start), bool)
        X, Y = self.m.X[win], self.m.Y[win]
        for b in self.pad_boxes.get(net, []):
            dx = np.maximum(np.maximum(b[0] - X, X - b[2]), 0)
            dy = np.maximum(np.maximum(b[1] - Y, Y - b[3]), 0)
            near |= np.hypot(dx, dy) < NEAR_PAD
        return near

    def components(self, net):
        nid = self.nid(net)
        roots = {}
        for c in range(1, len(self.parent)):
            if self.comp_net[c] == nid:
                roots.setdefault(self.find(c), []).append(c)
        return roots


def astar(sources, tgt, fr, xc, vi, vx, tb, nx, nn, nl, via_cost):
    """Cheapest path from any source cell to any target cell, on a window nx wide (nn cells per
    layer).  Cells are flat indexes layer*nn + y*nx + x; fr = passable, xc = extra cost of
    entering a cell, vi = via allowed (per x,y), vx = extra cost of a via there."""
    tx0, ty0, tx1, ty1 = tb
    best, prev, heap = {}, {}, []
    push, pop = heapq.heappush, heapq.heappop
    for i in sources:
        best[i] = 0
        r = i % nn
        x, y = r % nx, r // nx
        dx = tx0 - x if x < tx0 else x - tx1 if x > tx1 else 0
        dy = ty0 - y if y < ty0 else y - ty1 if y > ty1 else 0
        push(heap, (ORTH * max(dx, dy) + (DIAG - ORTH) * min(dx, dy), 0, i))
    moves = ((1, ORTH), (-1, ORTH), (nx, ORTH), (-nx, ORTH), (nx + 1, DIAG), (nx - 1, DIAG), (-nx + 1, DIAG), (-nx - 1, DIAG))
    others = [[b for b in range(nl) if b != a] for a in range(nl)]
    INF = 1 << 60
    limit = ASTAR_LIMIT
    while heap:
        f, g, i = pop(heap)
        if g > best.get(i, INF):
            continue
        if tgt[i]:
            path = [i]
            while i in prev:
                i = prev[i]
                path.append(i)
            return path[::-1]
        limit -= 1
        if limit < 0:
            return None
        li, r = divmod(i, nn)
        x, y = r % nx, r // nx
        for d, c in moves:
            rj = r + d
            if rj < 0 or rj >= nn:
                continue
            j = i + d
            xj = rj % nx
            if xj - x > 1 or x - xj > 1 or not fr[j]:
                continue
            ng = g + c + xc[j]
            if ng < best.get(j, INF):
                best[j] = ng
                prev[j] = i
                yj = rj // nx
                dx = tx0 - xj if xj < tx0 else xj - tx1 if xj > tx1 else 0
                dy = ty0 - yj if yj < ty0 else yj - ty1 if yj > ty1 else 0
                push(heap, (ng + ORTH * max(dx, dy) + (DIAG - ORTH) * min(dx, dy), ng, j))
        if vi[r]:
            for lj in others[li]:
                j = lj * nn + r
                if fr[j]:
                    ng = g + via_cost + xc[j] + vx[r]
                    if ng < best.get(j, INF):
                        best[j] = ng
                        prev[j] = i
                        dx = tx0 - x if x < tx0 else x - tx1 if x > tx1 else 0
                        dy = ty0 - y if y < ty0 else y - ty1 if y > ty1 else 0
                        push(heap, (ng + ORTH * max(dx, dy) + (DIAG - ORTH) * min(dx, dy), ng, j))
    return None


def _seg_dist(p, a, b):
    vx, vy = b[0] - a[0], b[1] - a[1]
    L2 = vx * vx + vy * vy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * vx + (p[1] - a[1]) * vy) / L2))
    return math.hypot(p[0] - a[0] - t * vx, p[1] - a[1] - t * vy)


def split(i):
    li, r = divmod(i, NN)
    return li, r // NX, r % NX


class Negotiator:
    """Negotiated-congestion routing (PathFinder)."""

    def __init__(self, data):
        self.data = data
        self.fx = Router(data)
        fx = self.fx
        self.cnt = {k: {L: np.zeros((NY, NX), np.int16) for L in LAYERS} for k in "TQWV"}
        self.vp = np.zeros((NY, NX), np.int16)
        self.hist = {L: np.zeros((NY, NX), np.int32) for L in LAYERS}
        self.routes = {}
        nets = sorted({p["net"] for f in data["footprints"] for p in f["pads"] if p["net"]})
        self.nets = [n for n in nets if len(fx.components(n)) > 1]
        self.static = {}
        self.state_file = None

    def prepare(self, net):
        """Per-net constants: the copper components to join, as sets of global cells, and the
        bounding box (cells) of the net."""
        if net in self.static:
            return self.static[net]
        fx = self.fx
        nid = fx.nid(net)
        roots = fx.components(net)
        full = (slice(0, NY), slice(0, NX))
        free, wide_ok, _ = fx.free_maps(net, full)
        comps = []
        plane_root = fx.find(fx.plane_comp[net]) if net in fx.plane_comp else None
        xs, ys = [], []
        for root, ids in sorted(roots.items()):
            cells = set()
            for li, L in enumerate(LAYERS):
                own = np.isin(fx.m.comp[L], ids)
                flat = np.flatnonzero(own & (free[L] | fx.m.deep[L]))
                if len(flat) == 0:
                    flat = np.flatnonzero(own)
                cells.update((flat + li * NN).tolist())
                if len(flat):
                    xs += [int((flat % NX).min()), int((flat % NX).max())]
                    ys += [int((flat // NX).min()), int((flat // NX).max())]
            comps.append((root == plane_root, cells))
        st = {"comps": comps, "nid": nid, "ids": [c for c in range(1, len(fx.parent)) if fx.comp_net[c] == nid], "special": net in POWER_CLASS or net in WIDE_SET, "wide": net in WIDE_SET,
              "bbox": (min(xs), min(ys), max(xs), max(ys)) if xs else (0, 0, NX - 1, NY - 1)}
        self.static[net] = st
        return st

    def counts_at(self, st, L, win, wide_ok):
        if st["wide"]:
            return np.where(wide_ok[L], self.cnt["W"][L][win], self.cnt["Q"][L][win])
        return self.cnt["Q" if st["special"] else "T"][L][win]

    def via_count(self, st, win):
        k = "W" if st["special"] else "V"
        total = self.vp[win].astype(np.int32)
        for L in LAYERS:
            total = total + self.cnt[k][L][win]
        return total

    def route(self, net, pres, hard=False, keep=()):
        st = self.prepare(net)
        left = None
        for margin in (MARGIN, 10 ** 6):
            left, conns = self._route(net, st, pres, hard, keep, margin)
            if left == 0:
                break
        self.routes[net] = {"conns": conns, "stamps": [], "open": left}
        self.stamp(net)
        return left

    def _route(self, net, st, pres, hard, keep, margin):
        fx = self.fx
        bx0, by0, bx1, by1 = st["bbox"]
        x0, y0 = max(0, bx0 - margin), max(0, by0 - margin)
        x1, y1 = min(NX, bx1 + margin + 1), min(NY, by1 + margin + 1)
        win = (slice(y0, y1), slice(x0, x1))
        nx, ny = x1 - x0, y1 - y0
        nn = nx * ny
        nl = len(LAYERS)

        def loc(i):
            li, y, x = split(i)
            if x0 <= x < x1 and y0 <= y < y1:
                return li * nn + (y - y0) * nx + (x - x0)
            return None

        def glob(j):
            li, r = divmod(j, nn)
            return li * NN + (r // nx + y0) * NX + (r % nx + x0)

        free, wide_ok, via = fx.free_maps(net, win)
        nid = st["nid"]
        passable, base = [], []
        near = fx.near_pads(net, win) if st["wide"] else None
        for L in LAYERS:
            passable.append(free[L] | (np.isin(fx.m.comp[L][win], st["ids"]) & fx.m.deep[L][win]))
            cost = fx.prox[L][win].astype(np.int32) + (INNER_COST if L == "I" else 0)
            lane = fx.lane[L][win]
            cost = cost + np.where((lane != 0) & (lane != nid), LANE_COST, 0)
            if st["wide"]:
                cost = cost + np.where(wide_ok[L] | near, 0, NARROW_COST)
            base.append(cost + self.hist[L][win])
        counts = [self.counts_at(st, L, win, wide_ok) for L in LAYERS]
        vcount = self.via_count(st, win)
        if hard:
            fr = b"".join((passable[li] & (counts[li] == 0)).astype(np.uint8).tobytes() for li in range(nl))
            xc = np.concatenate([base[li].ravel() for li in range(nl)]).tolist()
            via = via & (vcount == 0)
            vx = [0] * nn
        else:
            fr = b"".join(passable[li].astype(np.uint8).tobytes() for li in range(nl))
            xc = np.concatenate([(base[li] + pres * np.minimum(counts[li], 3)).ravel() for li in range(nl)]).tolist()
            vx = (pres * np.minimum(vcount, 3)).ravel().tolist()
        own_ban = np.zeros((ny, nx), bool)
        X, Y = fx.m.X[win], fx.m.Y[win]

        def ban(q):
            own_ban[np.hypot(X - q[0], Y - q[1]) < VIA_PITCH - 1e-9] = True

        comps = []
        for is_plane, cells in st["comps"]:
            lc = {j for j in (loc(i) for i in cells) if j is not None}
            comps.append((is_plane, lc))
        conns = []
        for c in keep:
            path = [loc(i) for i in c["cells"]]
            if any(j is None for j in path):
                continue
            ia = next((j for j, (_, cells) in enumerate(comps) if path[0] in cells), None)
            ib = next((j for j, (_, cells) in enumerate(comps) if path[-1] in cells), None)
            idx = sorted({j for j in (ia, ib) if j is not None})
            merged = (any(comps[j][0] for j in idx), set(path).union(*[comps[j][1] for j in idx]))
            comps = [c2 for j, c2 in enumerate(comps) if j not in idx] + [merged]
            conns.append(c)
            for q in c["vias"]:
                ban(q["at"])
        via_cost = WIDE_VIA_COST if st["wide"] else VIA_COST
        left = 0
        while len(comps) > 1:
            done = False
            for k in sorted((k for k in range(len(comps)) if not comps[k][0]), key=lambda k: (len(comps[k][1]), k)):
                src = comps[k][1]
                tgt = bytearray(nl * nn)
                allc = []
                for j, (_, cells) in enumerate(comps):
                    if j != k:
                        for i in cells:
                            tgt[i] = 1
                        allc.extend(cells)
                if not allc or not src:
                    continue
                rr = np.array(allc, dtype=np.int64) % nn
                tb = (int((rr % nx).min()), int((rr // nx).min()), int((rr % nx).max()), int((rr // nx).max()))
                vi = bytearray((via & ~own_ban).astype(np.uint8).tobytes())
                path = None
                for _ in range(6):
                    path = astar(src, tgt, fr, xc, vi, vx, tb, nx, nn, nl, via_cost)
                    if not path:
                        break
                    vcells = [path[i] % nn for i in range(1, len(path)) if path[i] // nn != path[i - 1] // nn]
                    clash = [b for i, b in enumerate(vcells) for a in vcells[:i]
                             if math.hypot((a % nx - b % nx), (a // nx - b // nx)) * GRID < VIA_PITCH]
                    if not clash:
                        break
                    for b in clash:
                        vi[b] = 0
                    path = None
                if not path:
                    continue
                gpath = [glob(j) for j in path]
                t, v = self.geometry(net, st, gpath)
                conns.append({"tracks": t, "vias": v, "cells": gpath})
                for q in v:
                    ban(q["at"])
                end = path[-1]
                j = next(j for j, (_, cells) in enumerate(comps) if j != k and end in cells)
                merged = (comps[k][0] or comps[j][0], comps[k][1] | comps[j][1] | set(path))
                comps = [c for idx, c in enumerate(comps) if idx not in (j, k)] + [merged]
                done = True
                break
            if not done:
                left = len(comps) - 1
                break
        return left, conns

    def wide_at(self, net, L, y, x):
        nid = self.fx.nid(net)
        w = self.fx.m.W[L][y, x]
        q = self.fx.m.Q[L][y, x]
        return (w == 0 or w == nid) and (q == 0 or q == nid)

    def geometry(self, net, st, path):
        narrow_w = WIDTHS.get(net, NARROW_POWER if st["special"] else TRACK)
        pts = []
        for i in path:
            li, y, x = split(i)
            L = LAYERS[li]
            w = WIDE if (st["wide"] and self.wide_at(net, L, y, x)) else narrow_w
            pts.append((L, y, x, w))
        tracks, vias = [], []

        def xy(p):
            return [coordx(p[2]), coordy(p[1])]

        def emit(a, b, w):
            if a[:3] != b[:3]:
                tracks.append({"net": net, "layer": a[0], "start": xy(a), "end": xy(b), "width": w})

        for end in (pts[0], pts[-1]):
            a = self.fx.anchor.get(end[:3])
            if a and (abs(a[0] - xy(end)[0]) > 1e-6 or abs(a[1] - xy(end)[1]) > 1e-6):
                tracks.append({"net": net, "layer": end[0], "start": list(a), "end": xy(end), "width": min(narrow_w, 0.3)})
        run_start = prev_p = pts[0]
        direction = run_w = None
        for p in pts[1:]:
            if p[0] != prev_p[0]:
                emit(run_start, prev_p, run_w or narrow_w)
                vias.append({"net": net, "at": xy(p), "dia": VIA_D, "drill": VIA_DRILL})
                run_start, direction, run_w = p, None, None
            else:
                d = (p[1] - prev_p[1], p[2] - prev_p[2])
                seg_w = min(p[3], prev_p[3])
                if direction is not None and (d != direction or seg_w != run_w):
                    emit(run_start, prev_p, run_w)
                    run_start = prev_p
                direction, run_w = d, seg_w
            prev_p = p
        emit(run_start, prev_p, run_w or narrow_w)
        return tracks, vias

    def stamp(self, net, sign=1):
        r = self.routes[net]
        if sign > 0 and not r["stamps"]:
            fx = self.fx
            special = net in POWER_CLASS or net in WIDE_SET
            for t in (t for c in r["conns"] for t in c["tracks"]):
                L = t["layer"][0]
                w = t["width"]
                if special:
                    w = NARROW_POWER if w <= NARROW_POWER else max(w, WIDE, VIA_D)
                for k, win, mask in fx.halo(L, "seg", (t["start"], t["end"], w), net):
                    r["stamps"].append((self.cnt[k][L], win, mask))
            for v in (v for c in r["conns"] for v in c["vias"]):
                p = v["at"]
                for L in LAYERS:
                    for k, win, mask in fx.halo(L, "via", (p, v["dia"]), net, extra=VIA_EXTRA):
                        r["stamps"].append((self.cnt[k][L], win, mask))
                win = fx.m.window((p[0], p[1], p[0], p[1]), VIA_PITCH)
                r["stamps"].append((self.vp, win, np.hypot(fx.m.X[win] - p[0], fx.m.Y[win] - p[1]) < VIA_PITCH - 1e-9))
        for arr, win, mask in r["stamps"]:
            arr[win][mask] += sign

    def conflicts(self, net, per_conn=False):
        r = self.routes[net]
        if not r["conns"]:
            return []
        st = self.prepare(net)
        self.stamp(net, -1)
        k_via = "W" if st["special"] else "V"
        out = []
        for c in r["conns"]:
            bad = []
            prev_li = None
            for i in c["cells"]:
                li, y, x = split(i)
                L = LAYERS[li]
                if st["wide"]:
                    k = "W" if self.wide_at(net, L, y, x) else "Q"
                else:
                    k = "Q" if st["special"] else "T"
                hit = self.cnt[k][L][y, x] > 0
                if not hit and prev_li is not None and prev_li != li:
                    hit = self.vp[y, x] > 0 or any(self.cnt[k_via][L2][y, x] > 0 for L2 in LAYERS)
                if hit:
                    bad.append((L, y, x))
                prev_li = li
            out.append(bad)
        self.stamp(net, 1)
        return out if per_conn else [q for bad in out for q in bad]

    def report_conflicts(self, order):
        seen = {}
        for n in order:
            for L, y, x in self.conflicts(n):
                px, py = coordx(x), coordy(y)
                for m in order:
                    if m == n:
                        continue
                    r = self.routes[m]
                    hit = any(t["layer"][0] == L and _seg_dist((px, py), t["start"], t["end"]) < 0.6 + t["width"] / 2 for c in r["conns"] for t in c["tracks"]) or \
                        any(math.hypot(v["at"][0] - px, v["at"][1] - py) < 0.9 for c in r["conns"] for v in c["vias"])
                    if hit:
                        seen.setdefault(tuple(sorted((n, m))), (L, px, py))
        for (a, b), (L, px, py) in sorted(seen.items(), key=lambda kv: (kv[1][1], kv[1][2])):
            print("  conflict %s x %s near %s (%.1f, %.1f)" % (a, b, L, px, py))

    def snapshot(self):
        return {n: {"conns": list(r["conns"]), "open": r["open"]} for n, r in self.routes.items()}

    def restore(self, snap):
        for n in list(self.routes):
            self.stamp(n, -1)
        self.routes = {n: dict(r, stamps=[]) for n, r in snap.items()}
        for n in self.routes:
            self.stamp(n)

    def renegotiate(self, order, seed, iterations=RENEG_ITERATIONS, pres=RENEG_PRES):
        todo = [seed]
        for _ in range(iterations):
            for n in todo:
                r = self.routes[n]
                per = self.conflicts(n, per_conn=True) if r["conns"] else []
                keep = [c for c, bad in zip(r["conns"], per) if not bad]
                self.stamp(n, -1)
                self.route(n, pres, keep=keep)
            conf = {n: self.conflicts(n) for n in order if n in self.routes}
            step = max(2, int(pres * HIST_FRACTION))
            for cells in conf.values():
                for L, y, x in cells:
                    self.hist[L][y, x] += step
            todo = [n for n in order if conf.get(n)]
            if not todo:
                return True
            pres = min(int(pres * PRES_GROWTH) + 1, PRES_MAX)
        return False

    def run(self):
        data = self.data

        def size(n):
            pts = [p["at"] for f in data["footprints"] for p in f["pads"] if p["net"] == n]
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            return (max(xs) - min(xs)) + (max(ys) - min(ys))

        first = [n for n in FIRST if n in self.nets]
        order = first + sorted((n for n in self.nets if n not in first), key=lambda n: (n in PLANES, size(n), n))
        seed = os.environ.get("AUTOROUTE_SEED")
        if seed:
            import random
            tail = order[len(first):]
            random.Random(int(seed)).shuffle(tail)
            order = order[:len(first)] + sorted(tail, key=lambda n: n in PLANES)
        print("autoroute: %d nets to route on a %d x %d grid" % (len(order), NX, NY), flush=True)
        t0 = time.time()
        resume = os.environ.get("AUTOROUTE_RESUME")
        if resume:
            self.routes = {}
            self.restore(pickle.load(open(resume, "rb")))
        else:
            pres = PRES_START
            todo = list(order)
            best = None
            for it in range(MAX_ITERATIONS):
                for n in todo:
                    if n in self.routes:
                        self.stamp(n, -1)
                    t1 = time.time()
                    self.route(n, pres)
                    if os.environ.get("AUTOROUTE_DEBUG"):
                        print("    %s: %.2f s, open %d" % (n, time.time() - t1, self.routes[n]["open"]), flush=True)
                conf = {n: self.conflicts(n) for n in order}
                step = max(2, int(pres * HIST_FRACTION))
                for cells in conf.values():
                    for L, y, x in cells:
                        self.hist[L][y, x] += step
                bad = [n for n in order if conf[n]]
                opened = [n for n in order if self.routes[n]["open"]]
                total = sum(len(c) for c in conf.values())
                print("iteration %d: %d nets rerouted, %d in conflict (%d cells), %d with open connections, %.0f s" % (
                    it + 1, len(todo), len(bad), total, len(opened), time.time() - t0), flush=True)
                score = (len(opened), total)
                if best is None or score < best[0]:
                    best = (score, it, self.snapshot())
                if not bad and not opened:
                    break
                if it - best[1] >= PATIENCE:
                    break
                todo = [n for n in order if conf[n] or self.routes[n]["open"]]
                pres = min(int(pres * PRES_GROWTH) + 1, PRES_MAX)
            self.restore(best[2])
            self.report_conflicts(order)
            ripped = set()
            while True:
                worst = None
                for n in order:
                    for k, bad in enumerate(self.conflicts(n, per_conn=True)):
                        key = (len(bad), -len(self.routes[n]["conns"][k]["cells"]))
                        if bad and (worst is None or key > worst[0]):
                            worst = (key, n, k)
                if worst is None:
                    break
                _, n, k = worst
                self.stamp(n, -1)
                r = self.routes[n]
                self.routes[n] = {"conns": r["conns"][:k] + r["conns"][k + 1:], "stamps": [], "open": r["open"] + 1}
                self.stamp(n)
                ripped.add(n)
            print("legalise: connections taken out of %d nets (%s)" % (len(ripped), ", ".join(sorted(ripped))), flush=True)
            for attempt in range(3):
                for n in sorted((n for n in order if self.routes[n]["open"]), key=lambda n: (size(n), n)):
                    self.stamp(n, -1)
                    self.route(n, 0, hard=True, keep=self.routes[n]["conns"])
                if not any(self.routes[n]["open"] for n in order):
                    break
            if self.state_file:
                pickle.dump(self.snapshot(), open(self.state_file, "wb"))

        def total_open():
            return sum(self.routes[n]["open"] for n in order)

        for round_no in range(RIPUP_ROUNDS):
            if total_open() == 0:
                break
            improved = False
            for n in sorted((n for n in order if self.routes[n]["open"]), key=lambda n: (size(n), n)):
                if not self.routes[n]["open"]:
                    continue
                snap, before = self.snapshot(), total_open()
                ok = self.renegotiate(order, n, pres=max(5, int(RENEG_PRES * RENEG_SCALE[round_no % len(RENEG_SCALE)])))
                if ok and total_open() < before:
                    improved = True
                    print("  %s: now %d open, %.0f s" % (n, total_open(), time.time() - t0), flush=True)
                else:
                    self.restore(snap)
            print("rip-up round %d: %d connections open, %.0f s" % (round_no + 1, total_open(), time.time() - t0), flush=True)
            if total_open() == 0 or (not improved and round_no >= len(RENEG_SCALE) - 1):
                break
        left = {n: self.routes[n]["open"] for n in order if self.routes[n]["open"]}
        tracks = [t for n in order for c in self.routes[n]["conns"] for t in c["tracks"]]
        vias = [v for n in order for c in self.routes[n]["conns"] for v in c["vias"]]
        return tracks, vias, left


def split_tees(tracks):
    """Split a segment where another segment of its net and layer ends on its interior: KiCad
    reports the end of a T junction in the middle of a segment as dangling."""
    ends = {}
    for t in tracks:
        for p in (t["start"], t["end"]):
            ends.setdefault((t["net"], t["layer"]), set()).add((round(p[0], 4), round(p[1], 4)))
    out = []
    for t in tracks:
        a, b = t["start"], t["end"]
        cuts = []
        for q in ends.get((t["net"], t["layer"]), ()):
            if q in ((round(a[0], 4), round(a[1], 4)), (round(b[0], 4), round(b[1], 4))):
                continue
            if _seg_dist(q, a, b) < 1e-4:
                vx, vy = b[0] - a[0], b[1] - a[1]
                cuts.append((((q[0] - a[0]) * vx + (q[1] - a[1]) * vy) / (vx * vx + vy * vy), q))
        pts = [a] + [list(q) for _, q in sorted(cuts)] + [b]
        for p0, p1 in zip(pts, pts[1:]):
            if p0 != p1:
                out.append(dict(t, start=p0, end=p1))
    return out


def main(src, dst):
    data = json.load(open(src))
    if os.environ.get("AUTOROUTE_FIXUP"):     # re-run only the output clean-up on an earlier result
        out = json.load(open(dst))
        out["tracks"] = split_tees(out["tracks"])
        json.dump(out, open(dst, "w"), indent=0)
        print("autoroute: %d tracks after splitting T junctions" % len(out["tracks"]))
        return
    r = Negotiator(data)
    r.state_file = dst + ".state"
    tracks, vias, left = r.run()
    tracks = split_tees(tracks)
    json.dump({"tracks": tracks, "vias": vias}, open(dst, "w"), indent=0)
    print("autoroute: %d tracks, %d vias, %d connections left open" % (len(tracks), len(vias), sum(left.values())))
    for n, k in sorted(left.items()):
        print("  open:", n, k)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
