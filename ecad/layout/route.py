"""Routing steps for redline_ecu that need KiCad's Python (run through layout/run.sh).

  route.py apply  <board> <tracks.json> locked|unlocked   add tracks and vias (preroute / autoroute)
  route.py finish <board>                                 pours, dangling vias out, silk, refill, save

Coordinates in the JSON are mm from the board centre (layout/board.py).
Adapted from drone-design ecad/layout/route.py.
"""
import json
import math
import os
import sys

import pcbnew

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from board import EDGE_INSET, ORIGIN, PLANE_NETS, WIDE_NETS  # noqa: E402
from place import add_zone, board_rect, pt  # noqa: E402

MM = pcbnew.FromMM
LAYER = {"F": pcbnew.F_Cu, "B": pcbnew.B_Cu, "I": pcbnew.In2_Cu, "F.Cu": pcbnew.F_Cu, "B.Cu": pcbnew.B_Cu,
         "In1.Cu": pcbnew.In1_Cu, "In2.Cu": pcbnew.In2_Cu}
POUR_HALF = 1.6       # a high-current pour reaches this far either side of its track centre
POUR_CLEARANCE = 0.25
GND_CLEARANCE = 0.25


def add_track(board, t, locked):
    tr = pcbnew.PCB_TRACK(board)
    tr.SetStart(pt(*t["start"]))
    tr.SetEnd(pt(*t["end"]))
    tr.SetWidth(MM(t["width"]))
    tr.SetLayer(LAYER[t["layer"]])
    tr.SetNetCode(board.GetNetcodeFromNetname(t["net"]))
    tr.SetLocked(locked)
    board.Add(tr)


def add_via(board, v, locked):
    via = pcbnew.PCB_VIA(board)
    via.SetPosition(pt(*v["at"]))
    via.SetViaType(pcbnew.VIATYPE_THROUGH)
    via.SetLayerPair(pcbnew.F_Cu, pcbnew.B_Cu)
    via.SetWidth(pcbnew.F_Cu, MM(v["dia"]))
    via.SetDrill(MM(v["drill"]))
    via.SetNetCode(board.GetNetcodeFromNetname(v["net"]))
    via.SetLocked(locked)
    board.Add(via)


def apply(board_path, json_path, locked="locked"):
    locked = locked == "locked"
    board = pcbnew.LoadBoard(board_path)
    if locked:
        assert not list(board.GetTracks()), "board already has tracks: run place first"
    data = json.load(open(json_path))
    for t in data["tracks"]:
        if t["start"] != t["end"]:
            add_track(board, t, locked)
    for v in data["vias"]:
        add_via(board, v, locked)
    pcbnew.ZONE_FILLER(board).Fill(board.Zones())
    board.Save(board_path)
    print("applied %d tracks, %d vias" % (len(data["tracks"]), len(data["vias"])))


def prune_vias(board, tracks, footprints):
    """Vias that ended up with copper on one layer only (not counting the ground plane): a drill hit."""
    segs = {}
    for t in tracks:
        if t.Type() != pcbnew.PCB_VIA_T:
            a, b = t.GetStart(), t.GetEnd()
            segs.setdefault(t.GetNetCode(), []).append((a.x, a.y, b.x, b.y, t.GetLayer()))
    pads = {}
    for fp in footprints:
        for p in fp.Pads():
            pads.setdefault(p.GetNetCode(), []).append(p)

    def touches(c, r, seg):
        ax, ay, bx, by, _ = seg
        vx, vy = bx - ax, by - ay
        L2 = vx * vx + vy * vy
        t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((c.x - ax) * vx + (c.y - ay) * vy) / L2))
        return (ax + t * vx - c.x) ** 2 + (ay + t * vy - c.y) ** 2 <= r * r

    gone = []
    for v in tracks:
        if v.Type() != pcbnew.PCB_VIA_T or v.GetNetname() in PLANE_NETS:
            continue
        c, r = v.GetPosition(), v.GetWidth(pcbnew.F_Cu) / 2
        layers = {seg[4] for seg in segs.get(v.GetNetCode(), []) if touches(c, r, seg)}
        for p in pads.get(v.GetNetCode(), []):
            if p.GetBoundingBox().Contains(c):
                layers.update(l for l in (pcbnew.F_Cu, pcbnew.B_Cu) if p.IsOnLayer(l))
        if len(layers) < 2:
            gone.append(v)
    names = ["%s @(%.2f, %.2f)" % (v.GetNetname(), pcbnew.ToMM(v.GetPosition().x) - ORIGIN[0], pcbnew.ToMM(v.GetPosition().y) - ORIGIN[1]) for v in gone]
    return gone, names


def capsule(a, b, r, n=6):
    """Polygon round a segment: a rectangle with half-octagon ends (mm, board-centre)."""
    ang = math.atan2(b[1] - a[1], b[0] - a[0])
    pts = []
    for k in range(n + 1):      # end cap at b
        t = ang - math.pi / 2 + math.pi * k / n
        pts.append((b[0] + r * math.cos(t), b[1] + r * math.sin(t)))
    for k in range(n + 1):      # end cap at a
        t = ang + math.pi / 2 + math.pi * k / n
        pts.append((a[0] + r * math.cos(t), a[1] + r * math.sin(t)))
    return pts


def finish(board_path):
    import silk
    board = pcbnew.LoadBoard(board_path)
    tracks, footprints, zones = list(board.GetTracks()), list(board.GetFootprints()), list(board.Zones())
    names = {z.GetZoneName() for z in zones}
    outline = board_rect(EDGE_INSET)
    # The routed layers connect the header ground pins solid, so that a track passing a pin cannot
    # starve its thermal relief; the In1 plane keeps thermal spokes on them (place.py) for soldering.
    for name, layer, tht in (("L1_GND", pcbnew.F_Cu, False), ("L3_GND", pcbnew.In2_Cu, False), ("L4_GND", pcbnew.B_Cu, False)):
        if name not in names:
            add_zone(board, "GND", layer, outline, name=name, clearance=GND_CLEARANCE, thermal_tht=tht)
    # High-current nets: a pour along each of their tracks, on the track's layer, above the
    # ground pour.  The fill keeps clearance to everything else, so it only adds copper where
    # there is room; each net has its own priority so that two such pours never overlap.
    pours = 0
    if not any(n.startswith("PWR_") for n in names):
        merged = {}
        for t in tracks:
            if t.Type() == pcbnew.PCB_VIA_T or t.GetNetname() not in WIDE_NETS:
                continue
            a = (pcbnew.ToMM(t.GetStart().x) - ORIGIN[0], pcbnew.ToMM(t.GetStart().y) - ORIGIN[1])
            b = (pcbnew.ToMM(t.GetEnd().x) - ORIGIN[0], pcbnew.ToMM(t.GetEnd().y) - ORIGIN[1])
            half = max(POUR_HALF, pcbnew.ToMM(t.GetWidth()) / 2)
            piece = pcbnew.SHAPE_POLY_SET()
            piece.NewOutline()
            for x, y in capsule(a, b, half):
                piece.Append(MM(ORIGIN[0] + x), MM(ORIGIN[1] + y))
            key = (t.GetNetname(), t.GetLayer())
            if key in merged:
                merged[key].BooleanAdd(piece)
            else:
                merged[key] = piece
        # one zone per connected piece of the union: overlapping zones of one net are a DRC error
        for (net, layer), poly in sorted(merged.items()):
            poly.Simplify()
            for i in range(poly.OutlineCount()):
                z = add_zone(board, net, layer, [], priority=10 + WIDE_NETS.index(net), name="PWR_%s" % net,
                             clearance=POUR_CLEARANCE, thermal_tht=False, min_width=0.25)
                o = z.Outline()
                o.RemoveAllContours()
                o.AddOutline(poly.Outline(i))
                for h in range(poly.HoleCount(i)):
                    o.AddHole(poly.Hole(i, h))
                pours += 1
    gone, pruned = prune_vias(board, tracks, footprints)
    silk.tidy(board)
    for v in gone:
        board.Remove(v)
    pcbnew.ZONE_FILLER(board).Fill(board.Zones())
    board.Save(board_path)
    print("finished: ground pours on F.Cu, In2.Cu, B.Cu; %d high-current pours" % pours)
    print("finished: %d dangling vias removed%s" % (len(pruned), (": " + ", ".join(pruned)) if pruned else ""))


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"apply": apply, "finish": finish}.get(cmd, lambda *a: sys.exit(__doc__))(*args)
