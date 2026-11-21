#!/usr/bin/env python3
"""Review pictures of the layout from layout/dump.py output (host Python + Pillow).

    python3 layout/plot.py <dump.json> <out_prefix> [drc.json]   -> <out_prefix>_top.png, _inner.png

_top: courtyards, pads, F.Cu (red) and B.Cu (blue) tracks, vias, designators.
_inner: the In2.Cu signal layer (orange) over the pads.  DRC markers (from drc.json) are
drawn as magenta rings on both.
"""
import json
import sys

from PIL import Image, ImageDraw, ImageFont

PX = 10  # pixels per mm
COL = {"F.Cu": (200, 50, 50), "B.Cu": (50, 80, 200), "In2.Cu": (220, 150, 0), "In1.Cu": (50, 160, 50)}


def draw(d, path, layers, drc=None):
    x0, y0, x1, y1 = d["edge"]
    pad = 2
    img = Image.new("RGB", (int((x1 - x0 + 2 * pad) * PX), int((y1 - y0 + 2 * pad) * PX)), "white")
    g = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    def P(x, y):
        return ((x - x0 + pad) * PX, (y - y0 + pad) * PX)

    for e in d["edges"]:
        g.line([P(*e["start"]), P(*e["end"])], fill="black", width=2)
    g.rectangle([P(x0, y0), P(x1, y1)], outline=(160, 160, 160))
    for z in d["zones"]:
        if z["keepout"]:
            pts = [P(*q) for q in z["outline"]]
            g.line(pts + [pts[0]], fill=(150, 150, 150), width=1)
    for f in d["footprints"]:
        if f["courtyard"]:
            c = f["courtyard"]
            g.rectangle([P(c[0], c[1]), P(c[2], c[3])], outline=(190, 190, 190))
        for p in f["pads"]:
            b = p["box"]
            col = (120, 120, 120) if p["net"] != "GND" else (90, 150, 90)
            if p["drill"] > 0:
                g.ellipse([P(b[0], b[1]), P(b[2], b[3])], outline=col, width=2)
            else:
                g.rectangle([P(b[0], b[1]), P(b[2], b[3])], fill=col)
    for t in d["tracks"]:
        if t["layer"] in layers:
            g.line([P(*t["start"]), P(*t["end"])], fill=COL[t["layer"]], width=max(1, int(t["width"] * PX)))
    for v in d["vias"]:
        r = v["dia"] / 2
        x, y = v["at"]
        g.ellipse([P(x - r, y - r), P(x + r, y + r)], fill=(40, 40, 40))
    for f in d["footprints"]:
        g.text(P(f["at"][0], f["at"][1]), f["ref"], fill=(0, 0, 0), font=font, anchor="mm")
    if drc:
        cx, cy = 100.0 + (x0 + x1) / 2, 100.0 + (y0 + y1) / 2  # dump is centred, DRC is in sheet mm
        for v in drc.get("violations", []) + drc.get("unconnected_items", []):
            if v.get("severity", "error") != "error" and "severity" in v:
                continue
            for i in v["items"]:
                x, y = i["pos"]["x"] - cx, i["pos"]["y"] - cy
                g.ellipse([P(x - 1, y - 1), P(x + 1, y + 1)], outline=(255, 0, 255), width=3)
    img.save(path)


def main(src, prefix, drc_path=None):
    d = json.load(open(src))
    drc = json.load(open(drc_path)) if drc_path else None
    draw(d, prefix + "_top.png", ("F.Cu", "B.Cu"), drc)
    draw(d, prefix + "_inner.png", ("In2.Cu",), drc)
    print("wrote %s_top.png, %s_inner.png" % (prefix, prefix))


if __name__ == "__main__":
    main(*sys.argv[1:])
