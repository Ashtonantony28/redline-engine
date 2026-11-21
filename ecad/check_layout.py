#!/usr/bin/env python3
"""Layout checks for redline_ecu: outline, KiCad DRC, ground plane, decoupling, power copper, fab outputs.

Verify:  cd ecad && layout/run.sh all fab && python -m pytest -q check_layout.py

The checks read the geometry dump written by layout/dump.py (mm, origin at the board centre,
x right, y down, seen from the top) and KiCad's own DRC report (layout/build/drc.json), so pytest
itself does not need KiCad.  Limits are stated here; never loosen one to make a layout pass:
fix the layout or report it.  Adapted from drone-design ecad/check_layout.py.
"""
import json
import math
import os
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "layout"))
import board as B  # noqa: E402  (layout/board.py: outline, holes, net classes)

PCB = HERE / "boards/redline_ecu/layout/redline_ecu/redline_ecu.kicad_pcb"
DUMP = HERE / "layout/build/board.json"
DRC = HERE / "layout/build/drc.json"
FAB = HERE / "fab"

TOL = 0.011
HOLE_DRILL = 3.2            # M3 clearance
MCU_DECAP_MAX = 2.0         # ours: 100 nF / VCAP capacitor pad to its MCU supply pin
MCU_DECAP_1U_MAX = 3.0      # ours: the 1 uF VDDA capacitor (second tier behind the 100 nF)
IC_DECAP_MAX = 4.0          # ours: other ICs, capacitor pad to the IC's supply pin
PLANE_FILL_MIN = 0.80       # In1.Cu ground plane: one piece covering >= 80 % of the board
WIDE_MIN = 0.3              # narrowest track of a high-current net (at its pins)
WIDE_SHARE_MIN = 0.5        # share of a high-current net's track length at the full 1.0 mm
ACT_NETS = ("ETB_M_P", "ETB_M_N", "INJ1", "INJ2", "INJ3", "INJ4", "VVT", "LSU_HTR_N")
FAB_FILES = ("gerbers/redline_ecu-F_Cu.gbr", "gerbers/redline_ecu-B_Cu.gbr", "gerbers/redline_ecu-In1_Cu.gbr",
             "gerbers/redline_ecu-In2_Cu.gbr", "gerbers/redline_ecu-Edge_Cuts.gbr", "gerbers/redline_ecu-PTH.drl", "redline_ecu-pos.csv", "redline_ecu.step")


def pad_gap(a, b):
    dx = max(a[0] - b[2], b[0] - a[2], 0)
    dy = max(a[1] - b[3], b[1] - a[3], 0)
    return math.hypot(dx, dy)


@pytest.fixture(scope="module")
def lay():
    if not DUMP.exists():
        pytest.skip("no layout dump: run `layout/run.sh dump` (needs docker + kicad/kicad:10.0.8)")
    return json.loads(DUMP.read_text())


def test_outline_and_mounting_holes(lay):
    x0, y0, x1, y1 = lay["edge"]   # bounding box of the 0.1 mm outline lines
    assert abs((x1 - x0 - 0.1) - B.W) < TOL and abs((y1 - y0 - 0.1) - B.H) < TOL, "board is not %g x %g mm" % (B.W, B.H)
    holes = [(f["at"], p) for f in lay["footprints"] if f["ref"].startswith("H") for p in f["pads"] if p["drill"] > 0]
    assert len(holes) == 4
    for (x, y), p in holes:
        assert abs(abs(x) - B.HOLE_XY[0]) < TOL and abs(abs(y) - B.HOLE_XY[1]) < TOL
        assert abs(p["drill"] - HOLE_DRILL) < TOL
    keepouts = [z for z in lay["zones"] if z["keepout"] and z["name"].startswith("HOLE")]
    assert len(keepouts) == 4 and all(len(z["layers"]) == 4 for z in keepouts)


def test_connectors_on_the_edge(lay):
    fps = {f["path"]: f for f in lay["footprints"]}
    for path in ("CONN.J1.X", "CONN.J2.X"):
        assert abs(fps[path]["courtyard"][3] - B.H / 2) < 0.05, "%s front is not on the bottom edge" % path
    assert abs(fps["COM.J3.X"]["courtyard"][1] + B.H / 2) < 0.05, "USB-C front is not on the top edge"


def test_all_on_top(lay):
    assert all(f["side"] == "F" for f in lay["footprints"]), "single-sided assembly: every part on top"


def test_drc_clean():
    if not DRC.exists():
        pytest.skip("no DRC report: run `layout/run.sh drc`")
    r = json.loads(DRC.read_text())
    errors = [v for v in r.get("violations", []) if v["severity"] == "error"]
    assert not errors, "%d DRC errors, first: %s" % (len(errors), errors[0]["description"])
    assert not r.get("unconnected_items"), "%d unconnected items" % len(r["unconnected_items"])


def test_drc_report_is_current():
    if not DRC.exists() or not PCB.exists():
        pytest.skip("no DRC report")
    assert DRC.stat().st_mtime >= PCB.stat().st_mtime - 1, "drc.json is older than the board: rerun `layout/run.sh drc`"


def test_ground_plane_is_solid(lay):
    planes = [z for z in lay["zones"] if z["name"] == "L2_GND"]
    assert len(planes) == 1 and planes[0]["layers"] == ["In1.Cu"] and planes[0]["net"] == "GND"
    filled = planes[0]["filled"]["In1.Cu"]
    assert filled["outlines"] == 1, "In1.Cu ground plane is in %d pieces" % filled["outlines"]
    assert filled["area"] >= PLANE_FILL_MIN * B.W * B.H, "ground plane covers %.0f mm2" % filled["area"]
    # nothing but ground on L2: no tracks there
    assert not [t for t in lay["tracks"] if t["layer"] == "In1.Cu"]


def test_mcu_decoupling(lay):
    fps = {f["path"]: f for f in lay["footprints"]}
    sys.path.insert(0, str(HERE / "layout"))
    from autoplace import MCU_DECAPS
    mcu = {p["number"]: p for p in fps["MCU.U1.X"]["pads"]}
    worst = []
    for path, pin in MCU_DECAPS.items():
        if "BULK" in path:
            continue    # the 10 uF bulk capacitor: anywhere on the 3V3 net near the MCU
        cap = fps[path]
        pin_pad = mcu[pin]
        gap = min(pad_gap(p["box"], pin_pad["box"]) for p in cap["pads"] if p["net"] == pin_pad["net"])
        worst.append((gap - (MCU_DECAP_1U_MAX if "1U" in path else MCU_DECAP_MAX), gap, path))
    over, gap, path = max(worst)
    assert over <= 0, "%s is %.2f mm from its MCU pin" % (path, gap)


def test_ic_decoupling(lay):
    fps = {f["path"]: f for f in lay["footprints"]}
    worst = []
    for path, f in fps.items():
        m = re.match(r"(\w+)\.C_(U\d+)_(VDD|VCC|VIO|IN|OUT)\w*\.C$", path)
        if not m:
            continue
        ic = fps.get("%s.%s.X" % (m.group(1), m.group(2)))
        if not ic:
            continue
        nets = {p["net"] for p in f["pads"]} - {"GND"}
        gaps = [pad_gap(p["box"], q["box"]) for p in f["pads"] for q in ic["pads"] if p["net"] in nets and q["net"] == p["net"]]
        if gaps:
            worst.append((min(gaps), path))
    assert len(worst) >= 15
    gap, path = max(worst)
    assert gap <= IC_DECAP_MAX, "%s is %.2f mm from its IC's supply pin" % (path, gap)


def test_high_current_copper(lay):
    """High-current nets: no track under WIDE_MIN, most of the length at the full width, and a pour."""
    pours = {z["net"] for z in lay["zones"] if z["name"].startswith("PWR_") and any(v["area"] > 0 for v in z["filled"].values())}
    for net in B.WIDE_NETS:
        ts = [t for t in lay["tracks"] if t["net"] == net]
        if not ts:
            continue
        assert min(t["width"] for t in ts) >= WIDE_MIN - TOL, "%s has a track under %.1f mm" % (net, WIDE_MIN)
        total = sum(t["len"] for t in ts)
        wide = sum(t["len"] for t in ts if t["width"] >= B.WIDE - TOL)
        if net in ACT_NETS:
            assert wide >= WIDE_SHARE_MIN * total, "%s: %.0f %% of its track at %.1f mm" % (net, 100 * wide / total, B.WIDE)
        assert net in pours, "%s has no filled pour" % net


def test_fab_outputs():
    if not FAB.exists():
        pytest.skip("no fab outputs: run `layout/run.sh fab`")
    missing = [f for f in FAB_FILES if not (FAB / f).exists()]
    assert not missing, "missing fab outputs: %s" % ", ".join(missing)
    pos = (FAB / "redline_ecu-pos.csv").read_text().splitlines()
    assert len(pos) - 1 >= 200, "pick-and-place file lists %d parts" % (len(pos) - 1)


def summary(lay):
    r = json.loads(DRC.read_text()) if DRC.exists() else {}
    return {
        "footprints": len(lay["footprints"]), "tracks": len(lay["tracks"]), "vias": len(lay["vias"]),
        "track_length_mm": round(sum(t["len"] for t in lay["tracks"]), 1),
        "drc_errors": sum(1 for v in r.get("violations", []) if v["severity"] == "error"),
        "drc_warnings": sum(1 for v in r.get("violations", []) if v["severity"] == "warning"),
        "unconnected": len(r.get("unconnected_items", [])),
    }


if __name__ == "__main__":
    print(json.dumps(summary(json.loads(DUMP.read_text())), indent=1))
    sys.exit(pytest.main(["-q", __file__]))
