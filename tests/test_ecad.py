"""ECU board (stretch goal, PLAN.md 3.7 / 6): schematic checks and layout checks.

Runs ecad/check_design.py (Zener build, firmware I/O map, ERC, datasheet limits) and
ecad/check_layout.py (outline, KiCad DRC, ground plane, decoupling, power copper, fab outputs).

Skips cleanly when the tools are missing:
  * the schematic checks need the `pcb` binary (built from source, docs/ecad.md);
  * the layout checks need the layout artefacts, which `ecad/layout/run.sh` makes with docker and
    the kicad/kicad:10.0.8 image; with docker available and RUN_ECAD_LAYOUT=1 this test runs the
    DRC step itself first (a minute), otherwise it checks the committed board as last dumped.

The board is a design, not hardware: nothing here is measured on a built ECU.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ECAD = Path(__file__).resolve().parent.parent / "ecad"
sys.path.insert(0, str(ECAD))


def _docker_ok():
    if not shutil.which("docker"):
        return False
    r = subprocess.run(["docker", "image", "inspect", "kicad/kicad:10.0.8"], capture_output=True)
    return r.returncode == 0


def test_schematic_checks():
    import check_design
    if not check_design.available():
        pytest.skip("pcb binary not found (see docs/ecad.md to build it)")
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "check_design.py"],
                       cwd=ECAD, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]


def test_layout_checks():
    board = ECAD / "boards/redline_ecu/layout/redline_ecu/redline_ecu.kicad_pcb"
    if not board.exists():
        pytest.skip("no KiCad board: run `ecad/layout/run.sh all` (docker)")
    if os.environ.get("RUN_ECAD_LAYOUT") == "1":
        if not _docker_ok():
            pytest.skip("docker or the kicad/kicad:10.0.8 image is not available")
        r = subprocess.run(["sh", "layout/run.sh", "dump", "drc"], cwd=ECAD, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    if not (ECAD / "layout/build/board.json").exists() or not (ECAD / "layout/build/drc.json").exists():
        if not _docker_ok():
            pytest.skip("no layout dump / DRC report and no docker to make them")
        r = subprocess.run(["sh", "layout/run.sh", "dump", "drc"], cwd=ECAD, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "check_layout.py"],
                       cwd=ECAD, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
