"""Build the three.js replay page: `python -m viz.build [--check]`.

Writes viz/replay/data.json (kinematics, pressure, volume, crank torque and valve lift over one
720-deg cycle) and viz/replay/index.html, a self-contained page that embeds that data and the CAD
glTF exports (crankshaft, conrod, piston assembly, block, head; base64) and loads three.js and its
GLTFLoader / OrbitControls from jsDelivr at a pinned version. Sources:

* geometry: out/cad/*.glb from `python -m cad.build` (built first if missing), placed as cad/assembly.py
  does: crank rotated by -theta about +x (ADR 0002), conrods at the crankpin rotated by the exact
  slider-crank rod angle, pistons at the exact wrist-pin height (kin.slider_crank);
* p(theta): the converged M4 1D engine at the viz operating point, per cylinder (dyn.pressure);
* crank torque: dyn.forces (gas + reciprocating inertia) on the cylinder-mean M4 trace;
* valve lift: kin.cam profiles (spec/valvetrain.yaml).

`--check` then runs the headless Playwright check (viz/check_replay.mjs): no console errors, the
animation advances, screenshot at out/viz/replay.png. ALL NUMBERS ARE MODELLED.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

from redline.spec import ROOT, geometry, load

VIZ_DIR = ROOT / "viz"
CONFIG = VIZ_DIR / "viz.yaml"
TEMPLATE = VIZ_DIR / "template.html"
DEG = math.pi / 180.0
PA_PER_BAR = 1.0e5
M3_PER_CC = 1.0e-6
MM_PER_M = 1.0e3
CYCLE_DEG = 720.0
BYTES_PER_MB = 1.0e6


def config() -> dict:
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def cad_placement():
    """Cylinder axial positions [m], firing phases [deg], crank radius and rod length [m] from cad.params
    (read-only use of the CAD parameter set that produced the glb files)."""
    from cad.params import params
    p = params()
    n = p.cylinders
    return ([p.to_si_length(p.cyl_x(i)) for i in range(n)], [p.firing_phase_deg(i) for i in range(n)],
            p.to_si_length(p.crank_r), p.to_si_length(p.rod_len))


def kinematics(theta_deg: np.ndarray, phases_deg, r: float, l: float) -> dict:
    """Per-cylinder crankpin (y, z), conrod angle about +x and wrist-pin height, exactly as cad/assembly.py."""
    from kin.slider_crank import piston_position
    out = {"pin_y": [], "pin_z": [], "rod_angle": [], "piston_z": []}
    for ph in phases_deg:
        a = np.mod(theta_deg - ph, CYCLE_DEG / 2) * DEG
        py, pz = r * np.sin(a), r * np.cos(a)
        out["pin_y"].append(py)
        out["pin_z"].append(pz)
        out["rod_angle"].append(np.arcsin(py / l))
        out["piston_z"].append(piston_position(a, r, l))
    return out


def compute_data(cfg: dict | None = None) -> dict:
    from dyn import forces, params as dparams
    from dyn.pressure import M4_SOURCE, PressureTrace, m4_cylinder_traces
    from gasdyn.spec import valvetrain_spec
    from kin.cam import design_cams
    from kin.slider_crank import cylinder_volume

    cfg = cfg or config()
    spec = load()
    geom = geometry(spec)
    rpm = float(cfg["operating_point"]["rpm"])
    step = float(cfg["frames"]["step_deg"])
    theta = np.arange(0.0, CYCLE_DEG, step)
    xs, phases, r, l = cad_placement()

    grid, p_cyl, met = m4_cylinder_traces(rpm)
    mean_tr = PressureTrace.from_degrees(grid, p_cyl.mean(axis=0), source=M4_SOURCE)
    w = dparams.omega_from_rpm(rpm)
    T_cyl, T = forces.multi_cylinder_torque(theta * DEG, mean_tr, w)

    cams = design_cams(spec, valvetrain_spec())
    p_bar, vol_cc, lift_in, lift_ex = [], [], [], []
    for c, ph in enumerate(phases):
        loc = np.mod(theta - ph, CYCLE_DEG)
        p_bar.append(np.interp(loc, grid, p_cyl[c]) / PA_PER_BAR)
        vol_cc.append(cylinder_volume(loc * DEG, geom) / M3_PER_CC)
        lift_in.append(cams["intake"].valve_lift_crank(loc * DEG) * MM_PER_M)
        lift_ex.append(cams["exhaust"].valve_lift_crank(loc * DEG) * MM_PER_M)

    kin = kinematics(theta, phases, r, l)
    rnd = lambda a, d=6: [round(float(v), d) for v in np.asarray(a).ravel()]
    return {
        "modelled": True,
        "note": "MODELLED: Redline I4-2.0 NA digital twin, no physical engine.",
        "rpm": rpm,
        "throttle": cfg["operating_point"]["throttle"],
        "sources": {"pressure": M4_SOURCE, "torque": "dyn.forces on the cylinder-mean M4 trace (gas + inertia)",
                    "kinematics": "kin.slider_crank (exact), placement as cad/assembly.py (ADR 0001/0002)",
                    "valve_lift": "kin.cam (spec/valvetrain.yaml)", "cad": "out/cad/*.glb (python -m cad.build)"},
        "metrics": {"torque_nm": met["torque_nm"], "power_kw": met["power_kw"], "imep_net_bar": met["imep_net_bar"],
                    "p_max_bar": met["p_max_bar"], "bsfc_g_per_kwh": met["bsfc_g_per_kwh"],
                    "ve_trapped": met["ve_trapped"]},
        "geometry": {"crank_radius": r, "conrod_length": l, "cylinder_x": xs, "phase_deg": phases,
                     "bore": geom.bore, "firing_order": spec["engine"]["firing_order"]},
        "theta_deg": rnd(theta, 3),
        "crank_rotation_x": rnd(-theta * DEG),
        "pin_y": [rnd(a, 7) for a in kin["pin_y"]],
        "pin_z": [rnd(a, 7) for a in kin["pin_z"]],
        "rod_angle": [rnd(a, 7) for a in kin["rod_angle"]],
        "piston_z": [rnd(a, 7) for a in kin["piston_z"]],
        "p_bar": [rnd(a, 4) for a in p_bar],
        "volume_cc": [rnd(a, 4) for a in vol_cc],
        "torque_total_nm": rnd(T, 3),
        "torque_cyl_nm": [rnd(a, 3) for a in T_cyl],
        "torque_mean_nm": float(T.mean()),
        "lift_intake_mm": [rnd(a, 4) for a in lift_in],
        "lift_exhaust_mm": [rnd(a, 4) for a in lift_ex],
    }


def ensure_glb(cfg: dict) -> dict:
    d = ROOT / cfg["cad"]["glb_dir"]
    missing = [p for p in cfg["cad"]["parts"] if not (d / f"{p}.glb").exists()]
    if missing:
        subprocess.run([sys.executable, "-m", "cad.build"], cwd=ROOT, check=True)
    return {p: base64.b64encode((d / f"{p}.glb").read_bytes()).decode("ascii") for p in cfg["cad"]["parts"]}


def render_html(data: dict, glb: dict, cfg: dict) -> str:
    t = cfg["three"]
    base = f"{t['cdn']}/three@{t['version']}"
    page_cfg = {"playback": cfg["playback"], "cad": {k: cfg["cad"][k] for k in
                                                     ("parts", "translucent", "hidden_by_default", "translucent_opacity")}}
    html = TEMPLATE.read_text()
    html = html.replace("__THREE_MODULE__", f"{base}/build/three.module.js")
    html = html.replace("__THREE_ADDONS__", f"{base}/examples/jsm/")
    html = html.replace("__THREE_VERSION__", t["version"])
    html = html.replace("__CONFIG_JSON__", json.dumps(page_cfg))
    html = html.replace("__GLB_JSON__", json.dumps(glb))
    html = html.replace("__DATA_JSON__", json.dumps(data, separators=(",", ":")))
    return html


def build(cfg: dict | None = None) -> dict:
    cfg = cfg or config()
    out = ROOT / cfg["output"]["dir"]
    out.mkdir(parents=True, exist_ok=True)
    data = compute_data(cfg)
    (out / cfg["output"]["data"]).write_text(json.dumps(data, separators=(",", ":")))
    glb = ensure_glb(cfg)
    html = render_html(data, glb, cfg)
    (out / cfg["output"]["html"]).write_text(html)
    size_mb = sum((out / cfg["output"][k]).stat().st_size for k in ("html", "data")) / BYTES_PER_MB
    if size_mb > cfg["output"]["max_total_mb"]:
        raise RuntimeError(f"replay output {size_mb:.1f} MB exceeds {cfg['output']['max_total_mb']} MB")
    return {"html": str(out / cfg["output"]["html"]), "data": str(out / cfg["output"]["data"]), "size_mb": size_mb,
            "rpm": data["rpm"], "torque_nm": data["metrics"]["torque_nm"], "p_max_bar": data["metrics"]["p_max_bar"]}


def fetch_vendor(cfg: dict) -> Path:
    """Download the pinned three.js files the page imports into out/viz/vendor (curl, honours the proxy)."""
    t = cfg["three"]
    root = ROOT / "out" / "viz" / "vendor"
    for f in cfg["check"]["vendor_files"]:
        dst = root / f"three@{t['version']}" / f
        if dst.exists() and dst.stat().st_size > 0:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sSfL", "-o", str(dst), f"{t['cdn']}/three@{t['version']}/{f}"], check=True)
    return root


def check(cfg: dict | None = None) -> dict:
    """Headless Playwright check through node (viz/check_replay.mjs). Returns its JSON report."""
    cfg = cfg or config()
    fetch_vendor(cfg)
    node = shutil.which("node")
    if node is None:
        raise RuntimeError("node not found: the replay check needs node + playwright")
    html = ROOT / cfg["output"]["dir"] / cfg["output"]["html"]
    shot = ROOT / cfg["output"]["screenshot"]
    shot.parent.mkdir(parents=True, exist_ok=True)
    c = cfg["check"]
    cmd = [node, str(VIZ_DIR / "check_replay.mjs"), str(html), str(shot), str(c["settle_ms"]), str(c["sample_gap_ms"]),
           str(c["viewport"][0]), str(c["viewport"][1]), cfg["three"]["cdn"], str(ROOT / "out" / "viz" / "vendor")]
    res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    lines = [ln for ln in res.stdout.splitlines() if ln.startswith("{")]
    if not lines:
        raise RuntimeError(f"replay check produced no report:\n{res.stdout}\n{res.stderr}")
    rep = json.loads(lines[-1])
    rep["exit_code"] = res.returncode
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="run the headless Playwright check afterwards")
    a = ap.parse_args(argv)
    info = build()
    print(f"Redline replay (MODELLED): {info['rpm']:.0f} rpm WOT, M4 torque {info['torque_nm']:.1f} N m, "
          f"p_max {info['p_max_bar']:.1f} bar")
    print(f"  wrote {info['html']} and {info['data']} ({info['size_mb']:.2f} MB total)")
    if a.check:
        rep = check()
        print("  headless check: " + json.dumps(rep))
        if not rep.get("ok"):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
