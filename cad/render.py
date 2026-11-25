"""Offscreen PBR renders of the CAD (three.js r160 in headless Chromium via Playwright, node).

Python writes one quantized glTF per part kind (cad.tess, smooth per-vertex normals from the exact surface) and a
scene description per view (instances with their 4x4 placements, PBR material keys, camera); cad/render.mjs serves
them to the page (cad/render.html) and screenshots each view at deviceScaleFactor 2 with a RoomEnvironment PMREM
light probe, ACES tone mapping and a shadow-casting key light. Geometry is the MODELLED Redline CAD.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np

from redline.spec import ROOT

CAD_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "out" / "viz" / "vendor"
THREE_VERSION = "0.160.0"
MM_PER_M = 1000.0
RENDER_LIN_MM = 0.05      # render meshes: absolute chordal deflection (mm), invisible at the render scale
RENDER_ANG_RAD = 0.12     # and angular deflection (rad); meshes are refined on the exact surface to both

# PBR materials (sRGB base colour, metalness, roughness, optional clearcoat). Families read apart under the
# neutral studio light: matte cast aluminium (block, head, pan, intake), bright machined steel (crank, rods, cams,
# valves), dark cast iron / black oxide (flywheel, liners, caps, bolts), copper-tin bearing shells, heat-tinted header.
MATERIALS = {
    "cast_al":     {"color": "#979ba0", "metalness": 0.35, "roughness": 0.66},
    "cast_al_dark": {"color": "#7b8087", "metalness": 0.35, "roughness": 0.62},
    "cast_al_warm": {"color": "#a09b92", "metalness": 0.35, "roughness": 0.64},
    "machined_al": {"color": "#c9ccd0", "metalness": 0.9, "roughness": 0.24},
    "forged_al":   {"color": "#bfc2c6", "metalness": 0.85, "roughness": 0.3},
    "steel":       {"color": "#9ca3ab", "metalness": 1.0, "roughness": 0.2},
    "steel_dark":  {"color": "#5f656c", "metalness": 0.95, "roughness": 0.32},
    "iron":        {"color": "#4a4d51", "metalness": 0.6, "roughness": 0.58},
    "bearing":     {"color": "#b8875a", "metalness": 1.0, "roughness": 0.28},
    "bolt":        {"color": "#26282c", "metalness": 0.6, "roughness": 0.42},
    "spring":      {"color": "#2f5d8a", "metalness": 0.4, "roughness": 0.35, "clearcoat": 0.6},
    "valve":       {"color": "#b4b9bf", "metalness": 1.0, "roughness": 0.18},
    "gasket":      {"color": "#55595e", "metalness": 0.6, "roughness": 0.45},
    "rubber":      {"color": "#1c1d20", "metalness": 0.0, "roughness": 0.8},
    "black_paint": {"color": "#2a2c30", "metalness": 0.2, "roughness": 0.5, "clearcoat": 0.3},
    "header":      {"color": "#5c524b", "metalness": 0.75, "roughness": 0.5},
    "red_paint":   {"color": "#9e2a25", "metalness": 0.3, "roughness": 0.4, "clearcoat": 0.5},
}


def loc_matrix(loc) -> list:
    """build123d Location -> column-major 4x4 (metres) for three.js Matrix4.fromArray."""
    t = loc.wrapped.Transformation()
    m = np.eye(4)
    for r in range(3):
        for c in range(4):
            m[r, c] = t.Value(r + 1, c + 1)
    m[:3, 3] /= MM_PER_M
    return m.T.ravel().tolist()


def write_part_glbs(shapes: dict, work: Path, lin=RENDER_LIN_MM, ang=RENDER_ANG_RAD) -> dict:
    """{key: shape} -> {key: glb filename}; quantized, smooth normals. Cached by key in `work`."""
    from cad.tess import triangulate, write_glb
    work.mkdir(parents=True, exist_ok=True)
    out = {}
    for key, shape in shapes.items():
        fn = f"{key}.glb"
        if not (work / fn).exists():   # cache: delete `work` after a geometry change
            write_glb({key: triangulate(shape, lin, ang)}, work / fn, quantize=True)
        out[key] = fn
    return out


def render(views: list, work: Path, out_dir: Path, width: int = 1600, height: int = 1000, scale: int = 2,
           timeout_s: int = 900) -> dict:
    """Render `views` ({name, objects: [{glb, material, matrix}], camera: {pos, target, fov}, ...}); returns
    {view name: png path}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    _vendor()
    shutil.copy(CAD_DIR / "render.html", work / "index.html")
    scene = {"materials": MATERIALS, "views": views, "width": width, "height": height}
    (work / "scene.json").write_text(json.dumps(scene))
    node = shutil.which("node")
    if node is None:
        raise RuntimeError("node not found (renders need node + playwright)")
    res = subprocess.run([node, str(CAD_DIR / "render.mjs"), str(work), str(out_dir), str(width), str(height), str(scale),
                          str(VENDOR)], capture_output=True, text=True, timeout=timeout_s, cwd=ROOT)
    lines = [ln for ln in res.stdout.splitlines() if ln.startswith("{")]
    if res.returncode != 0 or not lines:
        raise RuntimeError(f"render failed:\n{res.stdout[-3000:]}\n{res.stderr[-3000:]}")
    rep = json.loads(lines[-1])
    if rep.get("errors"):
        raise RuntimeError(f"render page errors: {rep['errors']}")
    return {v["name"]: str(out_dir / f"{v['name']}.png") for v in views}


def _vendor():
    files = ["build/three.module.js", "examples/jsm/environments/RoomEnvironment.js",
             "examples/jsm/loaders/GLTFLoader.js", "examples/jsm/utils/BufferGeometryUtils.js"]
    files += [f"examples/jsm/{f}" for f in (   # ambient occlusion: GTAO post-processing chain
        "postprocessing/EffectComposer.js", "postprocessing/RenderPass.js", "postprocessing/OutputPass.js",
        "postprocessing/GTAOPass.js", "postprocessing/Pass.js", "postprocessing/ShaderPass.js",
        "postprocessing/MaskPass.js", "shaders/GTAOShader.js", "shaders/PoissonDenoiseShader.js",
        "shaders/CopyShader.js", "shaders/OutputShader.js", "math/SimplexNoise.js")]
    for f in files:
        dst = VENDOR / f"three@{THREE_VERSION}" / f
        if dst.exists() and dst.stat().st_size:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sSfL", "-o", str(dst), f"https://cdn.jsdelivr.net/npm/three@{THREE_VERSION}/{f}"],
                       check=True)


def look(pos_dir, target, dist, fov=28.0):
    """Camera at `dist` mm from target along the unit direction pos_dir (CAD frame, z up)."""
    d = np.asarray(pos_dir, float)
    d /= np.linalg.norm(d)
    return {"pos": (np.asarray(target) + dist * d).tolist(), "target": list(target), "fov": fov}


# --------------------------------------------------------------------------------------------- audit views
def material_for(kind: str) -> str:
    k = kind
    table = [("exhaust_stay_bolt", "bolt"), ("exhaust_stay", "steel_dark"), ("exhaust_header", "header"), ("bracket_bolt", "bolt"), ("intake_bracket", "steel_dark"),
             ("throttle", "machined_al"),
             ("block", "cast_al"), ("head_gasket", "gasket"), ("head_bolt", "bolt"), ("head", "cast_al_warm"),
             ("liner", "iron"), ("main_cap_bolt", "bolt"), ("main_cap", "iron"), ("shell", "bearing"),
             ("crank", "steel"), ("conrod_bolt", "bolt"), ("conrod", "steel_dark"), ("piston", "forged_al"),
             ("wrist_pin", "steel"), ("ring", "steel_dark"), ("clip", "steel"), ("oil_pan", "cast_al_dark"),
             ("pan_bolt", "bolt"), ("flywheel_bolt", "bolt"), ("ring_gear", "steel_dark"), ("flywheel", "iron"),
             ("damper", "steel_dark"), ("sprocket", "steel"), ("pulley", "steel"), ("belt", "rubber"),
             ("tensioner", "steel"), ("idler", "steel"), ("cam_cap_bolt", "bolt"), ("cam_cap", "cast_al_warm"),
             ("cam", "steel"), ("bucket", "steel"), ("shim", "steel"),
             ("spring_seat", "steel_dark"), ("spring", "spring"), ("retainer", "steel"), ("collet", "steel_dark"), ("guide", "bearing"),
             ("seat_insert", "iron"), ("valve", "valve"), ("spark_plug", "steel"), ("manifold_bolt", "bolt"),
             ("intake", "cast_al"), ("throttle", "machined_al"), ("exhaust", "header"), ("header", "header")]
    for key, mat in table:
        if key in k:
            return mat
    return "steel"


def assembly_views(theta: float = 0.0, adv: float = 0.0, hide=(), p=None):
    """Objects for the whole assembly at a crank angle, sharing one glb per part kind."""
    from cad import assembly as asm
    from cad.params import params
    p = p or params()
    insts = [i for i in asm.instances(p) if not any(h in i.kind for h in hide)]
    mot = asm.motions(p)
    objs = []
    for it in insts:
        loc = mot[it.group].loc(theta, adv) * it.base
        objs.append({"glb": f"{it.kind}.glb", "material": material_for(it.kind), "matrix": loc_matrix(loc),
                     "name": it.name})
    return objs, {i.kind for i in insts}


def audit_views(tag: str, out_dir: Path) -> dict:
    """Quick evidence renders of the audited assembly (3/4 front, bottom, side)."""
    from cad import assembly as asm
    from cad.params import params
    p = params()
    work = ROOT / "out" / "cad" / "render_work" / f"audit_{tag}"
    src = asm.parts(p)
    objs, kinds = assembly_views(0.0, 0.0, p=p)
    write_part_glbs({k: src[k][0] for k in kinds}, work)
    ctr = [0.0, 0.0, 150.0]
    views = [
        {"name": f"audit_{tag}_front34", "objects": objs, "camera": look([-1.0, -1.15, 0.75], ctr, 1250)},
        {"name": f"audit_{tag}_bottom", "objects": objs, "camera": look([-0.6, -0.9, -0.8], [0, 0, 40], 1050)},
        {"name": f"audit_{tag}_side", "objects": objs, "camera": look([0.0, -1.0, 0.05], ctr, 1300)},
    ]
    return render(views, work, out_dir, width=1200, height=800, scale=1)


# --------------------------------------------------------------------------------------------- final renders
EXPLODE_MM = {   # exploded view: offset (x, y, z) per part family, mm (first match wins)
    "intake_manifold": (0, -140, 260), "throttle_body": (90, -140, 260), "intake_bracket": (120, 0, 120),
    "bracket_bolt_block": (170, 0, 120), "bracket_bolt_pad": (120, 0, 80),
    "manifold_bolt_intake": (0, -110, 0), "exhaust_header": (0, 170, 0), "exhaust_stay_bolt_block": (0, 135, 0),
    "exhaust_stay_bolt_pad": (0, 170, 40), "exhaust_stay": (0, 120, 0),
    "manifold_bolt_exhaust": (0, 110, 0), "cam_cap_bolt": (0, 0, 330), "cam_cap": (0, 0, 290),
    "camshaft": (0, 0, 250), "cam_pulley": (-60, 0, 250), "timing_belt": (-110, 0, 60), "idler": (-90, 0, 0),
    "tensioner": (-90, 0, 0), "bucket": (0, 0, 215), "shim": (0, 0, 215), "retainer": (0, 0, 200),
    "collet": (0, 0, 200), "spring_seat": (0, 0, 175), "spring": (0, 0, 185), "head_bolt": (0, 0, 300),
    "spark_plug": (0, 0, 260), "valve": (0, 0, 120), "guide": (0, 0, 160), "seat_insert": (0, 0, 160),
    "head_gasket": (0, 0, 80), "head": (0, 0, 160), "piston": (0, 0, 60), "wrist_pin": (0, 0, 60),
    "ring": (0, 0, 60), "clip": (0, 0, 60), "conrod_bolt": (0, 0, -150), "conrod_cap": (0, 0, -130),
    "big_end_shell_lower": (0, 0, -115), "conrod": (0, 0, 20), "big_end_shell_upper": (0, 0, 20),
    "main_cap_bolt": (0, 0, -260), "main_cap": (0, 0, -230), "main_shell_lower": (0, 0, -205),
    "crankshaft": (0, 0, -170), "crank_sprocket": (-60, 0, -170), "damper": (-110, 0, -170),
    "flywheel_bolt": (110, 0, -170), "flywheel": (80, 0, -170), "ring_gear": (80, 0, -170),
    "main_shell_upper": (0, 0, -150), "pan_bolt": (0, 0, -330), "oil_pan": (0, 0, -300),
}


def _explode(kind: str):
    for k, v in EXPLODE_MM.items():
        if kind.startswith(k):
            return v
    return (0, 0, 0)


def section_objects(theta: float, x_cut: float, keep_side: float, p, work: Path, adv: float = 0.0, hide=(),
                    tag: str = "sec"):
    """Exact CAD section: every placed instance crossing the plane x = x_cut is cut by the half space (build123d
    boolean), so the cut faces are real faces of the solids (water jacket, ports, bearing shells visible)."""
    from build123d import Box, Pos
    from cad import assembly as asm
    from cad.audit import bbox
    from cad.tess import triangulate, write_glb
    src, mot = asm.parts(p), asm.motions(p)
    big = 4000.0
    half = Pos(x_cut + keep_side * big / 2, 0, 0) * Box(big, big, big)
    objs, shapes = [], {}
    work.mkdir(parents=True, exist_ok=True)
    for it in asm.instances(p):
        if any(it.kind.startswith(h) for h in hide):
            continue
        loc = mot[it.group].loc(theta, adv) * it.base
        shp = src[it.kind][0].moved(loc)
        b = bbox(shp.wrapped)
        lo, hi = (b[0], b[3])
        if (keep_side > 0 and hi <= x_cut) or (keep_side < 0 and lo >= x_cut):
            continue
        if lo < x_cut < hi:
            cut = shp & half
            if cut.volume <= 0:
                continue
            fn = f"{tag}_{it.name}.glb"
            write_glb({it.name: triangulate(cut, RENDER_LIN_MM, RENDER_ANG_RAD)}, work / fn, quantize=True)
            objs.append({"glb": fn, "material": material_for(it.kind), "matrix": loc_matrix(asm.IDENTITY),
                         "name": it.name})
        else:
            objs.append({"glb": f"{it.kind}.glb", "material": material_for(it.kind), "matrix": loc_matrix(loc),
                         "name": it.name})
            shapes[it.kind] = src[it.kind][0]
    return objs, shapes


def clearance_annotations(p, theta: float, pairs) -> list:
    """Measured gap lines for a view: exact minimum distance (BRepExtrema) between the placed instances of each
    (kind_a, kind_b, label) in `pairs` at crank angle `theta`; returns [{a, b, label, mm}] in CAD mm."""
    from cad import assembly as asm
    from cad.audit import distance
    src, mot = asm.parts(p), asm.motions(p)
    insts = asm.instances(p)
    out = []
    for ka, kb, label, extra in pairs:
        best = None
        for A in (i for i in insts if i.kind == ka):
            sa = src[A.kind][0].wrapped.Moved((mot[A.group].loc(theta) * A.base).wrapped)
            for B in (i for i in insts if i.kind == kb):
                sb = src[B.kind][0].wrapped.Moved((mot[B.group].loc(theta) * B.base).wrapped)
                d, qa, qb = distance(sa, sb)
                if best is None or d < best[0]:
                    best = (d, qa, qb)
        out.append({"a": best[1], "b": best[2], "label": f"{label} {best[0]:.1f} mm", "mm": best[0], **extra})
    return out


ASSEMBLY_FRAME = {"width": 1400, "height": 1100}   # landscape: the header runs rearward past the flywheel


def final_renders(out_dir: Path | None = None, only=None) -> dict:
    """The documentation renders (1600 x 1000 at deviceScaleFactor 2, GTAO ambient occlusion, auto-framed): full
    assembly 3/4 front and 3/4 rear (exhaust side), section through cylinder 1, valve-train cutaway (head sectioned
    through cylinder 1's rear valve pair), intake cutaway (runner 2 plane: bellmouth, curl, port), bottom end with
    the pan removed, exhaust-to-pan/flywheel clearance view with the measured gaps, exploded view. `only`: a set
    of view names to render (the rest are skipped)."""
    from build123d import Pos
    from cad import assembly as asm
    from cad.params import params
    p = params()
    out_dir = out_dir or ROOT / "out" / "cad" / "renders"
    work = ROOT / "out" / "cad" / "render_work" / "final"
    src = asm.parts(p)
    # crank angle with cylinder 1 on its intake stroke near full intake lift (valves visibly open)
    theta_open = float(np.degrees(asm.vt.cam_profiles(p)["intake"].centre_crank))
    views, kinds = [], set()

    objs, k = assembly_views(0.0, 0.0, p=p)
    kinds |= k
    views.append({"name": "assembly_front34", "objects": objs, **ASSEMBLY_FRAME,
                  "camera": {**look([-1.0, -1.25, 0.75], [0, 0, 150], 2500, 26), "fit": 0.95}, "key": [-0.6, -0.9, 1.5]})
    views.append({"name": "assembly_rear34_exhaust", "objects": objs, **ASSEMBLY_FRAME,
                  "camera": {**look([1.0, 1.3, 0.55], [0, 0, 150], 2500, 26), "fit": 0.95}, "key": [0.5, 0.9, 1.5]})

    hide_sec = ("intake_manifold", "exhaust_header", "exhaust_stay", "manifold_bolt", "throttle_body", "intake_bracket",
                "bracket_bolt")
    sec, shapes = section_objects(theta_open, p.cyl_x(0), 1.0, p, work, hide=hide_sec, tag="sec")
    kinds |= set(shapes)
    views.append({"name": "section_cyl1", "objects": sec, "width": 1250, "height": 1300,
                  "camera": {**look([-1.0, -0.32, 0.18], [p.cyl_x(0), 0, 170], 1050, 30), "fit": 0.9},
                  "key": [-1.0, -0.5, 1.2]})

    # valve-train cutaway: everything cut at the axis plane of cylinder 1's rear valve pair, so the cam towers,
    # bucket bores, springs, guides, seats and ports under the caps and cams are visible in section
    x_vt = p.cyl_x(0) + max(xl for _k, _s, xl, _y, _d in p.valve_layout())
    hide_vt = hide_sec + ("timing_belt", "cam_pulley", "idler", "tensioner", "oil_pan", "pan_bolt", "crank", "conrod",
                          "main_", "big_end", "flywheel", "ring_gear", "damper", "piston", "wrist_pin", "ring_",
                          "clip", "liner")
    vt_objs, shapes = section_objects(theta_open, x_vt, 1.0, p, work, hide=hide_vt, tag="vt")
    kinds |= set(shapes)
    views.append({"name": "valvetrain_cutaway", "objects": vt_objs,
                  "camera": {**look([-1.0, -0.62, 0.5], [x_vt + 40, 0, 290], 560, 30)}, "key": [-0.7, -0.5, 1.6],
                  "ground": False})

    objs, k = assembly_views(90.0, 0.0, hide=("oil_pan", "pan_bolt") + hide_sec, p=p)
    kinds |= k
    views.append({"name": "bottom_end_pan_removed", "objects": objs,
                  "camera": {**look([-0.35, -0.55, -1.0], [0, 0, 10], 1150, 30), "fit": 0.9},
                  "key": [-0.3, -0.5, -1.4], "ground": False})

    objs, k = assembly_views(0.0, 0.0, hide=("intake_manifold", "throttle_body", "intake_bracket", "bracket_bolt"), p=p)
    kinds |= k
    # measured header gaps (the same exact distances as docs/CAD-AUDIT.md key clearances)
    # (the flywheel disc, 77 mm, sits just inboard of its ring gear, so the ring gear gap is the flywheel's minimum)
    ann = clearance_annotations(p, 0.0, [
        ("exhaust_header", "oil_pan", "header to oil pan", {"offset": [60, -40, -60], "size": 0.032}),
        ("exhaust_header", "pan_bolt", "header to pan bolt", {"offset": [-90, 40, -60], "size": 0.032}),
        ("exhaust_header", "ring_gear", "header to flywheel ring gear", {"offset": [80, 40, 40], "size": 0.032})])
    views.append({"name": "exhaust_clearance", "objects": objs, "annotations": ann, **ASSEMBLY_FRAME,
                  "camera": {**look([0.5, 0.75, -0.6], [150, 120, -60], 1600, 28), "fit": 0.95},
                  "key": [0.3, 0.8, 1.4], "ground": False})

    # intake cutaway: everything cut on runner 2's centre-line plane (the swan neck lies in that plane), so the bore
    # from the head port through the curl to the bellmouth lip inside the plenum is visible in section
    from cad import manifolds as mf
    x_in = p.cyl_x(1)
    g = mf._runner_y(p)
    pg = mf.plenum_geom(p)
    hide_in = ("exhaust_header", "exhaust_stay", "manifold_bolt_exhaust", "timing_belt", "cam_pulley", "idler", "tensioner", "oil_pan",
               "pan_bolt", "crank", "conrod", "main_", "big_end", "flywheel", "ring_gear", "damper", "piston",
               "wrist_pin", "ring_", "clip", "liner", "block", "head_gasket", "head_bolt", "intake_bracket",
               "bracket_bolt")
    in_objs, shapes = section_objects(0.0, x_in, 1.0, p, work, hide=hide_in, tag="in")
    kinds |= set(shapes)
    tgt = [x_in, 0.5 * (g["yc"] + pg["y1"]), 0.5 * (mf.port_flange_z(p, "intake") + pg["z1"])]
    views.append({"name": "intake_cutaway", "objects": in_objs,
                  "camera": {**look([-1.0, -0.55, 0.35], tgt, 900, 30), "fit": 0.92}, "key": [-0.8, -0.6, 1.5],
                  "ground": False})

    objs = []
    for it in asm.instances(p):
        loc = Pos(*_explode(it.kind)) * asm.motions(p)[it.group].loc(0.0) * it.base
        objs.append({"glb": f"{it.kind}.glb", "material": material_for(it.kind), "matrix": loc_matrix(loc),
                     "name": it.name})
        kinds.add(it.kind)
    views.append({"name": "exploded", "objects": objs, **ASSEMBLY_FRAME,
                  "camera": {**look([-1.0, -1.35, 0.55], [0, 0, 150], 2600, 26), "fit": 0.93}, "key": [-0.6, -0.9, 1.5]})

    if only:
        views = [v for v in views if v["name"] in only]
    write_part_glbs({k: src[k][0] for k in kinds}, work)
    return render(views, work, out_dir, width=1600, height=1000, scale=2, timeout_s=3600)


def save_doc_png(src: Path, dst: Path, colors: int = 256):
    """Committed copy of a render: 256-colour palette PNG (dithered), ~5x smaller than the RGB original."""
    from PIL import Image
    Image.open(src).convert("RGB").quantize(colors=colors, method=Image.Quantize.MEDIANCUT,
                                            dither=Image.Dither.FLOYDSTEINBERG).save(dst, optimize=True)


if __name__ == "__main__":   # python -m cad.render [view ...] -> out/cad/renders/*.png and docs/renders/*.png
    import sys
    docs = ROOT / "docs" / "renders"
    docs.mkdir(parents=True, exist_ok=True)
    for name, path in final_renders(only=set(sys.argv[1:]) or None).items():
        save_doc_png(Path(path), docs / Path(path).name)
        print(name, path)
