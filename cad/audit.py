"""CAD clash / clearance / support / tessellation audit: `python -m cad.audit [--tag final] [--jobs 4]`.

What it checks (all on the exact BRep solids, CAD units mm, MODELLED geometry):

(a) static clash: every part instance against every other at the static pose (crank 0 deg, VVT parked). Pairs
    whose bounding boxes (inflated by the contact tolerance) do not overlap cannot touch and are skipped; every
    other pair gets an exact boolean common volume (BRepAlgoAPI_Common) and an exact minimum distance
    (BRepExtrema_DistShapeShape). A pair is an OVERLAP if the common volume exceeds tolerances.overlap_mm3, a
    CONTACT if the distance is below tolerances.contact_mm, else a clearance.
(b) moving sweep: every pair of instances in different motion groups, at least one moving, whose swept bounding
    boxes come within the near distance. The minimum distance over the cycle is found by branch and bound on the
    crank angle with a rigorous Lipschitz bound: the distance between two rigid bodies changes by at most
    L = v_A + v_B per degree, v = largest surface speed of the body (mm/deg, evaluated on its bounding-box
    corners, which bound every surface point because |v(p)| is convex in p). An interval [a, b] with sampled
    distances d_a, d_b cannot contain a value below (d_a + d_b - L (b - a)) / 2, so intervals are bisected until
    that bound is within the stated resolution of the best sample (or the interval is below min_step_deg); the
    certified lower bound is reported with the minimum. Initial grid: grid_step_deg (2 deg). Valve groups moving
    against static parts are swept in valve lift (1 mm/mm) instead of crank angle. Pairs allowed to stay in
    contact while moving (bearings, sliding fits) are sampled on the grid and checked for zero penetration.
(c) support graph: parts joined by contacts or overlaps; any part not connected to the block is FLOATING.
(d) tessellation: every source solid meshed with the export settings (spec/cad.yaml export) and measured.
(e) renders: cad.render views of the audited model.

Intended contacts and overlaps come from spec/cad_fits.yaml (each with its build reason); any other overlap or
contact is a DEFECT. Outputs out/cad/audit_<tag>.json and .md.

The full run takes ~2.5 h on 4 cores. It is resumable: finished pairs are appended to out/cad/audit_<tag>.ckpt/
as they complete (progress.json there shows the stage and counts), and rerunning the same command after a kill
continues where it stopped as long as the source hash is unchanged (`--fresh` discards the checkpoint).
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import heapq
import json
import math
import re
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import yaml
from OCP.Bnd import Bnd_Box
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

from cad import assembly as asm
from cad.params import params
from redline.spec import ROOT

FITS_PATH = ROOT / "spec" / "cad_fits.yaml"
OUT_DIR = ROOT / "out" / "cad"
SRC_FILES = ("cad/*.py", "spec/cad.yaml", "spec/cad_fits.yaml", "spec/engine.yaml", "spec/valvetrain.yaml",
             "spec/gasdyn.yaml", "spec/dynamics.yaml", "spec/fea.yaml")
H_DEG = 0.05            # central-difference step for surface speeds (deg)
VEL_SAMPLES_DEG = 1.0   # crank-angle spacing of the surface-speed scan
VEL_SAFETY = 1.05       # margin on the sampled speed maximum


# ------------------------------------------------------------------------------------------------ geometry
def bbox(shape) -> np.ndarray:
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    lo, hi = b.CornerMin(), b.CornerMax()
    return np.array([lo.X(), lo.Y(), lo.Z(), hi.X(), hi.Y(), hi.Z()])


def boxes_gap(a: np.ndarray, b: np.ndarray) -> float:
    """Separation of two AABBs (0 if they overlap)."""
    d = np.maximum(0.0, np.maximum(a[:3] - b[3:], b[:3] - a[3:]))
    return float(np.linalg.norm(d))


def distance(a, b):
    d = BRepExtrema_DistShapeShape(a, b)
    d.Perform()
    if not d.IsDone() or d.NbSolution() == 0:
        return math.inf, None, None
    pa, pb = d.PointOnShape1(1), d.PointOnShape2(1)
    return d.Value(), [pa.X(), pa.Y(), pa.Z()], [pb.X(), pb.Y(), pb.Z()]


def common_volume(a, b) -> float:
    c = BRepAlgoAPI_Common(a, b)
    c.Build()
    if not c.IsDone():
        return math.nan
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(c.Shape(), g)
    return abs(g.Mass())


HASH_EXCLUDE = ("render.py",)   # cameras/materials of the documentation renders: no effect on audited geometry


def source_hash() -> str:
    h = hashlib.sha256()
    for pat in SRC_FILES:
        for f in sorted(ROOT.glob(pat)):
            if f.name in HASH_EXCLUDE:
                continue
            h.update(f.name.encode())
            h.update(f.read_bytes())
    return h.hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ fits
class Fits:
    def __init__(self, path: Path = FITS_PATH):
        d = yaml.safe_load(path.read_text())
        self.tol = {k: float(v) for k, v in d["tolerances"].items()}
        self.entries = d.get("fits", []) or []
        self.sweep = d.get("sweep", {})
        self.report = d.get("report_clearances", []) or []

    @staticmethod
    def _same(rule, a, b) -> bool:
        if rule is None:
            return True
        ka, kb = _key(a, rule), _key(b, rule)
        return ka is not None and ka == kb

    def match(self, a, b):
        """Fit entry covering instances a, b (order free), or None."""
        for e in self.entries:
            pa, pb = e["parts"]
            for x, y in ((a, b), (b, a)):
                if fnmatch.fnmatchcase(x.kind, pa) and fnmatch.fnmatchcase(y.kind, pb) and self._same(e.get("same"), x, y):
                    return e
        return None


def _key(inst, rule: str):
    """Grouping key for `same:` rules: 'cyl', or 'site' (the instance name with its kind prefix removed)."""
    if rule == "site":
        m = SITE_RE.search(inst.name)
        return m.group(1) if m else inst.name
    return getattr(inst, rule, None)


# site of an instance: the valve (cyl_kind_n), the cam side and journal (kind_j), the main bearing / bolt index, or
# the cam side alone (camshaft_intake, cam_pulley_intake)
SITE_RE = re.compile(r"_(\d+_(?:intake|exhaust)_\d+|(?:intake|exhaust)_\d+|\d+|intake|exhaust)[ab]?$")


# ------------------------------------------------------------------------------------------------ worker state
_S: dict = {}


def _init_state(tag_theta: float, adv: float):
    p = params()
    insts = asm.instances(p)
    src = asm.parts(p)
    mot = asm.motions(p)
    shapes = [src[i.kind][0].wrapped.Moved((mot[i.group].loc(tag_theta, adv) * i.base).wrapped) for i in insts]
    # springs deform: at the static pose each spring is rebuilt at the length its valve's lift gives
    from cad import valvetrain as vt
    for k, it in enumerate(insts):
        if _is_spring(it):
            m = SITE_RE.search(it.name).group(1)              # "<cyl>_<kind>_<n>"
            g = mot[f"valve_{m}"]
            kind = m.split("_")[1]
            lift = g.lift(tag_theta, adv)
            spring = asm.Pos(0, 0, vt.stack(p, kind)["spring_bottom"]) * vt.spring(
                p, round(vt.spring_located_length(p, kind, lift), 6))
            shapes[k] = spring.wrapped.Moved(it.base.wrapped)
    vvt_max = float(p.engine["valves"]["vvt_intake_range_deg"])
    _S.update(p=p, insts=insts, src=src, mot=mot, shapes=shapes, theta=tag_theta, adv=adv, vvt_max=vvt_max)


def _static_pair(ij):
    i, j = ij
    a, b = _S["shapes"][i], _S["shapes"][j]
    t0 = time.time()
    vol = common_volume(a, b)
    d, pa, pb = distance(a, b)
    return i, j, vol, d, pa, pb, time.time() - t0


# ------------------------------------------------------------------------------------------------ kinematics
def group_speed(p, mot, group: str, corners: np.ndarray, period: float, adv: float = 0.0) -> float:
    """Largest surface speed (mm/deg) of the bounding-box corners (local frame) of a body in `group`."""
    m = mot[group]
    if m.static:
        return 0.0
    vmax = 0.0
    for th in np.arange(0.0, period, VEL_SAMPLES_DEG):
        l0, l1 = m.loc(th - H_DEG, adv), m.loc(th + H_DEG, adv)
        for c in corners:
            a = _apply(l0, c)
            b = _apply(l1, c)
            vmax = max(vmax, float(np.linalg.norm(b - a)) / (2 * H_DEG))
    return vmax * VEL_SAFETY


def _apply(loc, c) -> np.ndarray:
    t = loc.wrapped.Transformation()
    from OCP.gp import gp_Pnt
    q = gp_Pnt(*c).Transformed(t)
    return np.array([q.X(), q.Y(), q.Z()])


def _corners(b: np.ndarray) -> np.ndarray:
    return np.array([[b[0 + 3 * i], b[1 + 3 * j], b[2 + 3 * k]] for i in (0, 1) for j in (0, 1) for k in (0, 1)])


# ------------------------------------------------------------------------------------------------ sweep
def _near_faces(shape, box: np.ndarray, pad: float):
    """Compound of the faces of `shape` whose bounding boxes come within `pad` of `box` (a static part only needs
    the faces a moving part can reach; penetration cannot start without touching one of them, and the static
    clash check rules out starting inside)."""
    from OCP.BRep import BRep_Builder
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS_Compound
    comp = TopoDS_Compound()
    bld = BRep_Builder()
    bld.MakeCompound(comp)
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    n = 0
    while exp.More():
        f = exp.Current()
        if boxes_gap(bbox(f), box) <= pad:
            bld.Add(comp, f)
            n += 1
        exp.Next()
    return comp if n else None


def _box_solid(b: np.ndarray):
    """Solid axis-aligned box from a bbox array (xmin, ymin, zmin, xmax, ymax, zmax)."""
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
    from OCP.gp import gp_Pnt as P
    return BRepPrimAPI_MakeBox(P(*b[:3]), P(*b[3:])).Solid()


def _sweep_pair(job):
    """Minimum distance of one moving pair over the period (see module docstring (b))."""
    i, j, period, L, contact_mode, adv, cfg, sbox_i, sbox_j, param = job
    insts, src, mot = _S["insts"], _S["src"], _S["mot"]
    A, B = insts[i], insts[j]
    sa, sb = src[A.kind][0].wrapped, src[B.kind][0].wrapped
    ma, mb = mot[A.group], mot[B.group]
    if param == "lift":   # parameter = valve lift (mm) instead of crank angle
        ma, mb = (_LiftMotion(ma) if ma.lift_loc else ma), (_LiftMotion(mb) if mb.lift_loc else mb)
        cfg = dict(cfg, grid_step_deg=cfg["lift_grid_mm"], min_step_deg=cfg["lift_min_step_mm"])
    # a static partner only needs its faces near the mover's swept box
    if ma.static:
        sa = _near_faces(sa.Moved(A.base.wrapped), sbox_j, cfg["near_mm"]) or sa.Moved(A.base.wrapped)
        a_base = asm.IDENTITY
    else:
        a_base = A.base
    if mb.static:
        sb = _near_faces(sb.Moved(B.base.wrapped), sbox_i, cfg["near_mm"]) or sb.Moved(B.base.wrapped)
        b_base = asm.IDENTITY
    else:
        b_base = B.base
    evals = 0
    t0 = time.time()

    def shapes(th):
        return (sa.Moved((ma.loc(th, adv) * a_base).wrapped), sb.Moved((mb.loc(th, adv) * b_base).wrapped))

    def full_shapes(th):   # whole solids (a static partner above may be reduced to its nearby faces)
        fa, fb = src[A.kind][0].wrapped, src[B.kind][0].wrapped
        return (fa.Moved((ma.loc(th, adv) * A.base).wrapped) if not ma.static else fa.Moved(A.base.wrapped),
                fb.Moved((mb.loc(th, adv) * B.base).wrapped) if not mb.static else fb.Moved(B.base.wrapped))

    def at(th):
        nonlocal evals
        evals += 1
        return distance(*shapes(th))

    if ma.static != mb.static and not contact_mode:
        # certified screen: exact distance from the static part to the mover's swept box (a solid box that holds the
        # mover at every parameter value) is a lower bound of the distance over the whole period
        st, box = (sa, sbox_j) if ma.static else (sb, sbox_i)
        lb = distance(st, _box_solid(box))[0]
        if lb > cfg["contact_mm"]:
            d0, pa, pb = at(0.0)
            return {"a": A.name, "b": B.name, "period_deg": period, "lipschitz_mm_per_deg": L, "mode": contact_mode,
                    "vvt_advance_deg": adv, "param": param, "min_mm": d0, "at_deg": 0.0, "point_a": pa, "point_b": pb,
                    "certified_lb_mm": lb, "screened": True, "evals": evals, "seconds": time.time() - t0}
    grid = [float(t) for t in np.arange(0.0, period + 1e-9, cfg["grid_step_deg"])]
    samples = {th: at(th) for th in grid}
    tol_c = cfg["contact_mm"]
    out = {"a": A.name, "b": B.name, "period_deg": period, "lipschitz_mm_per_deg": L, "mode": contact_mode,
           "vvt_advance_deg": adv, "param": param}
    if contact_mode:
        worst_v, worst_th = 0.0, None
        dth = min(samples, key=lambda t: samples[t][0])
        vstep = cfg["contact_volume_step_mm"] if param == "lift" else cfg["contact_volume_step_deg"]
        for th, (d, _pa, _pb) in samples.items():
            on_grid = abs(th / vstep - round(th / vstep)) < 1e-9 or th == dth
            if d <= tol_c and on_grid:
                v = common_volume(*full_shapes(th))
                if v > worst_v:
                    worst_v, worst_th = v, th
        out.update(min_mm=samples[dth][0], at_deg=dth, max_overlap_mm3=worst_v, overlap_at_deg=worst_th,
                   certified_lb_mm=None, evals=evals, seconds=time.time() - t0)
        return out
    # local refinement of the best grid sample (bisection on the bracketing intervals)
    best_th = min(samples, key=lambda t: samples[t][0])
    h = cfg["grid_step_deg"]
    while h > cfg["min_step_deg"]:
        h /= 2
        for th in (best_th - h, best_th + h):
            if 0.0 <= th <= period and th not in samples:
                samples[th] = at(th)
        best_th = min(samples, key=lambda t: samples[t][0])
    best = samples[best_th][0]
    # certification: no angle can be closer than `target`
    target = best - max(cfg["resolution_mm"], cfg["certify_fraction"] * best)
    keys = sorted(samples)
    heap = [(0.5 * (samples[a][0] + samples[b][0] - L * (b - a)), a, b) for a, b in zip(keys[:-1], keys[1:])]
    heapq.heapify(heap)
    unresolved = []
    while heap and evals < cfg["max_evals"]:
        lb, a, b = heapq.heappop(heap)
        if lb >= target:
            heapq.heappush(heap, (lb, a, b))
            break
        if b - a <= cfg["min_step_deg"]:
            unresolved.append(lb)
            continue
        m = 0.5 * (a + b)
        samples[m] = at(m)
        if samples[m][0] < best:
            best, best_th = samples[m][0], m
        for x, y in ((a, m), (m, b)):
            heapq.heappush(heap, (0.5 * (samples[x][0] + samples[y][0] - L * (y - x)), x, y))
    cert = min([hh[0] for hh in heap[:1]] + unresolved + [best]) if heap or unresolved else best
    if heap:
        cert = min(cert, min(hh[0] for hh in heap))
    d, pa, pb = samples[best_th]
    vol = 0.0
    if best <= tol_c:
        vol = common_volume(*full_shapes(best_th))
    out.update(min_mm=best, at_deg=best_th, point_a=pa, point_b=pb, certified_lb_mm=max(cert, 0.0),
               overlap_at_min_mm3=vol, evals=evals, seconds=time.time() - t0)
    return out


class _LiftMotion:
    """A valve group's motion parametrized by its lift u (mm) instead of the crank angle."""
    static = False

    def __init__(self, m):
        self.m = m

    def loc(self, u, adv=0.0):
        return self.m.lift_loc(u)


# ------------------------------------------------------------------------------------------------ main audit
class Checkpoint:
    """Resumable run state in out/cad/audit_<tag>.ckpt/: every finished static pair and sweep job is appended to a
    JSON-lines file as it completes (with the source hash of the run in meta.json), the tessellation report is
    stored once done, and progress.json carries the stage and counts. A restarted run with the same source hash
    skips everything already on disk; a different hash (or --fresh) starts over."""

    def __init__(self, out_dir: Path, tag: str, src_hash: str, fresh: bool = False):
        self.dir = out_dir / f"audit_{tag}.ckpt"
        meta = self.dir / "meta.json"
        if self.dir.exists() and (fresh or not meta.exists()
                                  or json.loads(meta.read_text()).get("source_hash") != src_hash):
            for f in self.dir.iterdir():
                f.unlink()
        self.dir.mkdir(parents=True, exist_ok=True)
        if not meta.exists():
            meta.write_text(json.dumps({"source_hash": src_hash, "tag": tag, "started": time.time()}))
        self.progress = {}

    def load(self, name: str) -> dict:
        """{key: record} of a JSON-lines stage file (a torn last line from a kill is ignored)."""
        out = {}
        f = self.dir / f"{name}.jsonl"
        if f.exists():
            for ln in f.read_text().splitlines():
                try:
                    r = json.loads(ln)
                except json.JSONDecodeError:
                    continue
                out[r["key"]] = r["rec"]
        return out

    def appender(self, name: str):
        return open(self.dir / f"{name}.jsonl", "a", buffering=1)

    @staticmethod
    def add(fh, key: str, rec):
        fh.write(json.dumps({"key": key, "rec": rec}, default=_jsonable) + "\n")
        fh.flush()
        os.fsync(fh.fileno())

    def get(self, name: str):
        f = self.dir / f"{name}.json"
        return json.loads(f.read_text()) if f.exists() else None

    def put(self, name: str, obj):
        tmp = self.dir / f"{name}.json.tmp"
        tmp.write_text(json.dumps(obj, default=_jsonable))
        tmp.replace(self.dir / f"{name}.json")

    def stage(self, **kw):
        self.progress.update(kw, updated=time.strftime("%Y-%m-%d %H:%M:%S"))
        self.put("progress", self.progress)


def run(tag: str = "final", jobs: int = 4, renders: bool = True, sweep: bool = True, out_dir: Path = OUT_DIR,
        fresh: bool = False) -> dict:
    t_start = time.time()
    src_hash = source_hash()   # hash of the sources this run is built from (taken before any geometry is made)
    ck = Checkpoint(out_dir, tag, src_hash, fresh)
    fits = Fits()
    tol = fits.tol
    _init_state(0.0, 0.0)
    p, insts, shapes, mot, src = _S["p"], _S["insts"], _S["shapes"], _S["mot"], _S["src"]
    n = len(insts)
    print(f"[audit:{tag}] {n} part instances")
    boxes = [bbox(s) for s in shapes]

    # ---- (a) static clash ------------------------------------------------------------------------------
    cand = [(i, j) for i in range(n) for j in range(i + 1, n) if boxes_gap(boxes[i], boxes[j]) <= tol["contact_mm"]]
    n_pairs = n * (n - 1) // 2
    print(f"[audit:{tag}] static: {n_pairs} pairs, {len(cand)} with touching bounding boxes", flush=True)
    done = ck.load("static")
    todo = [ij for ij in cand if f"{insts[ij[0]].name}|{insts[ij[1]].name}" not in done]
    print(f"[audit:{tag}] static: {len(cand) - len(todo)} pairs from the checkpoint, {len(todo)} to run", flush=True)
    ck.stage(stage="static", static_total=len(cand), static_done=len(cand) - len(todo))
    ctx = mp.get_context("fork")
    with ck.appender("static") as fh, ctx.Pool(jobs) as pool:
        for k, r in enumerate(pool.imap_unordered(_static_pair, todo, chunksize=4)):
            key = f"{insts[r[0]].name}|{insts[r[1]].name}"
            Checkpoint.add(fh, key, list(r))
            done[key] = list(r)
            if (k + 1) % 50 == 0 or k + 1 == len(todo):
                ck.stage(static_done=len(cand) - len(todo) + k + 1)
    res = [done[f"{insts[i].name}|{insts[j].name}"] for i, j in cand]
    print(f"[audit:{tag}] static done in {time.time() - t_start:.0f} s", flush=True)
    static = []
    for i, j, vol, d, pa, pb, dt in res:
        A, B = insts[i], insts[j]
        fit = fits.match(A, B)
        if vol > tol["overlap_mm3"]:
            state = "overlap"
        elif d <= tol["contact_mm"]:
            state = "contact"
        else:
            state = "clear"
        allowed = None
        if state == "overlap":
            allowed = bool(fit and fit["type"] == "overlap" and vol <= float(fit.get("max_overlap_mm3", math.inf)))
        elif state == "contact":
            allowed = bool(fit)
        static.append({"a": A.name, "b": B.name, "kind_a": A.kind, "kind_b": B.kind, "state": state,
                       "overlap_mm3": vol, "distance_mm": d, "point": pa, "allowed": allowed,
                       "reason": fit["reason"] if fit else None, "seconds": dt})
    # required contacts that are missing
    missing = []
    touching = {(s["a"], s["b"]) for s in static if s["state"] in ("contact", "overlap")}
    touching |= {(b, a) for a, b in touching}
    for e in fits.entries:
        if not e.get("required"):
            continue
        for A in insts:
            if not fnmatch.fnmatchcase(A.kind, e["parts"][0]):
                continue
            partners = [B for B in insts if B is not A and fnmatch.fnmatchcase(B.kind, e["parts"][1])
                        and Fits._same(e.get("same"), A, B)]
            if partners and not any((A.name, B.name) in touching for B in partners):
                missing.append({"part": A.name, "expected_contact_with": e["parts"][1], "reason": e["reason"]})

    # ---- (c) support graph -------------------------------------------------------------------------------
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    idx = {it.name: k for k, it in enumerate(insts)}
    for s in static:
        if s["state"] in ("contact", "overlap"):
            parent[find(idx[s["a"]])] = find(idx[s["b"]])
    root = find(idx["block"])
    floating = [it.name for k, it in enumerate(insts) if find(k) != root]
    # parts held only through defects (unlisted contacts): rebuild with allowed links only
    parent = list(range(n))
    for s in static:
        if s["state"] in ("contact", "overlap") and s["allowed"]:
            parent[find(idx[s["a"]])] = find(idx[s["b"]])
    root = find(idx["block"])
    floating_allowed_only = [it.name for k, it in enumerate(insts) if find(k) != root]

    # ---- (b) moving sweep ---------------------------------------------------------------------------------
    sweeps = []
    if sweep:
        sweeps = moving_sweep(fits, boxes, jobs, tag, ck)

    ck.stage(stage="valve_checks")
    valve = valve_checks(p, sweeps, fits.tol) if sweep else {}

    # ---- (d) tessellation ---------------------------------------------------------------------------------
    ck.stage(stage="tessellation")
    tess = ck.get("tessellation")
    if tess is None:
        tess = tessellation_report(p, src)
        ck.put("tessellation", tess)
    ck.stage(stage="key_clearances")

    # ---- key clearances: exact minimum distance of selected far-apart pairs (static pose) --------------------
    key = []
    for e in fits.report:
        pa, pb = e["parts"]
        best = None
        for i, a in enumerate(insts):
            if not fnmatch.fnmatchcase(a.kind, pa):
                continue
            for j, b in enumerate(insts):
                if j == i or not fnmatch.fnmatchcase(b.kind, pb):
                    continue
                d, qa, _qb = distance(shapes[i], shapes[j])
                if best is None or d < best[0]:
                    best = (d, a.name, b.name, qa)
        if best:
            key.append({"parts": [pa, pb], "min_mm": best[0], "a": best[1], "b": best[2], "at": best[3],
                        "reason": e.get("reason", "")})

    # ---- summary ------------------------------------------------------------------------------------------
    overl = [s for s in static if s["state"] == "overlap"]
    cont = [s for s in static if s["state"] == "contact"]
    defects = [s for s in static if s["state"] in ("overlap", "contact") and not s["allowed"]]
    sweep_bad = [s for s in sweeps if s.get("clash")]
    summary = {
        "instances": n, "pairs": n_pairs, "bbox_candidates": len(cand),
        "overlaps": len(overl), "overlaps_allowed": sum(1 for s in overl if s["allowed"]),
        "contacts": len(cont), "contacts_allowed": sum(1 for s in cont if s["allowed"]),
        "unexplained_overlaps": sum(1 for s in overl if not s["allowed"]),
        "unexplained_contacts": sum(1 for s in cont if not s["allowed"]),
        "missing_required_contacts": len(missing),
        "floating": len(floating), "floating_without_unlisted_contacts": len(floating_allowed_only),
        "moving_pairs": len(sweeps), "moving_clashes": len(sweep_bad) + len(valve.get("failures", [])),
        "tessellation_fail": sum(1 for t in tess["parts"].values() if not t["ok"]),
        "seconds": time.time() - t_start,
    }
    report = {"note": "MODELLED geometry (no physical engine). CAD units mm.", "tag": tag,
              "source_hash": src_hash, "pose": {"crank_deg": 0.0, "vvt_advance_deg": 0.0},
              "tolerances": tol, "summary": summary, "static": static, "missing_required_contacts": missing,
              "floating": floating, "floating_without_unlisted_contacts": floating_allowed_only,
              "defects": defects, "sweep": sweeps, "valve_checks": valve, "tessellation": tess,
              "key_clearances": key}
    ck.stage(stage="renders" if renders else "writing")
    if renders:
        try:
            from cad.render import audit_views
            report["renders"] = audit_views(tag, out_dir)
        except Exception as exc:  # renders are evidence, not a gate
            report["renders_error"] = repr(exc)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"audit_{tag}.json").write_text(json.dumps(report, indent=1, default=_jsonable))
    (out_dir / f"audit_{tag}.md").write_text(markdown(report))
    ck.stage(stage="done", seconds_this_run=time.time() - t_start)
    print(f"[audit:{tag}] " + json.dumps(summary))
    return report


def _jsonable(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def moving_sweep(fits: Fits, boxes, jobs: int, tag: str, ck: "Checkpoint | None" = None) -> list:
    cfg = dict(fits.sweep)
    cfg.update(contact_mm=fits.tol["contact_mm"], overlap_mm3=fits.tol["overlap_mm3"])
    p, insts, mot, src = _S["p"], _S["insts"], _S["mot"], _S["src"]
    n = len(insts)
    local_boxes = {k: bbox(src[k][0].wrapped) for k in {i.kind for i in insts}}
    # swept boxes and speeds per instance
    speed, sbox = {}, {}
    for k, it in enumerate(insts):
        m = mot[it.group]
        if m.static:
            speed[k], sbox[k] = 0.0, boxes[k]
            continue
        lb = local_boxes[it.kind]
        corners = _corners(lb)
        base_corners = np.array([_apply(it.base, c) for c in corners])
        bb = np.array([np.inf] * 3 + [-np.inf] * 3)
        for th in np.arange(0.0, m.period, cfg["swept_box_step_deg"]):
            pts = np.array([_apply(m.loc(th, 0.0), c) for c in base_corners])
            bb[:3] = np.minimum(bb[:3], pts.min(axis=0))
            bb[3:] = np.maximum(bb[3:], pts.max(axis=0))
        speed[k] = group_speed(p, mot, it.group, base_corners, m.period)
        pad = speed[k] * cfg["swept_box_step_deg"]
        sbox[k] = bb + np.array([-pad] * 3 + [pad] * 3)
    jobs_list = []
    for i in range(n):
        for j in range(i + 1, n):
            A, B = insts[i], insts[j]
            if A.group == B.group or (mot[A.group].static and mot[B.group].static):
                continue
            if _is_spring(A) or _is_spring(B):   # deforming part: checked at full lift in valve_checks()
                continue
            if boxes_gap(sbox[i], sbox[j]) > cfg["near_mm"]:
                continue
            ma, mb = mot[A.group], mot[B.group]
            fit = fits.match(A, B)
            contact_mode = bool(fit and fit.get("moving") == "contact")
            lifter = ma if ma.lift_loc is not None else mb if mb.lift_loc is not None else None
            if lifter is not None and (ma.static or mb.static):
                # valve part against a static part: the pose depends on the lift only -> sweep the lift range
                jobs_list.append((i, j, lifter.u_max, asm.LIFT_SPEED * VEL_SAFETY, contact_mode, 0.0, cfg, sbox[i],
                                  sbox[j], "lift"))
                continue
            period = max(ma.period, mb.period)
            advs = [0.0]
            if (ma.vvt or mb.vvt) and not (ma.vvt and mb.vvt):
                advs.append(_S["vvt_max"])   # a phased group against an unphased one: parked and full advance
            for adv in advs:
                jobs_list.append((i, j, period, speed[i] + speed[j], contact_mode, adv, cfg, sbox[i], sbox[j],
                                  "crank"))
    def jkey(jb):
        return f"{insts[jb[0]].name}|{insts[jb[1]].name}|{jb[9]}|{jb[5]:g}"
    done = ck.load("sweep") if ck else {}
    todo = [jb for jb in jobs_list if jkey(jb) not in done]
    # longest first (valve/cam pairs over 720 deg) so the pool does not end on one long straggler
    todo.sort(key=lambda jb: -jb[2])
    print(f"[audit:{tag}] sweep: {len(jobs_list)} moving pairs, {len(jobs_list) - len(todo)} from the checkpoint",
          flush=True)
    if ck:
        ck.stage(stage="sweep", sweep_total=len(jobs_list), sweep_done=len(jobs_list) - len(todo))
    ctx = mp.get_context("fork")
    fh = ck.appender("sweep") if ck else None
    try:
        with ctx.Pool(jobs) as pool:
            for k, r in enumerate(pool.imap_unordered(_sweep_pair, todo, chunksize=1)):
                key = f"{r['a']}|{r['b']}|{r['param']}|{r['vvt_advance_deg']:g}"
                done[key] = r
                if fh:
                    Checkpoint.add(fh, key, r)
                    ck.stage(sweep_done=len(jobs_list) - len(todo) + k + 1, last=key)
                print(f"[audit:{tag}]   {len(jobs_list) - len(todo) + k + 1}/{len(jobs_list)} {r['a']} / {r['b']}: "
                      f"min {r['min_mm']:.3f} mm ({r['evals']} evals, {r['seconds']:.0f} s)", flush=True)
    finally:
        if fh:
            fh.close()
    res = [done[jkey(jb)] for jb in jobs_list]
    for r in res:
        A = next(i for i in insts if i.name == r["a"])
        B = next(i for i in insts if i.name == r["b"])
        fit = fits.match(A, B)
        r["reason"] = fit["reason"] if fit else None
        if r["mode"]:
            r["clash"] = r["max_overlap_mm3"] > fits.tol["overlap_mm3"] and not (
                fit and fit["type"] == "overlap" and r["max_overlap_mm3"] <= float(fit.get("max_overlap_mm3", math.inf)))
        else:
            touching = r["min_mm"] <= fits.tol["contact_mm"]
            r["clash"] = bool(touching and not fit)
    return res


def valve_checks(p, sweeps: list, tol: dict) -> dict:
    """Spring coil bind and compressed-spring clearances at full lift (the spring is a deforming part, placed at its
    installed length in the assembly), and the CAD piston-to-valve minimum against the M3 model (kin/valvetrain.py)
    at VVT parked and full advance."""
    from cad import valvetrain as vt
    from kin.valvetrain import clearance_map
    out = {"springs": {}, "piston_valve": {}, "failures": []}
    insts, src, mot = _S["insts"], _S["src"], _S["mot"]
    sd = vt.spring_design(p)
    for kind in ("intake", "exhaust"):
        lift = p.max_lift(kind)
        st = vt.stack(p, kind)
        Lc = vt.spring_located_length(p, kind, lift)
        margin = Lc - p.L(sd["solid_height_m"])
        g = asm.valve_group(0, kind, 1)
        frame = mot[g].lift_loc(0.0)
        spring = (asm.Pos(0, 0, st["spring_bottom"]) * vt.spring(p, Lc)).wrapped.Moved(frame.wrapped)
        near = {}
        for it in insts:
            if it.cyl != 1 or not (it.name.endswith(f"_1_{kind}_1") or it.name == "head"):
                continue
            if it.kind.startswith("spring_"):
                continue
            shp = src[it.kind][0].wrapped
            loc = mot[it.group].lift_loc(lift) * it.base if mot[it.group].lift_loc else mot[it.group].loc(0.0) * it.base
            d = distance(spring, shp.Moved(loc.wrapped))[0]
            near[it.name] = d
        lim = p.L(p.vt["spring"]["min_coil_bind_margin"])
        ok = margin >= lim
        bad = {k: v for k, v in near.items() if v <= tol["contact_mm"] and not any(
            k.startswith(a) for a in ("retainer_", "spring_seat_"))}
        out["springs"][kind] = {"compressed_length_mm": Lc, "solid_height_mm": p.L(sd["solid_height_m"]),
                                "coil_bind_margin_mm": margin, "min_margin_mm": lim, "distances_at_full_lift_mm": near,
                                "ok": ok and not bad}
        if not (ok and not bad):
            out["failures"].append(f"spring {kind}: margin {margin:.2f} mm, contacts {bad}")
    # piston-to-valve: CAD sweep minimum vs M3 (flat crown + pocket model along the valve axis). The CAD pocket is a
    # cylinder on the valve axis, valve_relief_margin/2 wider than the head all round, with its floor where M3 puts
    # it: piston material inside the cylinder lies below the floor (>= the M3 gap) and material outside it is at
    # least the wall gap from the head, so min(M3, wall) <= CAD <= M3 (up to the sweep resolution).
    wall = p.L(p.section("piston")["valve_relief_margin"]) / 2
    from kin.cam import design_cams
    prof = design_cams(p.engine, p.vt)
    for kind in ("intake", "exhaust"):
        advs = [0.0] + ([_S["vvt_max"]] if kind == "intake" else [])
        for adv in advs:
            rows = [w for w in sweeps if w.get("param") == "crank" and abs(w["vvt_advance_deg"] - adv) < 1e-9
                    and ((w["a"].startswith(f"valve_") and f"_{kind}_" in w["a"] and w["b"].startswith("piston_"))
                         or (w["b"].startswith("valve_") and f"_{kind}_" in w["b"] and w["a"].startswith("piston_")))]
            _th, gaps = clearance_map(prof[kind], kind, p.engine, p.vt, [math.radians(adv)])
            m3 = p.L(float(gaps.min()))
            cad = min((w["min_mm"] for w in rows), default=None)
            out["piston_valve"][f"{kind}_adv{adv:g}"] = {"cad_min_mm": cad, "m3_min_mm": m3,
                                                        "difference_mm": None if cad is None else cad - m3,
                                                        "pocket_wall_mm": wall,
                                                        "expected_range_mm": [min(m3, wall), m3],
                                                        "pairs": len(rows)}
    return out


def _is_spring(inst) -> bool:
    return inst.kind.startswith("spring_") and not inst.kind.startswith("spring_seat")


def tessellation_report(p, src) -> dict:
    from cad.tess import triangulate
    ex = p.section("export") if "export" in p.cad else {"linear_deflection": 0.001, "angular_deflection_rad": 0.1,
                                                         "relative": True, "max_chord_mm": 0.02,
                                                         "max_facet_angle_rad": 0.1}
    rel = bool(ex.get("relative", False))
    lin, ang = float(ex["linear_deflection"]), float(ex["angular_deflection_rad"])
    if not rel:
        lin = p.L(lin)
    out = {"settings": {"linear_deflection": lin, "relative": rel, "angular_deflection_rad": ang},
           "limits": {"max_chord_mm": float(ex["max_chord_mm"]), "max_facet_angle_rad": float(ex["max_facet_angle_rad"])},
           "parts": {}}
    for kind, (shape, _m) in src.items():
        _mesh, q = triangulate(shape, lin, ang, measure=True, relative=rel)
        q["ok"] = q["chord_max_mm"] <= out["limits"]["max_chord_mm"] and \
            q["facet_angle_max_rad"] <= out["limits"]["max_facet_angle_rad"]
        out["parts"][kind] = q
    vals = list(out["parts"].values())
    out["max_chord_mm"] = max(v["chord_max_mm"] for v in vals)
    out["max_facet_angle_rad"] = max(v["facet_angle_max_rad"] for v in vals)
    out["triangles_total"] = sum(v["triangles"] for v in vals)
    return out


# ------------------------------------------------------------------------------------------------ markdown
def _f(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isinf(x):
            return "inf"
        return f"{x:.{nd}f}"
    return str(x)


def markdown(r: dict) -> str:
    s = r["summary"]
    L = [f"# Redline CAD audit: {r['tag']} (MODELLED geometry)", "",
         f"Source hash `{r['source_hash']}`; static pose crank {r['pose']['crank_deg']} deg, VVT parked. "
         f"Tolerances: contact <= {r['tolerances']['contact_mm']} mm, overlap > {r['tolerances']['overlap_mm3']} mm^3. "
         f"Generated by `python -m cad.audit --tag {r['tag']}` in {s['seconds']:.0f} s.", "",
         "## Summary", "", "| Check | Count |", "|---|---|"]
    for k in ("instances", "pairs", "bbox_candidates", "overlaps", "overlaps_allowed", "unexplained_overlaps",
              "contacts", "contacts_allowed", "unexplained_contacts", "missing_required_contacts", "floating",
              "floating_without_unlisted_contacts", "moving_pairs", "moving_clashes", "tessellation_fail"):
        L.append(f"| {k.replace('_', ' ')} | {s[k]} |")
    L += ["", "## Defects (unexplained overlaps and contacts, static pose)", ""]
    if r["defects"]:
        L += ["| Part A | Part B | State | Overlap mm^3 | Distance mm | At (mm) |", "|---|---|---|---|---|---|"]
        for d in sorted(r["defects"], key=lambda d: -d["overlap_mm3"]):
            pt = "-" if d["point"] is None else "(" + ", ".join(f"{v:.1f}" for v in d["point"]) + ")"
            L.append(f"| {d['a']} | {d['b']} | {d['state']} | {_f(d['overlap_mm3'])} | {_f(d['distance_mm'], 4)} | {pt} |")
    else:
        L.append("None.")
    L += ["", "## Missing required contacts", ""]
    L += [f"- {m['part']} does not touch any {m['expected_contact_with']} ({m['reason']})" for m in r["missing_required_contacts"]] or ["None."]
    L += ["", "## Floating parts (not connected to the block through contacts)", ""]
    L.append(", ".join(r["floating"]) if r["floating"] else "None.")
    if r["floating_without_unlisted_contacts"] != r["floating"]:
        L += ["", "Held only through unlisted (defect) contacts: " +
              ", ".join(sorted(set(r["floating_without_unlisted_contacts"]) - set(r["floating"])))]
    L += ["", "## Moving clearances (full cycle)", ""]
    if r["sweep"]:
        L += ["Param `crank`: full 720 deg cycle (at VVT advance shown); `lift`: valve lift range (at = lift mm).", "",
              "| Moving pair | Param | VVT deg | Mode | Min mm | At | Certified >= mm | Overlap mm^3 | Evals | Result |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for w in sorted(r["sweep"], key=lambda w: w["min_mm"]):
            ov = w.get("max_overlap_mm3", w.get("overlap_at_min_mm3"))
            res = "CLASH" if w["clash"] else ("contact (allowed)" if w["min_mm"] <= r["tolerances"]["contact_mm"] else "clear")
            L.append(f"| {w['a']} / {w['b']} | {w.get('param', 'crank')} | {_f(w.get('vvt_advance_deg', 0.0), 0)} | "
                     f"{'contact' if w['mode'] else 'b&b'} | {_f(w['min_mm'], 3)} | {_f(w['at_deg'], 2)} | "
                     f"{_f(w.get('certified_lb_mm'), 3)} | {_f(ov, 3)} | {w['evals']} | {res} |")
    else:
        L.append("Not run.")
    v = r.get("valve_checks") or {}
    if v:
        L += ["", "## Valve springs at full lift and piston-to-valve vs M3", "",
              "| Spring | Length at full lift mm | Solid mm | Coil-bind margin mm | Limit mm | OK |", "|---|---|---|---|---|---|"]
        for k, x in v["springs"].items():
            L.append(f"| {k} | {x['compressed_length_mm']:.2f} | {x['solid_height_mm']:.2f} | {x['coil_bind_margin_mm']:.2f} | "
                     f"{x['min_margin_mm']:.2f} | {'yes' if x['ok'] else 'NO'} |")
        L += ["", "| Valve / VVT | CAD min mm | M3 model min mm | CAD - M3 mm | Pairs |", "|---|---|---|---|---|"]
        for k, x in v["piston_valve"].items():
            L.append(f"| {k} | {_f(x['cad_min_mm'])} | {_f(x['m3_min_mm'])} | {_f(x['difference_mm'])} | {x['pairs']} |")
    t = r["tessellation"]
    L += ["", "## Tessellation", "",
          f"Settings: linear deflection {t['settings']['linear_deflection']} ({'relative' if t['settings']['relative'] else 'absolute, mm'}), "
          f"angular {t['settings']['angular_deflection_rad']} rad. Limits: chord <= {t['limits']['max_chord_mm']} mm, "
          f"facet angle <= {t['limits']['max_facet_angle_rad']} rad. Worst: chord {t['max_chord_mm']:.4f} mm, facet "
          f"{t['max_facet_angle_rad']:.4f} rad, {t['triangles_total']} triangles over all source solids.", "",
          "| Part | Triangles | Max chord mm | Max facet angle rad | OK |", "|---|---|---|---|---|"]
    for k, v in t["parts"].items():
        L.append(f"| {k} | {v['triangles']} | {v['chord_max_mm']:.4f} | {v['facet_angle_max_rad']:.4f} | {'yes' if v['ok'] else 'NO'} |")
    L += ["", "## Static pairs (touching bounding boxes)", "",
          "| Part A | Part B | State | Overlap mm^3 | Distance mm | Allowed | Reason |", "|---|---|---|---|---|---|---|"]
    for d in sorted(r["static"], key=lambda d: (d["state"] != "overlap", d["state"] != "contact", d["a"], d["b"])):
        L.append(f"| {d['a']} | {d['b']} | {d['state']} | {_f(d['overlap_mm3'])} | {_f(d['distance_mm'], 4)} | "
                 f"{'-' if d['allowed'] is None else ('yes' if d['allowed'] else 'NO')} | {d['reason'] or ''} |")
    if r.get("renders"):
        L += ["", "## Renders", ""] + [f"![{k}]({Path(v).name})" for k, v in r["renders"].items()]
    return "\n".join(L) + "\n"


def _kind_family(kind: str) -> str:
    return re.sub(r"_(\d+)$", "", kind)


SIMPLIFICATIONS = (
    "Threads are plain cylinders at the minor/major diameter with a line-to-line fit (listed as contacts), no helix.",
    "No head water jacket, oil drain-backs, cam cover or timing cover; the block jacket is an open-deck annulus.",
    "Valve springs are static solids at the installed length in the assembly; the audit rebuilds each spring at the "
    "lift of its valve for the static pose and checks coil bind and neighbour distances at full lift instead of "
    "sweeping the deforming coil.",
    "Head bolts sit under the cam axes in closed pockets below the journals (bolts cannot be removed with the cams "
    "fitted, as on some production heads); the pocket roof is head_bolt_pocket_roof below the journal bore.",
    "Intake: swan-neck runners curl over the head into a plenum above the cam area; the 550 mm gasdyn runner "
    "length (port + flange + runner to the bellmouth lip) is met exactly, and there is no cam cover modelled under "
    "the plenum. The throttle body is held by flange contact (no flange bolts modelled) and its blade is fixed at "
    "the rendered angle.",
    "Exhaust: 4-2-1 long-tube header. Primaries leave the ports, bend down beside the block (pair 1+4 in the inboard "
    "plane, 2+3 outboard) and merge at x = 0 beside the skirt; the secondaries bend rearward and run beside the sump "
    "to a collector whose axis is 10 mm below the sump floor, outlet flange facing rearward. Every duct centre line "
    "is solved to the gasdyn length (primaries 450 mm flange face to merge mouth, secondaries 450 mm cone outlet to "
    "collector mouth; tests/test_cad_audit.py checks them) with bends >= 1.5 x tube OD. 450 + 450 mm cannot be "
    "folded inside the block length, so the collector and outlet sit behind the flywheel (under where the gearbox "
    "would be; no bellhousing, starter or downpipe is modelled). One pressed-steel stay ties the 1+4 secondary to "
    "the skirt rib of bulkhead 4. Merges are plain cones (no formed merge spike); the gasdyn collector volume is "
    "not matched by the CAD collector.",
    "Bearings: spec/dynamics.yaml now uses the CAD journal diameters (main 52 mm, big end 45 mm, spec/cad.yaml "
    "crankshaft; tests/test_m5_dynamics.py::test_bearing_diameters_match_cad keeps them equal), so the oil-film model "
    "and the 3D shells describe the same bearings. The shells are plain half tubes (no oil grooves, holes or "
    "locating tangs).",
    "12 clearance pairs keep a constant gap while they move: crank vs conrod and conrod cap (1.0 mm side gaps) and "
    "conrod vs piston (0.5 mm small-end side gap). The faces involved are planes normal to the common rotation "
    "axis, so the gap is exactly constant and the grid minimum is the true minimum; the Lipschitz certificate "
    "cannot see that (it bounds the surface speed, not the gap's derivative) and reports a lower bound of 0 for "
    "these pairs (Certified column).",
    "Running-contact pairs (bearing shells, sliding fits, cam lobe on bucket) are checked for penetration by exact "
    "common volume on a subsampled grid (spec/cad_fits.yaml sweep.contact_volume_step_deg/_mm), not certified between "
    "samples; clearance pairs carry the Lipschitz certificate.",
    "Belts, the oil pan gasket and sealants are not modelled as compliant parts; the belt is a rigid band at its "
    "static path.",
)


DOC_RENDERS = (
    ("assembly_front34", "full assembly, 3/4 front (intake side): swan-neck runners, plenum, throttle body, stay"),
    ("assembly_rear34_exhaust", "full assembly, 3/4 rear (exhaust side): 4-2-1 header and its stay, flywheel"),
    ("exhaust_clearance", "header past the oil pan, pan bolts and flywheel, from below (gaps in the key clearances table above)"),
    ("section_cyl1", "exact section through cylinder 1 near full intake lift"),
    ("valvetrain_cutaway", "head sectioned through cylinder 1's rear valve pair: cams, caps, buckets, springs, guides, seats"),
    ("intake_cutaway", "plenum sectioned through the runner axes: bellmouths inside the plenum, runner curl, port"),
    ("bottom_end_pan_removed", "bottom end with the pan removed: main caps, rods, crank"),
    ("exploded", "exploded assembly"),
)


def docs_markdown(final: dict, baseline: dict | None) -> str:
    """docs/CAD-AUDIT.md: method, baseline vs final, every touching-box part pair aggregated by part kinds (with the
    reason of each allowed contact), the moving-pair minima, valve checks and tessellation."""
    s = final["summary"]
    L = ["# Redline CAD audit (MODELLED geometry)", "",
         "All numbers describe the parametric CAD model; there is no physical engine. Regenerate with "
         "`python -m cad.audit --tag final` (writes out/cad/audit_final.{json,md}) and `python -m cad.audit --docs`.",
         "", "## Method", "",
         "* **Static clash** (crank 0 deg, VVT parked): every pair of the "
         f"{s['instances']} placed part instances ({s['pairs']} pairs) is screened with exact (optimal) bounding boxes; "
         f"the {s['bbox_candidates']} pairs whose boxes come within the contact tolerance get an exact BRep boolean "
         "common volume and an exact minimum distance (BRepExtrema). Overlap > "
         f"{final['tolerances']['overlap_mm3']} mm^3 = overlap, distance <= {final['tolerances']['contact_mm']} mm = "
         "contact. Every contact/overlap must be listed in spec/cad_fits.yaml with its build reason; anything else is "
         "a defect.",
         "* **Moving sweep**: every pair with at least one moving part whose swept boxes come within "
         "5 mm is swept. Crank-driven pairs: 2 deg grid over the full period (720 deg when a cam or valve is involved, at "
         "VVT parked and at full intake advance when only one side is phased), local refinement of the minimum down to "
         "0.01 deg, then a branch and bound with a rigorous Lipschitz bound (largest surface speed of the bodies) "
         "certifies a lower bound. Valve parts against static parts are swept over the valve lift (0.5 mm grid, refined) "
         "instead of the crank angle. Running-contact pairs (bearings, sliding fits, cam on bucket) are checked for "
         "penetration (common volume) instead of clearance.",
         "* **Support graph**: union-find over the contact/overlap pairs, rooted at the block; a part not connected is "
         "floating (and connections through unlisted contacts are reported separately).",
         "* **Tessellation**: every source solid is meshed with the export settings and the chordal deviation of every "
         "triangle from the exact surface (exact point projection of the worst samples) and the largest angle between "
         "exact vertex normals inside a triangle are measured. OCCT's mesher alone does not honour the angular "
         "deflection across tori and swept B-spline faces, so cad/tess.py refines every curved face afterwards on the "
         "exact surface (red-green splits at surface midpoints) until the chord and facet-angle targets hold.", ""]
    if baseline:
        b = baseline["summary"]
        L += ["## Baseline (before the polish) vs final", "", "| Check | Baseline | Final |", "|---|---|---|"]
        for k in ("instances", "bbox_candidates", "overlaps", "unexplained_overlaps", "contacts", "unexplained_contacts",
                  "floating", "moving_pairs", "moving_clashes", "tessellation_fail"):
            L.append(f"| {k.replace('_', ' ')} | {b.get(k, '-')} | {s.get(k, '-')} |")
        bt, ft = baseline["tessellation"], final["tessellation"]
        L.append(f"| worst chord mm | {bt['max_chord_mm']:.4f} | {ft['max_chord_mm']:.4f} |")
        L.append(f"| worst facet angle rad | {bt['max_facet_angle_rad']:.4f} | {ft['max_facet_angle_rad']:.4f} |")
        L += ["", "Baseline defects: " + ("; ".join(f"{d['a']} / {d['b']} {d['state']} at "
                                                    f"{_f(d['distance_mm'], 3)} mm" for d in baseline["defects"]) or "none")
              + ". Baseline floating: " + (", ".join(baseline["floating"]) or "none") + "; held only through defect "
              "contacts: " + (", ".join(sorted(set(baseline["floating_without_unlisted_contacts"]) - set(baseline["floating"])))
                              or "none") + ".", ""]
    L += ["## Final summary", "", "| Check | Count |", "|---|---|"]
    for k, v in s.items():
        if k != "seconds":
            L.append(f"| {k.replace('_', ' ')} | {v} |")
    agg = {}
    for d in final["static"]:
        key = tuple(sorted((_kind_family(d["kind_a"]), _kind_family(d["kind_b"])))) + (d["state"],)
        a = agg.setdefault(key, {"n": 0, "dmin": math.inf, "vmax": 0.0, "allowed": True, "reason": d["reason"]})
        a["n"] += 1
        a["dmin"] = min(a["dmin"], d["distance_mm"])
        a["vmax"] = max(a["vmax"], d["overlap_mm3"])
        a["allowed"] = a["allowed"] and d["allowed"] is not False
    L += ["", "## Static part pairs (every pair with touching bounding boxes, aggregated by part kind)", "",
          f"All other {s['pairs'] - s['bbox_candidates']} pairs have exact bounding boxes more than "
          f"{final['tolerances']['contact_mm']} mm apart (clear).", "",
          "| Part kinds | State | Pairs | Min distance mm | Max overlap mm^3 | Allowed | Reason |", "|---|---|---|---|---|---|---|"]
    order = {"overlap": 0, "contact": 1, "clear": 2}
    for (ka, kb, st), a in sorted(agg.items(), key=lambda kv: (order[kv[0][2]], kv[0][0], kv[0][1])):
        allowed = "-" if st == "clear" else ("yes" if a["allowed"] else "NO")
        L.append(f"| {ka} / {kb} | {st} | {a['n']} | {_f(a['dmin'], 3)} | {_f(a['vmax'], 3)} | {allowed} | "
                 f"{a['reason'] or ''} |")
    mv = {}
    for w in final["sweep"]:
        key = (re.sub(r"_\d.*$", "", w["a"]), re.sub(r"_\d.*$", "", w["b"]), w.get("param", "crank"),
               w.get("vvt_advance_deg", 0.0), w["mode"])
        a = mv.setdefault(key, {"n": 0, "min": math.inf, "where": None, "cert": math.inf, "clash": False,
                                "reason": w.get("reason")})
        a["n"] += 1
        if w["min_mm"] < a["min"]:
            a["min"], a["where"] = w["min_mm"], f"{w['a']} / {w['b']} at {_f(w['at_deg'], 2)}"
        if w.get("certified_lb_mm") is not None:
            a["cert"] = min(a["cert"], w["certified_lb_mm"])
        a["clash"] = a["clash"] or bool(w.get("clash"))
    if final.get("key_clearances"):
        L += ["", "## Key clearances (selected pairs, exact minimum distance in the static pose)", "",
              "Pairs whose boxes are far apart never reach the exact check above; these are measured anyway because "
              "they matter for packaging (spec/cad_fits.yaml report_clearances). Rotating partners (flywheel, ring "
              "gear) are axisymmetric, so the static distance holds over the cycle. The header-to-pan, header-to-pan-bolt "
              "and header-to-flywheel gaps are shown in docs/renders/exhaust_clearance.png.", "",
              "| Pair (kinds) | Min distance mm | Closest instances | Why it matters |", "|---|---|---|---|"]
        for k in final["key_clearances"]:
            L.append(f"| {k['parts'][0]} / {k['parts'][1]} | {_f(k['min_mm'], 1)} | {k['a']} / {k['b']} | {k['reason']} |")
    L += ["", "## Moving pairs (minimum over the cycle or the lift range, aggregated by part family)", "",
          "| Family A | Family B | Param | VVT deg | Mode | Pairs | Min mm | Worst pair at (deg or lift mm) | Certified >= mm | Result |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for (fa, fb, prm, adv, mode), a in sorted(mv.items(), key=lambda kv: kv[1]["min"]):
        res = "CLASH" if a["clash"] else ("running contact, no penetration" if mode else "clear")
        cert = "-" if math.isinf(a["cert"]) else _f(a["cert"], 3)
        L.append(f"| {fa} | {fb} | {prm} | {adv:g} | {'contact' if mode else 'clearance'} | {a['n']} | {_f(a['min'], 3)} | "
                 f"{a['where']} | {cert} | {res} |")
    v = final.get("valve_checks") or {}
    if v:
        L += ["", "## Valve springs and piston-to-valve", ""]
        for k, x in v["springs"].items():
            L.append(f"* {k} spring at full lift: {x['compressed_length_mm']:.2f} mm vs solid {x['solid_height_mm']:.2f} mm, "
                     f"coil-bind margin {x['coil_bind_margin_mm']:.2f} mm (limit {x['min_margin_mm']:.2f} mm); nearest "
                     "parts at full lift: " + ", ".join(f"{n} {d:.2f} mm" for n, d in sorted(
                         x["distances_at_full_lift_mm"].items(), key=lambda t: t[1])[:4]) + ".")
        L += ["", "The CAD minimum is the 3D distance over the 720-deg sweep; M3 is the gap along the valve axis to the "
              "pocket floor. The pocket wall is valve_relief_margin/2 from the head all round, so the CAD minimum must "
              "lie in [min(M3, wall), M3]: below M3 when the head passes the pocket rim or wall before the floor.", "",
              "| Valve / VVT advance | CAD min mm | M3 model min mm | Pocket wall mm | Expected range mm |",
              "|---|---|---|---|---|"]
        for k, x in v["piston_valve"].items():
            lo, hi = x.get("expected_range_mm", [None, None])
            L.append(f"| {k} | {_f(x['cad_min_mm'])} | {_f(x['m3_min_mm'])} | {_f(x.get('pocket_wall_mm'))} | "
                     f"{_f(lo)} to {_f(hi)} |")
    t = final["tessellation"]
    L += ["", "## Tessellation", "",
          f"Export settings: absolute chordal deflection {t['settings']['linear_deflection']} mm, angular "
          f"{t['settings']['angular_deflection_rad']} rad. Worst measured chordal deviation {t['max_chord_mm']:.4f} mm "
          f"(limit {t['limits']['max_chord_mm']}), worst facet angle {t['max_facet_angle_rad']:.4f} rad (limit "
          f"{t['limits']['max_facet_angle_rad']}), {t['triangles_total']} triangles over all source solids. glTF "
          "exports carry exact per-vertex surface normals (smooth shading)."]
    if baseline:
        bt = baseline["tessellation"]
        L += ["", f"Baseline (build123d exporters, relative deflection {bt['settings']['linear_deflection']}, "
              f"{bt['settings']['angular_deflection_rad']} rad): worst chord {bt['max_chord_mm']:.4f} mm, facet angle "
              f"{bt['max_facet_angle_rad']:.4f} rad, no exported normals."]
    L += ["", "## Renders (docs/renders/, `python -m cad.render`)", ""] + [f"* `{n}.png`: {c}" for n, c in DOC_RENDERS]
    L += ["", "## Simplifications and known gaps (honest list)", ""] + [f"* {x}" for x in SIMPLIFICATIONS]
    return "\n".join(L) + "\n"


DOCS_JSON = ROOT / "docs" / "cad-audit.json"   # committed compact copy (out/ is not versioned)
SWEEP_KEYS = ("a", "b", "param", "vvt_advance_deg", "mode", "min_mm", "at_deg", "certified_lb_mm", "max_overlap_mm3",
              "clash")


def compact(r: dict) -> dict:
    """The audit result without the per-pair static records of clear pairs and without sweep internals: what
    tests/test_cad_audit.py checks and docs/CAD-AUDIT.md reports, small enough to commit."""
    keep = ("note", "tag", "source_hash", "pose", "tolerances", "summary", "defects", "missing_required_contacts",
            "floating", "floating_without_unlisted_contacts", "valve_checks", "tessellation", "key_clearances")
    out = {k: r[k] for k in keep if k in r}
    out["static_contacts"] = [{k: d[k] for k in ("a", "b", "state", "distance_mm", "overlap_mm3", "allowed", "reason")
                               if k in d} for d in r["static"] if d["state"] != "clear"]
    out["sweep"] = [{k: w.get(k) for k in SWEEP_KEYS} for w in r["sweep"]]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tag", default="final")
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--no-renders", action="store_true")
    ap.add_argument("--no-sweep", action="store_true")
    ap.add_argument("--fresh", action="store_true", help="discard the checkpoint in out/cad/audit_<tag>.ckpt/")
    ap.add_argument("--docs", action="store_true",
                    help="only write docs/CAD-AUDIT.md and docs/cad-audit.json from out/cad/audit_*.json")
    a = ap.parse_args(argv)
    if a.docs:
        final = json.loads((OUT_DIR / "audit_final.json").read_text())
        bp = OUT_DIR / "audit_baseline.json"
        base = json.loads(bp.read_text()) if bp.exists() else None
        (ROOT / "docs" / "CAD-AUDIT.md").write_text(docs_markdown(final, base))
        DOCS_JSON.write_text(json.dumps(compact(final), indent=1) + "\n")
        return 0
    r = run(a.tag, a.jobs, renders=not a.no_renders, sweep=not a.no_sweep, fresh=a.fresh)
    s = r["summary"]
    bad = s["unexplained_overlaps"] + s["unexplained_contacts"] + s["floating"] + s["moving_clashes"] + s["tessellation_fail"]
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
