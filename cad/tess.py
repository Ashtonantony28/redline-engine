"""Tessellation, mesh export (glTF binary / STL) and tessellation-quality metrics (CAD units, mm).

build123d's exporters mesh with a *relative* linear deflection and reuse any triangulation already cached on the
shape, so the delivered chordal error is neither known nor controlled. Here every export is meshed afresh with an
ABSOLUTE linear deflection (mm) and an angular deflection (rad) from spec/cad.yaml `export`, and the result is
measured (`quality`): chordal deviation of every triangle from the exact BRep surface and the angle between the
normals of adjacent triangles (what shows up as faceting on silhouettes). glTF files carry per-vertex normals
evaluated on the exact surface (smooth shading inside each BRep face, hard edges between faces).
"""
from __future__ import annotations

import math
import struct
from dataclasses import dataclass

import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepLib import BRepLib_ToolTriangulatedShape
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.GeomAbs import GeomAbs_Plane
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Pnt

MM_PER_M = 1000.0
Y_UP_QUAT = [-math.sqrt(0.5), 0.0, 0.0, math.sqrt(0.5)]   # -90 deg about x: CAD Z-up -> glTF Y-up
REFINE_TOP = 8            # per face: worst UV-midpoint samples re-measured by exact point projection
INT16_MAX = 32767
INT8_MAX = 127
REFINE_MAX_PASSES = 40    # per-face bisection passes (each splits the worst edge of every failing triangle)
EDGE_KEY = np.int64(1) << 31   # edge key a * EDGE_KEY + b (vertex ids < 2^31)
REFINE_MIN_EDGE_MM = 1e-3  # never split an edge shorter than this (singular points: poles, apexes)


@dataclass
class Mesh:
    """Triangle soup with per-vertex normals; CAD units (mm). face_id labels the BRep face of each triangle."""
    pos: np.ndarray       # (n, 3) float64
    nrm: np.ndarray       # (n, 3) float64, unit
    tri: np.ndarray       # (m, 3) int64
    face_id: np.ndarray   # (m,) int64

    @property
    def n_tri(self) -> int:
        return len(self.tri)


def _fresh_copy(shape):
    """Deep copy of the TopoDS shape (no cached triangulation shared with the caller's shape)."""
    w = shape.wrapped if hasattr(shape, "wrapped") else shape
    cp = BRepBuilderAPI_Copy(w, True, False).Shape()
    BRepTools.Clean_s(cp)
    return cp


def triangulate(shape, lin_tol: float, ang_tol: float, measure: bool = False, relative: bool = False,
                refine: bool | None = None):
    """Mesh `shape` with absolute deflection lin_tol (mm) and angular deflection ang_tol (rad).

    OCCT's mesher honours the linear deflection on most faces but not the angular one across the narrow direction of
    tori and swept B-spline faces (fillets, helical spring wire, swept ports), so with refine (default: on for an
    absolute deflection) every curved face's triangulation is refined afterwards on the exact surface (`_refine_face`)
    until each triangle's estimated chordal deviation is <= lin_tol and its facet angle is <= ang_tol.

    relative=True reproduces build123d's exporters (deflection relative to each edge's size). Returns Mesh, plus (with measure=True) a dict of quality metrics: max chordal deviation (mm) and max angle
    the exact surface normal turns across one triangle (rad)."""
    if refine is None:
        refine = not relative
    w = _fresh_copy(shape)
    BRepMesh_IncrementalMesh(w, lin_tol, relative, ang_tol, True)
    pos, nrm, tri, fid = [], [], [], []
    base = 0
    chord_max, chord_at = 0.0, None
    exp = TopExp_Explorer(w, TopAbs_FACE)
    k = 0
    while exp.More():
        face = TopoDS.Face(exp.Current())
        exp.Next()
        loc = TopLoc_Location()
        t = BRep_Tool.Triangulation_s(face, loc)
        if t is None:
            k += 1
            continue
        if not t.HasNormals():
            BRepLib_ToolTriangulatedShape.ComputeNormals_s(face, t)
        trsf = loc.Transformation()
        n = t.NbNodes()
        P = np.empty((n, 3))
        N = np.empty((n, 3))
        for i in range(1, n + 1):
            pt = t.Node(i).Transformed(trsf)
            P[i - 1] = (pt.X(), pt.Y(), pt.Z())
            d = t.Normal(i).Transformed(trsf)
            N[i - 1] = (d.X(), d.Y(), d.Z())
        T = np.array([t.Triangle(i).Get() for i in range(1, t.NbTriangles() + 1)], dtype=np.int64) - 1
        if face.Orientation() == TopAbs_REVERSED:
            N = -N
            T = T[:, ::-1]
        UV = np.array([(t.UVNode(i).X(), t.UVNode(i).Y()) for i in range(1, n + 1)]) if t.HasUVNodes() else None
        if refine and UV is not None and len(T) and BRepAdaptor_Surface(face).GetType() != GeomAbs_Plane:
            P, N, T, UV = _refine_face(BRep_Tool.Surface_s(face), P, N, T, UV, lin_tol, ang_tol)
            n = len(P)
        if measure and len(T):
            c, at = _face_chord(face, UV, P, T)
            if c > chord_max:
                chord_max, chord_at = c, at
        pos.append(P)
        nrm.append(N)
        tri.append(T + base)
        fid.append(np.full(len(T), k))
        base += n
        k += 1
    mesh = Mesh(np.vstack(pos) if pos else np.zeros((0, 3)), np.vstack(nrm) if nrm else np.zeros((0, 3)),
                np.vstack(tri) if tri else np.zeros((0, 3), dtype=np.int64),
                np.concatenate(fid) if fid else np.zeros(0, dtype=np.int64))
    if not measure:
        return mesh
    return mesh, {"chord_max_mm": chord_max, "chord_at_mm": chord_at, "facet_angle_max_rad": _adjacent_angle(mesh),
                  "triangles": mesh.n_tri, "vertices": len(mesh.pos)}


def _tri_angle(N, T):
    """Per triangle: largest angle (rad) between the unit normals at its three vertices."""
    a = np.zeros(len(T))
    for i, j in ((0, 1), (1, 2), (2, 0)):
        a = np.maximum(a, np.arccos(np.clip(np.einsum("ij,ij->i", N[T[:, i]], N[T[:, j]]), -1.0, 1.0)))
    return a


class _Cache:
    """Surface samples keyed by int64 (edge or triangle keys), so unchanged edges are evaluated once."""

    def __init__(self):
        self.keys = np.zeros(0, np.int64)
        self.pts = np.zeros((0, 3))
        self.nrm = np.zeros((0, 3))

    def get(self, keys, uv_fn, nref_fn, ev):
        idx = np.searchsorted(self.keys, keys)
        idx_c = np.minimum(idx, max(len(self.keys) - 1, 0))
        have = (len(self.keys) > 0) & (self.keys[idx_c] == keys) if len(self.keys) else np.zeros(len(keys), bool)
        miss = ~have
        if miss.any():
            p, n = ev(uv_fn(miss), nref_fn(miss))
            k = np.concatenate([self.keys, keys[miss]])
            o = np.argsort(k, kind="stable")
            self.keys, self.pts, self.nrm = k[o], np.vstack([self.pts, p])[o], np.vstack([self.nrm, n])[o]
            idx_c = np.searchsorted(self.keys, keys)
        return self.pts[idx_c], self.nrm[idx_c]


def _refine_face(surf, P, N, T, UV, chord_tol: float, ang_tol: float):
    """Adaptive bisection of one face's triangulation on its exact surface.

    Each pass evaluates the surface at the UV midpoint of every edge and the UV centroid of every triangle; a triangle
    is marked when one of those surface points is further than chord_tol from the corresponding point of the flat
    triangle (an upper bound of the projection distance) or when its vertex normals differ by more than ang_tol. Each
    marked triangle splits its worst edge (largest deviation or normal turn relative to its limit; the longest edge
    for a centroid-only failure) at the exact surface midpoint (point and normal). Neighbours sharing a split edge are
    split into two or three triangles, so the mesh stays conforming inside the face and anisotropic faces (thin
    swept wire, fillets) are only refined across their curved direction. The face boundary stays on the BRep edges."""
    P, N, T, UV = P.copy(), N.copy(), T.copy(), UV.copy()

    def ev(uv, nref):
        pts = np.empty((len(uv), 3))
        nrm = np.empty((len(uv), 3))
        for k, (u, v) in enumerate(uv):
            pr = GeomLProp_SLProps(surf, float(u), float(v), 1, 1e-9)
            q = pr.Value()
            pts[k] = (q.X(), q.Y(), q.Z())
            if pr.IsNormalDefined():
                d = pr.Normal()
                nv = np.array((d.X(), d.Y(), d.Z()))
                nrm[k] = nv if nv @ nref[k] >= 0.0 else -nv
            else:
                nrm[k] = nref[k]
        ln = np.linalg.norm(nrm, axis=1)
        nrm[ln > 0] /= ln[ln > 0, None]
        return pts, nrm

    ecache = _Cache()
    cen_bad = None                                   # per triangle: centroid chord failure (None = not evaluated)
    for _ in range(REFINE_MAX_PASSES):
        E = np.sort(np.stack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]], axis=1), axis=2)   # (m, 3, 2)
        ue, inv = np.unique(E.reshape(-1, 2), axis=0, return_inverse=True)
        inv = inv.reshape(-1, 3)
        a, b = ue[:, 0], ue[:, 1]
        mid_p, mid_n = ecache.get(a * EDGE_KEY + b, lambda m: 0.5 * (UV[a[m]] + UV[b[m]]), lambda m: N[a[m]] + N[b[m]], ev)
        length = np.linalg.norm(P[a] - P[b], axis=1)
        splittable = length > REFINE_MIN_EDGE_MM
        dev = np.linalg.norm(mid_p - 0.5 * (P[a] + P[b]), axis=1)
        turn = np.arccos(np.clip(np.einsum("ij,ij->i", N[a], N[b]), -1.0, 1.0))
        score = np.where(splittable, np.maximum(dev / chord_tol, turn / ang_tol), -1.0)
        if cen_bad is None:
            cen_bad = np.ones(len(T), bool)
            fresh = np.ones(len(T), bool)
        else:
            fresh = np.isnan(cen_bad_f)
        if fresh.any():
            cp, _ = ev(UV[T[fresh]].mean(axis=1), N[T[fresh]].sum(axis=1))
            cen_bad[fresh] = np.linalg.norm(cp - P[T[fresh]].mean(axis=1), axis=1) > chord_tol
        bad_c = cen_bad
        tri_score = score[inv]                                   # (m, 3)
        bad_t = (tri_score.max(axis=1) > 1.0) | ((_tri_angle(N, T) > ang_tol) | bad_c) & splittable[inv].any(axis=1)
        if not bad_t.any():
            break
        pick = np.where(tri_score.max(axis=1) > 1.0, tri_score.argmax(axis=1),
                        np.where(splittable[inv], length[inv], -1.0).argmax(axis=1))
        split = np.zeros(len(ue), bool)
        split[inv[bad_t, pick[bad_t]]] = True
        split &= splittable
        if not split.any():
            break
        newid = np.full(len(ue), -1)
        newid[split] = len(P) + np.arange(split.sum())
        P = np.vstack([P, mid_p[split]])
        N = np.vstack([N, mid_n[split]])
        UV = np.vstack([UV, 0.5 * (UV[a[split]] + UV[b[split]])])
        M = newid[inv]                       # (m, 3) midpoint of edge (0-1), (1-2), (2-0); -1 if not split
        s = M >= 0
        cnt = s.sum(axis=1)
        keep = cnt == 0
        out = [T[keep]]
        r = cnt == 3
        A, B, C, m0, m1, m2 = T[r, 0], T[r, 1], T[r, 2], M[r, 0], M[r, 1], M[r, 2]
        out += [np.stack(x, axis=1) for x in ((A, m0, m2), (m0, B, m1), (m2, m1, C), (m0, m1, m2))]
        for e in range(3):
            g = (cnt == 1) & s[:, e]          # split edge e (vertex e -> e+1)
            A, B, C, m = T[g, e], T[g, (e + 1) % 3], T[g, (e + 2) % 3], M[g, e]
            out += [np.stack((A, m, C), axis=1), np.stack((m, B, C), axis=1)]
            g = (cnt == 2) & ~s[:, e]         # edge e unsplit: rotate it to (2 -> 0), split edges are then 0 and 1
            A, B, C = T[g, (e + 1) % 3], T[g, (e + 2) % 3], T[g, e]
            m0, m1 = M[g, (e + 1) % 3], M[g, (e + 2) % 3]
            out += [np.stack((m0, B, m1), axis=1), np.stack((A, m0, m1), axis=1), np.stack((A, m1, C), axis=1)]
        T = np.vstack(out)
        cen_bad_f = np.full(len(T), np.nan)          # kept triangles first (out[0]); new ones are re-evaluated
        cen_bad_f[:keep.sum()] = cen_bad[keep]
        cen_bad = np.zeros(len(T), bool)
        cen_bad[:keep.sum()] = cen_bad_f[:keep.sum()].astype(bool)
    return P, N, T, UV


def _face_chord(face, UV, Pl, T):
    """Max distance from triangle centroids / edge midpoints to the exact surface (global coordinates)."""
    if BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane:
        return 0.0, None
    surf = BRep_Tool.Surface_s(face)
    samples = [Pl[T].mean(axis=1), 0.5 * (Pl[T[:, 0]] + Pl[T[:, 1]]), 0.5 * (Pl[T[:, 1]] + Pl[T[:, 2]]),
               0.5 * (Pl[T[:, 2]] + Pl[T[:, 0]])]
    pts = np.vstack(samples)
    if UV is not None:
        uvs = np.vstack([UV[T].mean(axis=1), 0.5 * (UV[T[:, 0]] + UV[T[:, 1]]), 0.5 * (UV[T[:, 1]] + UV[T[:, 2]]),
                         0.5 * (UV[T[:, 2]] + UV[T[:, 0]])])
        est = np.empty(len(pts))
        for i, (u, v) in enumerate(uvs):
            s = surf.Value(float(u), float(v))
            est[i] = math.dist((s.X(), s.Y(), s.Z()), pts[i])
        order = np.argsort(est)[::-1][:REFINE_TOP]
    else:
        order = np.arange(len(pts))
    best, at = 0.0, None
    for i in order:
        # the extremum search can miss the nearest foot point on swept B-spline faces; the UV-mapped surface point is
        # an exact upper bound of the distance as well, so the smaller of the two is reported
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*pts[i]), surf)
        d = pr.LowerDistance() if pr.NbPoints() else math.inf
        if UV is not None:
            d = min(d, float(est[i]))
        if not math.isfinite(d):
            continue
        if d > best:
            best, at = d, pts[i].tolist()
    return best, at


def _adjacent_angle(mesh: Mesh) -> float:
    """Max angle (rad) the exact surface normal turns across one triangle: the largest angle between the normals at
    its three vertices (evaluated on the BRep surface). This is the facet angle a flat-shaded or silhouette view
    shows, and unlike triangle-normal differences it is not polluted by sliver triangles on narrow faces."""
    if mesh.n_tri == 0:
        return 0.0
    return float(_tri_angle(mesh.nrm, mesh.tri).max())


# ------------------------------------------------------------------------------------------------ writers
def write_stl(mesh: Mesh, path) -> None:
    P, T = mesh.pos, mesh.tri
    fn = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    ln = np.linalg.norm(fn, axis=1)
    fn[ln > 0] /= ln[ln > 0, None]
    rec = np.zeros(len(T), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = fn, P[T]
    with open(path, "wb") as f:
        f.write(b"redline binary STL (mm)".ljust(80, b" "))
        f.write(struct.pack("<I", len(T)))
        f.write(rec.tobytes())


def write_glb(meshes: dict, path, quantize: bool = False) -> int:
    """Write {name: Mesh} (CAD frame, mm) as one glTF binary, metres, root node rotated Z-up -> Y-up.

    quantize=True stores positions as normalized int16 (dequantized by the mesh node's scale/translation) and
    normals as normalized int8 (KHR_mesh_quantization); otherwise float32. Returns the file size in bytes."""
    import pygltflib as g

    blob = bytearray()
    views, accessors, gl_meshes, nodes = [], [], [], []

    def add_view(data: bytes, target=None, stride=None):
        while len(blob) % 4:
            blob.append(0)
        views.append(g.BufferView(buffer=0, byteOffset=len(blob), byteLength=len(data), target=target,
                                  byteStride=stride))
        blob.extend(data)
        return len(views) - 1

    children = []
    for name, m in meshes.items():
        P = m.pos / MM_PER_M
        lo, hi = P.min(axis=0), P.max(axis=0)
        node_kw = {}
        if quantize:
            ctr, half = 0.5 * (lo + hi), np.maximum(0.5 * (hi - lo), 1e-9)
            q = np.round((P - ctr) / half * INT16_MAX).astype("<i2")
            qp = np.zeros((len(P), 4), dtype="<i2")
            qp[:, :3] = q
            pv = add_view(qp.tobytes(), g.ARRAY_BUFFER, 8)
            accessors.append(g.Accessor(bufferView=pv, componentType=g.SHORT, normalized=True, count=len(P),
                                        type=g.VEC3, min=[-INT16_MAX / INT16_MAX] * 3, max=[1.0] * 3))
            node_kw = {"translation": ctr.tolist(), "scale": half.tolist()}
            qn = np.zeros((len(P), 4), dtype="<i1")
            qn[:, :3] = np.round(m.nrm * INT8_MAX).astype("<i1")
            nv = add_view(qn.tobytes(), g.ARRAY_BUFFER, 4)
            accessors.append(g.Accessor(bufferView=nv, componentType=g.BYTE, normalized=True, count=len(P),
                                        type=g.VEC3))
        else:
            pv = add_view(P.astype("<f4").tobytes(), g.ARRAY_BUFFER)
            accessors.append(g.Accessor(bufferView=pv, componentType=g.FLOAT, count=len(P), type=g.VEC3,
                                        min=lo.tolist(), max=hi.tolist()))
            nv = add_view(m.nrm.astype("<f4").tobytes(), g.ARRAY_BUFFER)
            accessors.append(g.Accessor(bufferView=nv, componentType=g.FLOAT, count=len(P), type=g.VEC3))
        idx_t, ct = ("<u2", g.UNSIGNED_SHORT) if len(P) < 65535 else ("<u4", g.UNSIGNED_INT)
        iv = add_view(m.tri.astype(idx_t).ravel().tobytes(), g.ELEMENT_ARRAY_BUFFER)
        accessors.append(g.Accessor(bufferView=iv, componentType=ct, count=m.tri.size, type=g.SCALAR))
        a0 = len(accessors) - 3
        gl_meshes.append(g.Mesh(name=name, primitives=[g.Primitive(attributes=g.Attributes(POSITION=a0, NORMAL=a0 + 1),
                                                                   indices=a0 + 2)]))
        nodes.append(g.Node(name=name, mesh=len(gl_meshes) - 1, **node_kw))
        children.append(len(nodes) - 1)
    nodes.append(g.Node(name="root", rotation=Y_UP_QUAT, children=children))
    gltf = g.GLTF2(asset=g.Asset(generator="redline cad.tess"), scene=0, scenes=[g.Scene(nodes=[len(nodes) - 1])],
                   nodes=nodes, meshes=gl_meshes, accessors=accessors, bufferViews=views,
                   buffers=[g.Buffer(byteLength=len(blob))])
    if quantize:
        gltf.extensionsUsed = ["KHR_mesh_quantization"]
        gltf.extensionsRequired = ["KHR_mesh_quantization"]
    gltf.set_binary_blob(bytes(blob))
    gltf.save_binary(str(path))
    import os
    return os.path.getsize(path)


def transformed(mesh: Mesh, loc) -> Mesh:
    """Mesh moved by a build123d Location (rigid: normals rotated)."""
    T = loc.wrapped.Transformation()
    R = np.array([[T.Value(i, j) for j in (1, 2, 3)] for i in (1, 2, 3)])
    t = np.array([T.Value(i, 4) for i in (1, 2, 3)])
    return Mesh(mesh.pos @ R.T + t, mesh.nrm @ R.T, mesh.tri.copy(), mesh.face_id.copy())


def merge(meshes: list) -> Mesh:
    """One triangle soup from several meshes (indices offset)."""
    off, P, N, T, F = 0, [], [], [], []
    for m in meshes:
        P.append(m.pos)
        N.append(m.nrm)
        T.append(m.tri + off)
        F.append(m.face_id)
        off += len(m.pos)
    return Mesh(np.vstack(P), np.vstack(N), np.vstack(T), np.concatenate(F))
