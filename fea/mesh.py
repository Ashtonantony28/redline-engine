"""gmsh meshing to second-order tetrahedra (C3D10) and boundary-face bookkeeping.

`mesh_step` imports a STEP file (CAD units, mm per spec/cad.yaml units) and returns a `Mesh` in SI
metres (ADR 0001). Boundary faces are extracted from the tet connectivity, so load and support
regions are picked geometrically on the mesh (centroid / normal predicates) instead of relying on
fragile CAD face indices.

Node ordering: gmsh's 10-node tetrahedron numbers the edge nodes (0-1, 1-2, 2-0, 3-0, 3-2, 3-1)
while CalculiX/Abaqus C3D10 uses (1-2, 2-3, 3-1, 1-4, 2-4, 3-4): the last two edge nodes swap
(`GMSH_TO_CCX_TET10`). `Mesh.check_midside` verifies the mapping on the actual elements.
Face numbering follows CalculiX for C3D10: S1 = (1,2,3), S2 = (1,4,2), S3 = (2,4,3), S4 = (3,4,1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

GMSH_TET10 = 11                        # gmsh element type id of the 10-node tetrahedron
GMSH_TO_CCX_TET10 = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
# corner and mid-side node positions (0-based, CCX order) of each C3D10 face S1..S4
FACE_CORNERS = np.array([[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]])
FACE_MIDS = np.array([[4, 5, 6], [7, 8, 4], [8, 9, 5], [6, 9, 7]])
OPPOSITE = np.array([3, 2, 0, 1])      # vertex not on face k
MIN_SCALED_JACOBIAN = 0.01            # reject meshes with (nearly) inverted curved tets
MIDSIDE_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]


@dataclass
class Mesh:
    nodes: np.ndarray            # (N, 3) m
    tets: np.ndarray             # (E, 10) 0-based node indices, CalculiX C3D10 order
    meta: dict = field(default_factory=dict)
    _bf: tuple | None = None

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_elems(self) -> int:
        return len(self.tets)

    def check_midside(self, tol: float = 0.05, quantile: float = 0.5) -> float:
        """Quantile of |mid-node - edge midpoint| / edge length over all edges. Interior edges are straight
        (offset 0); only boundary edges on curved faces bend, so a wrong permutation (offset ~0.5 on almost
        every edge) is caught while a few strongly curved boundary edges are tolerated."""
        offs = []
        for k, (a, b) in enumerate(MIDSIDE_EDGES):
            pa, pb, pm = self.nodes[self.tets[:, a]], self.nodes[self.tets[:, b]], self.nodes[self.tets[:, 4 + k]]
            offs.append(np.linalg.norm(pm - 0.5 * (pa + pb), axis=1) / np.linalg.norm(pb - pa, axis=1))
        q = float(np.quantile(np.concatenate(offs), quantile))
        if q > tol:
            raise ValueError(f"C3D10 mid-side node mapping inconsistent ({quantile:.0%} offset {q:.3f})")
        return q

    def boundary_faces(self):
        """(elem, face_no 1..4, corner node ids (F,3), all 6 node ids (F,6), centroid (F,3),
        unit outward normal (F,3), area (F,))."""
        if self._bf is not None:
            return self._bf
        E = self.n_elems
        corners = self.tets[:, FACE_CORNERS]                  # (E, 4, 3)
        keys = np.sort(corners.reshape(-1, 3), axis=1)
        _, inv, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
        once = counts[inv.ravel()] == 1
        idx = np.nonzero(once)[0]
        elem, face = idx // 4, idx % 4
        cn = corners.reshape(-1, 3)[idx]
        mids = self.tets[elem[:, None], FACE_MIDS[face]]
        six = np.concatenate([cn, mids], axis=1)
        p = self.nodes[cn]
        c = p.mean(axis=1)
        n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        area = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2 * area[:, None])
        opp = self.nodes[self.tets[elem, OPPOSITE[face]]]
        flip = np.einsum("ij,ij->i", n, c - opp) < 0
        n[flip] *= -1
        self._bf = (elem, face + 1, cn, six, c, n, area)
        return self._bf

    def select_faces(self, pred) -> np.ndarray:
        """Indices into boundary_faces() for which pred(centroid (F,3), normal (F,3)) is True."""
        _e, _f, _cn, _six, c, n, _a = self.boundary_faces()
        return np.nonzero(pred(c, n))[0]

    def face_nodes(self, sel) -> np.ndarray:
        """Unique node ids (all 6 per face) of the selected boundary faces."""
        return np.unique(self.boundary_faces()[3][sel].ravel())

    def select_nodes(self, pred) -> np.ndarray:
        return np.nonzero(pred(self.nodes))[0]

    def volume(self) -> float:
        p = self.nodes[self.tets[:, :4]]
        return float(np.abs(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                                      p[:, 3] - p[:, 0])).sum() / 6.0)


def _gmsh_options(gmsh, h_max: float, h_min: float, curvature_pts: int, threads: int):
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curvature_pts)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)   # untangles curved 2nd-order elements
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("General.NumThreads", 1)   # single-threaded: reproducible meshes
    gmsh.option.setNumber("Mesh.RandomSeed", 1)


GMSH_OVERRIDES_ENV = "REDLINE_GMSH_OPTIONS"   # "Name=value;Name=value" (debugging aid)


def mesh_step(step_path: str | Path, length_scale: float, h_max: float, h_min: float,
              curvature_pts: int, threads: int = 1, refine=None) -> Mesh:
    """Mesh a STEP solid with 2nd-order tets.

    length_scale: metres per CAD unit (1/cad_length_per_m). h_max/h_min in CAD units.
    refine: optional list of distance-threshold refinements (see below).
    """
    import gmsh

    gmsh.initialize(["-noenv"])
    try:
        _gmsh_options(gmsh, h_max, h_min, curvature_pts, threads)
        import os
        for kv in filter(None, os.environ.get(GMSH_OVERRIDES_ENV, "").split(";")):
            k, v = kv.split("=")
            gmsh.option.setNumber(k.strip(), float(v))
        gmsh.model.add("part")
        gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()
        vols = gmsh.model.getEntities(3)
        if len(vols) > 1:  # fuse stray multi-volume imports into one conforming solid
            gmsh.model.occ.fragment(vols[:1], vols[1:])
            gmsh.model.occ.synchronize()
        if refine:
            # each entry: {"expr": distance expression in x, y, z (CAD units, gmsh MathEval syntax),
            #              "size": size at distance <= dist_min, "dist_min", "dist_max"}; size grows
            # linearly to h_max at dist_max (gmsh Threshold field)
            fields = []
            for r in refine:
                fd = gmsh.model.mesh.field.add("MathEval")
                gmsh.model.mesh.field.setString(fd, "F", r["expr"])
                ft = gmsh.model.mesh.field.add("Threshold")
                gmsh.model.mesh.field.setNumber(ft, "InField", fd)
                gmsh.model.mesh.field.setNumber(ft, "SizeMin", r["size"])
                gmsh.model.mesh.field.setNumber(ft, "SizeMax", h_max)
                gmsh.model.mesh.field.setNumber(ft, "DistMin", r["dist_min"])
                gmsh.model.mesh.field.setNumber(ft, "DistMax", r["dist_max"])
                fields.append(ft)
            fmin = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        nodes = coords.reshape(-1, 3) * length_scale
        etypes, _etags, enodes = gmsh.model.mesh.getElements(3)
        conn = None
        min_sj = float("nan")
        for t, en, et in zip(etypes, enodes, _etags):
            if t == GMSH_TET10:
                min_sj = float(np.min(gmsh.model.mesh.getElementQualities(et, "minSJ")))
        for t, en in zip(etypes, enodes):
            if t == GMSH_TET10:
                conn = np.asarray(en, dtype=np.int64).reshape(-1, 10)
        if conn is None:
            raise RuntimeError("gmsh produced no 10-node tetrahedra")
    finally:
        gmsh.finalize()
    # compact node numbering to the nodes used by the tets
    lookup = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    lookup[tags.astype(np.int64)] = np.arange(len(tags))
    conn = lookup[conn][:, GMSH_TO_CCX_TET10]
    used = np.unique(conn)
    remap = np.full(len(tags), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    if not min_sj > MIN_SCALED_JACOBIAN:
        raise RuntimeError(f"invalid curved elements: min scaled Jacobian {min_sj:.3g}")
    m = Mesh(nodes[used], remap[conn], meta={"h_max": h_max, "h_min": h_min, "curvature_pts": curvature_pts,
                                             "min_scaled_jacobian": min_sj})
    m.check_midside()
    return m
