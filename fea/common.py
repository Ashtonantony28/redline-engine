"""Shared machinery for the part analyses: mesh-level options, bore bearing pressures, ROI masks,
nodal Goodman evaluation and plotting. ADR 0001 units (SI) inside, CAD units (mm) only for meshing."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from fea import ccx
from fea.fatigue import goodman_sf, signed_von_mises
from fea.mesh import Mesh
from redline.spec import ROOT

PLANE_TOL = 1.0e-7    # m, node-on-plane test (mesh coordinates are exact copies of the CAD planes)


@dataclass
class LevelResult:
    level: int
    n_nodes: int
    n_elems: int
    peak: float                      # monitored peak stress [Pa] for the convergence study
    sf: dict                         # named safety factors
    data: dict = field(default_factory=dict)   # reported numbers
    mesh: Mesh | None = None
    fields: dict = field(default_factory=dict) # nodal fields for plots


def mesh_options(mspec: dict, level: int) -> tuple:
    h = float(mspec["h_max"][level])
    return h, h * float(mspec["h_min_ratio"]), int(mspec["curvature_pts"][level])


def workdir(fspec: dict, name: str) -> Path:
    return ROOT / fspec["solver"]["workdir"] / name


def bore_faces(mesh: Mesh, axis_point, radius: float, side_dir, tol: float):
    """Boundary faces on a bore along x through axis_point (y, z) of `radius`, on the half facing
    side_dir (unit (y, z)); returns (face indices, cos(angle from side_dir))."""
    _e, _f, _cn, _six, c, n, _a = mesh.boundary_faces()
    dy, dz = c[:, 1] - axis_point[0], c[:, 2] - axis_point[1]
    r = np.hypot(dy, dz)
    cosang = (dy * side_dir[0] + dz * side_dir[1]) / np.maximum(r, 1e-30)
    sel = np.nonzero((np.abs(r - radius) < tol) & (cosang > 0.0))[0]
    return sel, cosang[sel]


def cosine_pressure(mesh: Mesh, sel, weight, force_dir, force: float):
    """Pressures p_i = p0 * weight_i on faces sel, scaled so that the resultant along force_dir (unit
    xyz) equals `force` (pressure acts against the outward normal). Returns [(elem, face, p)]."""
    e, f, _cn, _six, _c, n, a = mesh.boundary_faces()
    w = np.asarray(weight, dtype=float)
    unit = -(w * a[sel])[:, None] * n[sel]
    resultant = unit.sum(axis=0) @ np.asarray(force_dir, dtype=float)
    p0 = force / resultant
    return [(int(e[s]), int(f[s]), float(p0 * wi)) for s, wi in zip(sel, w)], float(p0)


def plane_nodes(mesh: Mesh, axis: int, value: float = 0.0):
    return np.nonzero(np.abs(mesh.nodes[:, axis] - value) < PLANE_TOL)[0]


def roi_mask(mesh: Mesh, exclude_points, radius: float) -> np.ndarray:
    """True for nodes farther than `radius` from every point in exclude_points (St Venant zone)."""
    mask = np.ones(mesh.n_nodes, dtype=bool)
    for pt in exclude_points:
        mask &= np.linalg.norm(mesh.nodes - np.asarray(pt), axis=1) > radius
    return mask


def surface_mask(mesh: Mesh) -> np.ndarray:
    m = np.zeros(mesh.n_nodes, dtype=bool)
    m[np.unique(mesh.boundary_faces()[3].ravel())] = True
    return m


def goodman_nodes(s_alt: np.ndarray, s_mean: np.ndarray, se, sut):
    """Nodal Goodman SF from alternating and mean stress tensors (N, 6)."""
    sa = np.abs(signed_von_mises(s_alt))
    sm = signed_von_mises(s_mean)
    return goodman_sf(sa, sm, se, sut), sa, sm


def plot_field(mesh: Mesh, values: np.ndarray, path: Path, title: str, label: str, view=(20, -60),
               mirror=None):
    """Boundary triangles coloured by a nodal field (face mean), saved as PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    _e, _f, cn, _six, _c, _n, _a = mesh.boundary_faces()
    tri = mesh.nodes[cn] * 1.0e3
    val = values[cn].mean(axis=1)
    fig = plt.figure(figsize=(7.5, 6.0), facecolor="#fcfcfb")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#fcfcfb")
    cmap = plt.get_cmap("viridis")
    vmin, vmax = float(np.percentile(values, 1)), float(values.max())
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    tris = [tri]
    vals = [val]
    if mirror:
        for ax_idx in mirror:
            new = []
            for t in tris:
                tt = t.copy()
                tt[..., ax_idx] *= -1
                new.append(tt)
            tris += new
            vals += vals[: len(new)]
    for t, v in zip(tris, vals):
        pc = Poly3DCollection(t, facecolors=cmap(norm(v)), edgecolors="none", linewidths=0)
        ax.add_collection3d(pc)
    allp = np.concatenate([t.reshape(-1, 3) for t in tris])
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    ctr, span = (lo + hi) / 2, (hi - lo).max() / 2
    ax.set_xlim(ctr[0] - span, ctr[0] + span)
    ax.set_ylim(ctr[1] - span, ctr[1] + span)
    ax.set_zlim(ctr[2] - span, ctr[2] + span)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.view_init(*view)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.7)
    cb.set_label(label)
    ax.set_title(title + "\n(MODELLED, no physical engine)", fontsize=10, loc="left")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, facecolor="#fcfcfb")
    plt.close(fig)


def solve(deck: ccx.Deck, fspec: dict, name: str, tag: str):
    return ccx.run(deck, workdir(fspec, name), tag, threads=int(fspec["solver"]["threads"]))
