"""CalculiX (ccx) input writer, runner and .frd result reader.

Units are SI throughout (m, N, Pa, K, W). Node and element ids in the deck are 1-based
(index + 1). A `Deck` collects model data (mesh, material, sets, couplings, boundary conditions)
and a list of steps; each step replaces the previous loads (OP=NEW) and requests nodal results
(U, S or NT) in the .frd file, which `read_frd` returns per step.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import os
import shutil
import subprocess

import numpy as np

from fea.mesh import Mesh

CCX_ENV = "REDLINE_CCX"
SOLVER_THREADS = 1


def ccx_executable() -> str | None:
    exe = os.environ.get(CCX_ENV) or shutil.which("ccx")
    return exe


@dataclass
class Material:
    name: str
    E: float
    nu: float
    rho: float = 0.0
    alpha: float = 0.0          # 1/K
    k: float = 0.0              # W/(m K)


@dataclass
class Step:
    kind: str = "static"                        # static | heat
    cloads: list = field(default_factory=list)  # (node_idx, dof, value)
    dloads: list = field(default_factory=list)  # (elem_idx, face_no, pressure)  pressure + into face
    centrif: tuple | None = None                # (omega^2, point, axis)
    temperatures: np.ndarray | None = None      # (N,) nodal temperatures for thermal stress
    films: list = field(default_factory=list)   # (elem_idx, face_no, T_sink, h)
    dfluxes: list = field(default_factory=list) # (elem_idx, face_no, q)  heat flux into the body
    fixed_temps: list = field(default_factory=list)  # (node_idx, T)
    boundaries: list = field(default_factory=list)   # step-local (node_idx, dof, value)


@dataclass
class Deck:
    mesh: Mesh
    material: Material
    boundaries: list = field(default_factory=list)   # (node_idx, dof_first, dof_last, value)
    extra_nodes: list = field(default_factory=list)  # (xyz) reference nodes appended after mesh nodes
    couplings: list = field(default_factory=list)    # (ref_node_idx, face_sel_indices, kind)
    equations: list = field(default_factory=list)    # [(node, dof, coef), ...] = 0
    steps: list = field(default_factory=list)
    t_ref: float = 0.0
    solver: str = "SPOOLES"                          # or "ITERATIVE CHOLESKY"

    def add_ref_node(self, xyz) -> int:
        self.extra_nodes.append(np.asarray(xyz, dtype=float))
        return self.mesh.n_nodes + len(self.extra_nodes) - 1

    def write(self, path: Path):
        m = self.mesh
        mat = self.material
        lines = ["*HEADING", "Redline M6 FEA (modelled)", "*NODE"]
        allnodes = m.nodes if not self.extra_nodes else np.vstack([m.nodes, np.array(self.extra_nodes)])
        lines += [f"{i + 1},{x:.12e},{y:.12e},{z:.12e}" for i, (x, y, z) in enumerate(allnodes)]
        lines.append("*ELEMENT, TYPE=C3D10, ELSET=EALL")
        t = m.tets + 1
        lines += [f"{e + 1}," + ",".join(str(v) for v in row) for e, row in enumerate(t)]
        lines += [f"*MATERIAL, NAME={mat.name}", "*ELASTIC", f"{mat.E:.6e},{mat.nu:.6f}"]
        if mat.rho:
            lines += ["*DENSITY", f"{mat.rho:.6e}"]
        if mat.alpha:
            lines += [f"*EXPANSION, ZERO={self.t_ref:.6f}", f"{mat.alpha:.6e}"]
        if mat.k:
            lines += ["*CONDUCTIVITY", f"{mat.k:.6e}"]
        lines.append(f"*SOLID SECTION, ELSET=EALL, MATERIAL={mat.name}")
        bf = m.boundary_faces()
        for ci, (ref, sel, kind) in enumerate(self.couplings):
            if kind.lower() == "rigid":   # sel = node ids; ref = (ref node, rot node)
                lines.append(f"*NSET, NSET=RB{ci}")
                lines += [str(n + 1) for n in sel]
                lines.append(f"*RIGID BODY, NSET=RB{ci}, REF NODE={ref[0] + 1}, ROT NODE={ref[1] + 1}")
                continue
            lines.append(f"*SURFACE, NAME=SC{ci}, TYPE=ELEMENT")
            lines += [f"{bf[0][s] + 1},S{bf[1][s]}" for s in sel]
            lines.append(f"*COUPLING, REF NODE={ref + 1}, SURFACE=SC{ci}, CONSTRAINT NAME=C{ci}")
            lines.append(f"*{kind.upper()}")
            lines.append("1,3")
        if self.equations:
            lines.append("*EQUATION")
            for eq in self.equations:
                lines.append(str(len(eq)))
                lines.append(",".join(f"{n + 1},{d},{c:.12e}" for n, d, c in eq))
        if any(s.kind == "heat" for s in self.steps) or any(s.temperatures is not None for s in self.steps):
            lines += ["*INITIAL CONDITIONS, TYPE=TEMPERATURE"]
            lines += [f"{i + 1},{self.t_ref:.6f}" for i in range(len(allnodes))]
        if self.boundaries:
            lines.append("*BOUNDARY")
            lines += [f"{n + 1},{a},{b},{v:.12e}" for n, a, b, v in self.boundaries]
        for st in self.steps:
            lines += self._step(st, bf)
        path.write_text("\n".join(lines) + "\n")

    def _step(self, st: Step, bf) -> list:
        out = ["*STEP"]
        if st.kind == "heat":
            out += ["*HEAT TRANSFER, STEADY STATE", "1.,1."]
            if st.fixed_temps:
                out.append("*BOUNDARY")
                out += [f"{n + 1},11,11,{T:.6f}" for n, T in st.fixed_temps]
            if st.films:
                out.append("*FILM, OP=NEW")
                out += [f"{bf_e + 1},F{f},{Ts:.6f},{h:.6e}" for bf_e, f, Ts, h in st.films]
            if st.dfluxes:
                out.append("*DFLUX, OP=NEW")
                out += [f"{e + 1},S{f},{q:.6e}" for e, f, q in st.dfluxes]
            out += ["*NODE FILE", "NT", "*END STEP"]
            return out
        out += [f"*STATIC, SOLVER={self.solver}", "*CLOAD, OP=NEW"]
        out += [f"{n + 1},{d},{v:.12e}" for n, d, v in st.cloads]
        out.append("*DLOAD, OP=NEW")
        out += [f"{e + 1},P{f},{p:.12e}" for e, f, p in st.dloads]
        if st.centrif is not None:
            w2, pt, ax = st.centrif
            out.append(f"EALL,CENTRIF,{w2:.12e},{pt[0]:.9e},{pt[1]:.9e},{pt[2]:.9e},{ax[0]:.9e},{ax[1]:.9e},{ax[2]:.9e}")
        if st.boundaries:
            out.append("*BOUNDARY, OP=NEW")
            out += [f"{n + 1},{d},{d},{v:.12e}" for n, d, v in st.boundaries]
        if st.temperatures is not None:
            out.append("*TEMPERATURE, OP=NEW")
            out += [f"{i + 1},{T:.6f}" for i, T in enumerate(st.temperatures)]
            # reference nodes (if any) keep the reference temperature
            out += [f"{self.mesh.n_nodes + j + 1},{self.t_ref:.6f}" for j in range(len(self.extra_nodes))]
        else:
            if any(s.temperatures is not None for s in self.steps):
                out.append("*TEMPERATURE, OP=NEW")
                out += [f"{i + 1},{self.t_ref:.6f}" for i in range(self.mesh.n_nodes + len(self.extra_nodes))]
        out += ["*NODE FILE", "U,RF", "*EL FILE", "S", "*END STEP"]
        return out


def run(deck: Deck, workdir: Path, name: str = "job", threads: int = 4, timeout: float = 3600.0) -> list[dict]:
    """Write <workdir>/<name>.inp, run ccx, return read_frd results (one dict per step)."""
    exe = ccx_executable()
    if exe is None:
        raise RuntimeError("CalculiX (ccx) not found; install calculix-ccx or set REDLINE_CCX")
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / f"{name}.inp"
    deck.write(inp)
    # The multi-threaded SPOOLES factorisation in ccx 2.21 occasionally returns a wrong solution (observed:
    # identical decks, different displacements); the equation solver therefore runs single-threaded,
    # element/stress loops keep `threads`.
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), CCX_NPROC_EQUATION_SOLVER=str(SOLVER_THREADS),
               CCX_NPROC_STIFFNESS=str(threads), CCX_NPROC_RESULTS=str(threads))
    proc = subprocess.run([exe, "-i", name], cwd=workdir, env=env, capture_output=True, text=True,
                          timeout=timeout)
    log = proc.stdout + proc.stderr
    (workdir / f"{name}.log").write_text(log)
    if proc.returncode != 0 or "*ERROR" in log:
        raise RuntimeError(f"ccx failed for {inp}:\n" + "\n".join(log.splitlines()[-30:]))
    n_total = deck.mesh.n_nodes + len(deck.extra_nodes)
    return read_frd(workdir / f"{name}.frd", n_total)


def read_frd(path: Path, n_nodes: int) -> list[dict]:
    """Nodal result blocks of an ASCII .frd file: [{'DISP': (N,3), 'STRESS': (N,6), 'NDTEMP': (N,),
    'FORC': (N,3)}, ...] grouped by step. Stress order: xx, yy, zz, xy, yz, zx."""
    steps: list[dict] = []
    cur_name, cur_step, data, ncomp = None, None, None, 0
    with open(path) as f:
        for line in f:
            if line.startswith("    1PSTEP"):
                continue
            if line.startswith("  100C"):
                cur_step = int(line[58:68]) if len(line) > 68 and line[58:68].strip() else None
                continue
            if line.startswith(" -4"):
                parts = line.split()
                cur_name, ncomp = parts[1], int(parts[2])
                width = {"DISP": 3, "STRESS": 6, "NDTEMP": 1, "FORC": 3}.get(cur_name, ncomp)
                data = np.zeros((n_nodes, width))
                continue
            if cur_name is None:
                continue
            if line.startswith(" -5"):
                continue
            if line.startswith(" -1"):
                node = int(line[3:13])
                vals = [float(line[13 + 12 * i: 25 + 12 * i]) for i in range(data.shape[1])]
                if node <= n_nodes:
                    data[node - 1] = vals
                continue
            if line.startswith(" -3"):
                key = cur_name
                if not steps or key in steps[-1]:
                    steps.append({})
                steps[-1][key] = data[:, 0] if data.shape[1] == 1 else data
                cur_name = None
    return steps


def von_mises(s: np.ndarray) -> np.ndarray:
    """Von Mises stress of (..., 6) tensors (xx, yy, zz, xy, yz, zx)."""
    xx, yy, zz, xy, yz, zx = (s[..., i] for i in range(6))
    return np.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2) + 3 * (xy**2 + yz**2 + zx**2))


def principal(s: np.ndarray) -> np.ndarray:
    """Principal stresses (..., 3), ascending."""
    xx, yy, zz, xy, yz, zx = (s[..., i] for i in range(6))
    T = np.stack([np.stack([xx, xy, zx], -1), np.stack([xy, yy, yz], -1), np.stack([zx, yz, zz], -1)], -2)
    return np.linalg.eigvalsh(T)
