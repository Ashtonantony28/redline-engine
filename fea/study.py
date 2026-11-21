"""Mesh-convergence study and design-iteration loop (same pattern as drone-design's frame loop).

For each part: analyse at the medium mesh; if a safety factor misses its target, grow the part's
iterate parameters in spec/cad.yaml (in memory) by `step`, keep the M1 mass budget (rebalancing where the
part names a rebalance parameter), and repeat. The accepted geometry then gets the full three-level
convergence study; the reported SFs come from the finest mesh. If the finest mesh misses a target the
loop resumes. `apply_to_spec` writes accepted parameter values back to spec/cad.yaml (reviewed change).
"""
from __future__ import annotations

import importlib
import json
import re
import time
from dataclasses import dataclass, field

from fea.config import CAD_SPEC_PATH, fea_spec, with_overrides
from fea.fatigue import richardson

PARTS = ("conrod", "crankshaft", "piston", "main_cap")
ITERATE_LEVEL = 1
MASS_BUDGET_TOL = 0.03          # M1 tolerance (tests/test_m1_cad.py BUDGET_TOL)


def part_module(name: str):
    return importlib.import_module(f"fea.parts.{name}")


def targets_for(name: str, fspec: dict) -> dict:
    t = fspec["targets"]
    out = {"fatigue": t["fatigue_sf"]}
    if name == "conrod":
        out["buckling"] = t["buckling_sf"]
    if name == "main_cap":
        out["bolt_fatigue"] = t["fatigue_sf"]
    if name == "piston":
        out["static"] = t["piston_static_sf"]
    return out


def passes(sf: dict, targets: dict) -> bool:
    return all(sf[k] >= v for k, v in targets.items())


def part_loads(name: str, p, fspec):
    from fea import loads as L
    if name == "conrod":
        return L.conrod_loads(p, fspec)
    if name == "crankshaft":
        return L.crank_histories(p, fspec)
    if name == "piston":
        return L.piston_loads(fspec)
    return L.main_bearing_peak(p, fspec)


def mass_state(p) -> dict:
    from fea.loads import cad_masses
    m = cad_masses(p)
    b = p.engine["masses"]
    return {"piston_assembly_kg": m["piston_assembly"], "conrod_kg": m["conrod"],
            "piston_rel_err": m["piston_assembly"] / b["piston_assembly"] - 1.0,
            "conrod_rel_err": m["conrod"] / b["conrod"] - 1.0, "reciprocating_kg": m["reciprocating"]}


def part_mass(p, which: str) -> float:
    from cad.conrod import conrod
    from cad.massprops import mass_props
    from cad.piston import piston_parts
    if which == "conrod":
        return mass_props(conrod(p), p.density("conrod"), p.k).mass
    return sum(mass_props(s, p.density(m), p.k).mass for s, m in piston_parts(p).values())


REBALANCE_STEPS = 24            # bisection steps (relative resolution ~ 1e-7 of the bracket)
REBALANCE_BRACKET = (0.25, 2.0) # search range as multiples of the current value


def rebalance(p, section: str, param: str, target_mass: float, which: str):
    """Bisection on `param` so the CAD mass of `which` returns to target_mass (kg); mass grows with param."""
    lo, hi = (f * p.cad[section][param] for f in REBALANCE_BRACKET)
    for _ in range(REBALANCE_STEPS):
        mid = 0.5 * (lo + hi)
        m = part_mass(with_overrides(p, section, {param: mid}), which)
        lo, hi = (mid, hi) if m < target_mass else (lo, mid)
    return with_overrides(p, section, {param: round(0.5 * (lo + hi), 6)})


@dataclass
class PartStudy:
    name: str
    targets: dict
    iterations: list = field(default_factory=list)
    levels: list = field(default_factory=list)
    final_params: dict = field(default_factory=dict)
    changed: dict = field(default_factory=dict)
    convergence: dict = field(default_factory=dict)
    finest: object = None
    seconds: float = 0.0

    @property
    def sf(self) -> dict:
        return self.finest.sf if self.finest else {}

    @property
    def passed(self) -> bool:
        # the piston also needs its thermal state consistent with the thermo wall temperature (ADR 0006)
        thermal_ok = self.finest.data.get("thermal", {}).get("consistent", True) if self.finest else False
        return passes(self.sf, self.targets) and self.convergence.get("ok", False) and thermal_ok


def convergence_summary(levels: list, fspec: dict) -> dict:
    peaks = [r.peak for r in levels]
    sizes = [r.n_nodes for r in levels]
    change = abs(peaks[-1] - peaks[-2]) / abs(peaks[-1])
    # mean refinement ratio in element size from the node counts (3D: N ~ h^-3, two refinement steps)
    ratio = (sizes[-1] / sizes[0]) ** (1.0 / 6.0) if len(sizes) > 2 else None
    p_obs, f_ext, gci = richardson(peaks, ratio) if len(peaks) == 3 else (None, None, None)
    ok = change < fspec["convergence"]["max_change"]
    return {"peaks_Pa": peaks, "nodes": sizes, "change_finest_two": change, "ok": ok, "refinement_ratio": ratio,
            "richardson_order": p_obs, "richardson_extrapolated_Pa": f_ext, "gci_fine": gci,
            "sf_by_level": [r.sf for r in levels]}


def run_part(name: str, p, fspec=None, log=print) -> PartStudy:
    fspec = fspec or fea_spec()
    mod = part_module(name)
    fd = fspec[mod.SECTION]
    it_spec = fd["iterate"]
    targets = targets_for(name, fspec)
    st = PartStudy(name, targets)
    st.stop_reason = ""
    t0 = time.time()
    mass0 = {"conrod": part_mass(p, "conrod"), "piston_assembly": part_mass(p, "piston_assembly")}
    cur = p
    n_iter = 0
    while True:
        loads = part_loads(name, cur, fspec)
        r = mod.analyse(cur, fspec, ITERATE_LEVEL, loads, tag=f"_it{n_iter}")
        entry = {"iter": n_iter, "params": {k: cur.cad[mod.SECTION][k] for k in it_spec["params"]},
                 "sf": r.sf, "peak_Pa": r.peak, "level": ITERATE_LEVEL, "mass": mass_state(cur)}
        if "rebalance" in it_spec:
            entry["params"][it_spec["rebalance"]] = cur.cad[mod.SECTION][it_spec["rebalance"]]
        st.iterations.append(entry)
        log(f"  [{name}] iteration {n_iter}: SF {fmt_sf(r.sf)} (targets {targets})")
        nxt = None
        if not passes(r.sf, targets):
            if n_iter >= it_spec["max_iter"]:
                st.stop_reason = "iteration limit reached before the SF target"
            else:
                nxt = next_design(cur, mod.SECTION, it_spec, mass0)
                if nxt is None:
                    st.stop_reason = (f"mass budget exhausted: rebalancing {it_spec['rebalance']} would go below "
                                      f"{it_spec['rebalance_min']} m")
        if nxt is not None:
            cur, n_iter = nxt, n_iter + 1
            continue
        # accepted (or final) design: three-level convergence study, SFs from the finest mesh
        levels = []
        for lev in range(fspec["convergence"]["levels"]):
            rr = r if lev == ITERATE_LEVEL else mod.analyse(cur, fspec, lev, loads, tag=f"_it{n_iter}")
            levels.append(rr)
            log(f"  [{name}] level {lev}: {rr.n_nodes} nodes, peak {rr.peak / 1e6:.1f} MPa, SF {fmt_sf(rr.sf)}")
        st.levels, st.finest = levels, levels[-1]
        st.convergence = convergence_summary(levels, fspec)
        if passes(st.finest.sf, targets) or st.stop_reason or n_iter >= it_spec["max_iter"]:
            break
        nxt = next_design(cur, mod.SECTION, it_spec, mass0)     # medium passed, finest did not
        if nxt is None:
            st.stop_reason = "mass budget exhausted after the convergence study"
            break
        cur, n_iter = nxt, n_iter + 1
    st.final_params = {k: cur.cad[mod.SECTION][k] for k in st.iterations[-1]["params"]}
    st.changed = {k: v for k, v in st.final_params.items() if v != p.cad[mod.SECTION][k]}
    st.mass_final = mass_state(cur)
    st.mass_ok = abs(st.mass_final["piston_rel_err"]) <= MASS_BUDGET_TOL and \
        abs(st.mass_final["conrod_rel_err"]) <= MASS_BUDGET_TOL
    st.seconds = time.time() - t0
    return st


def next_design(cur, section, it_spec, mass0):
    """Grow the iterate parameters by one step and rebalance the mass; None if the rebalance is infeasible."""
    upd = {k: round(cur.cad[section][k] + it_spec["step"], 6) for k in it_spec["params"]}
    nxt = with_overrides(cur, section, upd)
    if "rebalance" in it_spec:
        which = it_spec["rebalance_mass"]
        nxt = rebalance(nxt, section, it_spec["rebalance"], mass0[which], which)
        if nxt.cad[section][it_spec["rebalance"]] < it_spec["rebalance_min"]:
            return None
    return nxt


def fmt_sf(sf: dict) -> str:
    return ", ".join(f"{k} {v:.2f}" for k, v in sf.items())


def apply_to_spec(section: str, values: dict, path=CAD_SPEC_PATH) -> None:
    """Rewrite `key: value` lines inside `section:` of spec/cad.yaml, keeping comments and layout."""
    lines = path.read_text().splitlines(keepends=True)
    in_sec = False
    done = set()
    for i, line in enumerate(lines):
        if re.match(r"^\S", line):
            in_sec = line.startswith(f"{section}:")
            continue
        m = re.match(r"^(\s+)(\w+):(\s*)([^#\n]*?)(\s*)(#.*)?$", line.rstrip("\n"))
        if in_sec and m and m.group(2) in values:
            key = m.group(2)
            val = f"{values[key]}"
            comment = m.group(6) or ""
            pad = " " * max(1, len(m.group(4)) + len(m.group(5)) - len(val)) if comment else ""
            lines[i] = f"{m.group(1)}{key}:{m.group(3)}{val}{pad}{comment}\n"
            done.add(key)
    missing = set(values) - done
    if missing:
        raise KeyError(f"{missing} not found in section {section} of {path}")
    path.write_text("".join(lines))


def dump_iterations(studies, path):
    out = {s.name: s.iterations for s in studies}
    path.write_text(json.dumps(out, indent=2))
