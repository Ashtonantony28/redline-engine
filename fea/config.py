"""Spec access for the M6 FEA package.

spec/engine.yaml stays the single source of truth (geometry, masses, materials); spec/cad.yaml holds the
detailed geometry the FEA reconstructs; spec/fea.yaml holds FEA-only parameters (mesh sizes, boundary
condition assumptions, targets, material supplements). ADR 0001 units throughout.
"""
from __future__ import annotations

import copy
import hashlib
import json

import numpy as np

from cad.params import CAD_SPEC_PATH, Params, params
from redline.spec import ROOT, SPEC_PATH, load

FEA_SPEC_PATH = ROOT / "spec" / "fea.yaml"
SPEC_FILES = ("engine.yaml", "cad.yaml", "fea.yaml", "thermo.yaml", "dynamics.yaml")


def fea_spec(path: str | None = None) -> dict:
    return load(str(path or FEA_SPEC_PATH))


def base_params() -> Params:
    return params()


def with_overrides(p: Params, section: str, overrides: dict) -> Params:
    """New Params with spec/cad.yaml values replaced in memory (used by the iteration loop)."""
    cad = copy.deepcopy(p.cad)
    cad[section].update(overrides)
    return Params(p.engine, cad, p.k)


def spec_hash() -> str:
    """Hash of every spec file the FEA results depend on (stale-result detection)."""
    h = hashlib.sha256()
    for name in SPEC_FILES:
        h.update((ROOT / "spec" / name).read_bytes())
    return h.hexdigest()[:16]


def material(engine: dict, key: str) -> dict:
    """engine.yaml material with floats forced (PyYAML may hand back strings for bare exponents);
    temperature tables (lists) become float arrays."""
    return {k: (np.asarray(v, dtype=float) if isinstance(v, list) else float(v))
            for k, v in engine["materials"][key].items()}


def strength_at(T, m: dict):
    """(S_ut, S_e) [Pa] of material m (from `material`) at temperature(s) T [K]: linear interpolation of its
    cited engine.yaml tables T_K / sut_T / se_T, clamped at the ends. A material without tables has no
    temperature data and is refused (no assumed derating)."""
    if "T_K" not in m:
        raise KeyError("material has no temperature tables (T_K, sut_T, se_T) in spec/engine.yaml")
    T = np.asarray(T, dtype=float)
    return np.interp(T, m["T_K"], m["sut_T"]), np.interp(T, m["T_K"], m["se_T"])


def dumps(obj) -> str:
    def conv(o):
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    return json.dumps(obj, indent=2, default=conv)


__all__ = ["CAD_SPEC_PATH", "SPEC_PATH", "fea_spec", "base_params", "with_overrides", "spec_hash",
           "material", "strength_at", "dumps"]
