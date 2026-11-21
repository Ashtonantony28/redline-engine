"""Load spec/gasdyn.yaml (M4 parameters; spec/engine.yaml stays the single source of truth)."""
from __future__ import annotations

import copy

import yaml

from redline.spec import ROOT

GASDYN_SPEC = ROOT / "spec" / "gasdyn.yaml"
VALVETRAIN_SPEC = ROOT / "spec" / "valvetrain.yaml"


def gasdyn_spec(path: str | None = None) -> dict:
    with open(path or GASDYN_SPEC) as f:
        return yaml.safe_load(f)


def valvetrain_spec(path: str | None = None) -> dict:
    with open(path or VALVETRAIN_SPEC) as f:
        return yaml.safe_load(f)


def with_overrides(spec: dict, overrides: dict | None) -> dict:
    """Deep copy of `spec` with dotted-key overrides, e.g. {"intake.runner.length": 0.5}."""
    out = copy.deepcopy(spec)
    for key, value in (overrides or {}).items():
        node = out
        parts = key.split(".")
        for part in parts[:-1]:
            node = node[part]
        if parts[-1] not in node:
            raise KeyError(key)
        node[parts[-1]] = value
    return out
