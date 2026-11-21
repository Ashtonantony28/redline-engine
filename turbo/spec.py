"""Load spec/turbo.yaml and build the Redline I4-2.0T override layer (spec/engine.yaml stays untouched)."""
from __future__ import annotations

import copy

import yaml

from gasdyn.spec import gasdyn_spec, with_overrides
from redline.spec import ROOT, load

TURBO_SPEC = ROOT / "spec" / "turbo.yaml"


def turbo_spec(path: str | None = None) -> dict:
    with open(path or TURBO_SPEC) as f:
        return yaml.safe_load(f)


def engine_spec(tspec: dict | None = None) -> dict:
    """spec/engine.yaml with the variant's dotted-key overrides (deep copy: load() is cached)."""
    tspec = tspec or turbo_spec()
    return with_overrides(copy.deepcopy(load()), tspec["variant"]["engine_overrides"])


def gasdyn_variant_spec(tspec: dict | None = None) -> dict:
    tspec = tspec or turbo_spec()
    return with_overrides(gasdyn_spec(), tspec["variant"].get("gasdyn_overrides") or {})
