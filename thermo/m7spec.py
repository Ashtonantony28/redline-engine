"""Load spec/combustion_ext.yaml (M7: knock, flame, emissions, spark, maps). spec/engine.yaml stays the
single source of truth; this file only adds model constants (CLAUDE.md rule 3)."""
from __future__ import annotations

import yaml

from redline.spec import ROOT

COMBUSTION_SPEC = ROOT / "spec" / "combustion_ext.yaml"


def combustion_spec(path: str | None = None) -> dict:
    with open(path or COMBUSTION_SPEC) as f:
        return yaml.safe_load(f)
