"""Redline M9: int8 knock virtual sensor (Xmlacc) and PPO transient-fueling controller. MODELLED."""
from functools import lru_cache

import yaml

from redline.spec import ROOT

SPEC_PATH = ROOT / "spec" / "ml.yaml"
POLICY_DIR = ROOT / "ml" / "policies"
FIXTURE_DIR = ROOT / "ml" / "fixtures"


@lru_cache(maxsize=None)
def ml_spec() -> dict:
    with open(SPEC_PATH) as f:
        return yaml.safe_load(f)
