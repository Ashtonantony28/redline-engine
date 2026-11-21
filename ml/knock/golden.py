"""numpy golden model of the int8 knock detector: the bit-exact reference for the C / Xmlacc kernel
(ecu/targets/rv64/knock_kernel.c) on the Vortex-V RTL.

Contract (the integer part of NeuroFlight's int8 contract, with ReLU in place of the tanh LUT), per
output j of a layer:

    acc = sat_int32(b_q[j] + sum_i W_q[j, i] * x_q[i])
    hidden:  y = clamp(rshift_round(acc * M[j], S[j]), -128, 127);   a = max(y, 0)
    last:    acc is the output (one unit: the knock logit in units of logit_scale)

    rshift_round(v, s) = (v + 2^(s-1)) >> s        arithmetic shift, round half up, 1 <= s <= 63

which is exactly Vortex-V's mlacc.rqs (docs/XMLACC-SPEC-v1.md) and NeuroFlight's requantize rule."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ml import POLICY_DIR

INT8_MIN, INT8_MAX = -128, 127
INT32_MIN, INT32_MAX = -(1 << 31), (1 << 31) - 1
SHIFT_MIN, SHIFT_MAX = 1, 63
POLICY_PATH = POLICY_DIR / "knock_int8.json"


def rshift_round(v, s):
    """(v + 2^(s-1)) >> s on int64 arrays (s >= 1); the caller guarantees no int64 overflow."""
    v = np.asarray(v, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    return (v + (np.int64(1) << (s - 1))) >> s


def requant(acc, m, s):
    """clamp(rshift_round(acc * M, S), -128, 127)."""
    return np.clip(rshift_round(np.asarray(acc, dtype=np.int64) * np.asarray(m, dtype=np.int64), s),
                   INT8_MIN, INT8_MAX)


def load_policy(path: Path | str = POLICY_PATH) -> dict:
    return json.loads(Path(path).read_text())


def forward_int8(policy: dict, x_q: np.ndarray):
    """x_q: int8 [n, n_in] (or [n_in]). Returns (logit int32 [n], [hidden int8 [n, h] ...])."""
    a = np.atleast_2d(np.asarray(x_q)).astype(np.int64)
    hidden = []
    for L in policy["layers"]:
        W = np.asarray(L["W_q"], dtype=np.int64)
        acc = np.clip(a @ W.T + np.asarray(L["b_q"], dtype=np.int64), INT32_MIN, INT32_MAX)
        if "M" not in L:
            return acc[:, 0].astype(np.int32), hidden
        a = np.maximum(requant(acc, L["M"], L["S"]), 0)
        hidden.append(a.astype(np.int8))
    raise ValueError("policy has no output layer")


def detect(policy: dict, x_q: np.ndarray) -> np.ndarray:
    return forward_int8(policy, x_q)[0] > int(policy["threshold_q"])
