"""Post-training int8 quantization of the float knock MLP, and the float-vs-int8 error bound.

Per layer, with input scale s_x (real units per int8 code) and one weight scale per output row:
    W_q[j] = round(W[j] / w_j),  w_j = max|W[j]| / 127        b_q[j] = round(b[j] / (s_x w_j))
    hidden: M[j] 2^-S[j] ~= s_x w_j / s_a,  M normalised to `multiplier_bits` bits
            s_a = (calibration percentile of the float activations) / 127
The last layer keeps its int32 accumulator; logit = acc * s_x * w_0."""
from __future__ import annotations

import math

import numpy as np

from ml.knock import golden

INT8_MAX = golden.INT8_MAX


def float_forward(layers, x):
    """Float reference MLP (float64): returns (logit [n], [post-ReLU activations ...])."""
    a = np.atleast_2d(np.asarray(x, dtype=np.float64))
    acts = []
    for i, (W, b) in enumerate(layers):
        z = a @ np.asarray(W, dtype=np.float64).T + np.asarray(b, dtype=np.float64)
        if i == len(layers) - 1:
            return z[:, 0], acts
        a = np.maximum(z, 0.0)
        acts.append(a)
    raise ValueError("empty model")


def multiplier(r: float, bits: int):
    """(M, S) with M 2^-S ~= r, 2^(bits-1) <= M < 2^bits."""
    if not r > 0.0:
        raise ValueError("requantize ratio must be positive")
    S = bits - 1 - math.floor(math.log2(r))
    M = int(round(r * 2.0 ** S))
    if M >= 1 << bits:
        M >>= 1
        S -= 1
    if not golden.SHIFT_MIN <= S <= golden.SHIFT_MAX:
        raise ValueError(f"shift {S} outside mlacc.rqs range")
    return M, S


def quantize(layers, s_in: float, x_cal: np.ndarray, bits: int, percentile: float) -> dict:
    """layers: [(W, b), ...] float; x_cal: calibration inputs in real units (already on the int8 grid)."""
    _, acts = float_forward(layers, x_cal)
    out = []
    s_x = s_in
    for i, (W, b) in enumerate(layers):
        W = np.asarray(W, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64)
        w = np.abs(W).max(axis=1) / INT8_MAX
        w[w == 0.0] = 1.0
        d = {"W_q": np.rint(W / w[:, None]).astype(int).tolist(),
             "b_q": np.clip(np.rint(b / (s_x * w)), golden.INT32_MIN, golden.INT32_MAX).astype(int).tolist(),
             "w_scale": w.tolist(), "in_scale": s_x}
        if i < len(layers) - 1:
            s_a = float(np.percentile(acts[i], percentile)) / INT8_MAX
            ms = [multiplier(s_x * wj / s_a, bits) for wj in w]
            d.update(M=[m for m, _ in ms], S=[s for _, s in ms], out_scale=s_a)
            s_x = s_a
        out.append(d)
    return {"layers": out, "in_scale": s_in, "logit_scale": out[-1]["in_scale"] * out[-1]["w_scale"][0]}


def error_bound(layers, policy: dict, x_q: np.ndarray):
    """Worst-case |int8 logit (real units) - float logit| per sample, for the float model fed the
    same int8-grid input. Sources, propagated through the layers:
      weight rounding  <= w_j / 2 per weight      bias rounding  <= s_x w_j / 2
      requantize       <= s_a / 2 (round half up) + |z| |M 2^-S s_a / (s_x w_j) - 1| (multiplier)
      saturation       float activation above 127 s_a is clipped (added explicitly)
    Returns (bound [n], float logit [n])."""
    x_q = np.atleast_2d(x_q).astype(np.int64)
    a_f = x_q * float(policy["in_scale"])
    a_q = x_q.copy()
    e_a = np.zeros_like(a_f)
    for (W, b), L in zip(layers, policy["layers"]):
        W = np.asarray(W, dtype=np.float64)
        w = np.asarray(L["w_scale"])
        s_x = float(L["in_scale"])
        z_f = a_f @ W.T + np.asarray(b, dtype=np.float64)
        acc = np.clip(a_q @ np.asarray(L["W_q"], dtype=np.int64).T + np.asarray(L["b_q"], dtype=np.int64),
                      golden.INT32_MIN, golden.INT32_MAX)
        e_z = e_a @ np.abs(W).T + (np.abs(a_q).sum(axis=1, keepdims=True) * s_x) * (w / 2.0) + s_x * w / 2.0
        if "M" not in L:
            return e_z[:, 0], z_f[:, 0]
        s_a = float(L["out_scale"])
        z_q = acc * (s_x * w)
        delta = np.abs(np.asarray(L["M"]) * 2.0 ** -np.asarray(L["S"], dtype=float) * s_a / (s_x * w) - 1.0)
        a_f = np.maximum(z_f, 0.0)
        e_a = e_z + s_a / 2.0 + delta * np.abs(z_q) + np.maximum(a_f - INT8_MAX * s_a, 0.0)
        a_q = np.maximum(golden.requant(acc, L["M"], L["S"]), 0)
    raise ValueError("policy has no output layer")
