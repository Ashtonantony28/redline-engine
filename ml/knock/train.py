"""Train, quantize and freeze the int8 knock detector (`python -m ml.knock.train`). MODELLED.

Writes ml/policies/knock_int8.json (int8 weights, scales, float weights for the bound test,
metrics) and ml/fixtures/knock_int8_vectors.json (bit-exact vectors from the golden model).
Both are committed; regenerating them is a reviewed change. Seeded: torch CPU, one thread."""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from ml import FIXTURE_DIR, ml_spec
from ml.knock import golden, metrics, quant, synth

FIXTURE_PATH = FIXTURE_DIR / "knock_int8_vectors.json"
SEED_VAL, SEED_TEST, SEED_FIX = 1, 2, 3          # offsets from knock.seed for the data splits


def datasets(k: dict):
    seed = int(k["seed"])
    out = {}
    for name, n, off in (("train", k["n_train"], 0), ("val", k["n_val"], SEED_VAL), ("test", k["n_test"], SEED_TEST)):
        x, y, meta = synth.make_windows(int(n), seed + off, k)
        out[name] = (x, synth.quantize_input(x, k), y, meta)
    return out


def train_float(k: dict, xq_train, y_train):
    import torch
    torch.manual_seed(int(k["seed"]))
    torch.set_num_threads(1)
    s_in = synth.input_scale(k)
    X = torch.tensor(xq_train.astype(np.float32) * s_in)
    Y = torch.tensor(y_train.astype(np.float32))
    dims = [int(k["window_samples"]), *[int(h) for h in k["hidden"]], 1]
    mods = []
    for i in range(len(dims) - 1):
        mods.append(torch.nn.Linear(dims[i], dims[i + 1]))
        if i < len(dims) - 2:
            mods.append(torch.nn.ReLU())
    net = torch.nn.Sequential(*mods)
    opt = torch.optim.Adam(net.parameters(), lr=float(k["learning_rate"]), weight_decay=float(k["weight_decay"]))
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, int(k["epochs"]))
    loss_fn = torch.nn.BCEWithLogitsLoss()
    g = torch.Generator().manual_seed(int(k["seed"]))
    bs = int(k["batch_size"])
    for _ in range(int(k["epochs"])):
        perm = torch.randperm(len(X), generator=g)
        for i in range(0, len(X), bs):
            idx = perm[i:i + bs]
            opt.zero_grad()
            loss_fn(net(X[idx])[:, 0], Y[idx]).backward()
            opt.step()
        sched.step()
    return [(m.weight.detach().numpy().astype(np.float64), m.bias.detach().numpy().astype(np.float64))
            for m in net if isinstance(m, torch.nn.Linear)]


def evaluate(policy: dict, layers, data, k: dict) -> dict:
    far = float(k["false_alarm_rate"])
    s_in = float(policy["in_scale"])
    _, xq_v, y_v, _ = data["val"]
    _, xq_t, y_t, meta_t = data["test"]
    strong = meta_t["excess_deg"] >= float(k["strong_excess_deg"])
    out = {}
    scores = {"int8": (golden.forward_int8(policy, xq_v)[0], golden.forward_int8(policy, xq_t)[0]),
              "float": (quant.float_forward(layers, xq_v * s_in)[0], quant.float_forward(layers, xq_t * s_in)[0]),
              "bandpass_energy": (synth.bandpass_energy(xq_v * s_in, k), synth.bandpass_energy(xq_t * s_in, k))}
    for name, (sv, st) in scores.items():
        thr = metrics.threshold_at_far(sv[y_v == 0], far)
        det, fa = metrics.rates(st, y_t, thr)
        out[name] = {"auc": metrics.roc_auc(st, y_t), "detection": det, "false_alarm": fa,
                     "detection_strong": float((st[strong] > thr).mean()), "threshold": float(thr)}
    return out


def build(k: dict | None = None):
    k = k or ml_spec()["knock"]
    data = datasets(k)
    layers = train_float(k, data["train"][1], data["train"][2])
    s_in = synth.input_scale(k)
    policy = quant.quantize(layers, s_in, data["train"][1][: int(k["n_val"])] * s_in,
                            int(k["multiplier_bits"]), float(k["activation_percentile"]))
    m = evaluate(policy, layers, data, k)
    policy["threshold_q"] = int(m["int8"]["threshold"])
    policy["float_layers"] = [{"W": W.tolist(), "b": b.tolist()} for W, b in layers]
    policy["metrics"] = m
    policy["source"] = {"generator": "python -m ml.knock.train", "spec": "spec/ml.yaml knock",
                        "seed": int(k["seed"]), "note": "MODELLED: synthetic accelerometer windows",
                        "false_alarm_rate": float(k["false_alarm_rate"])}
    return policy, data


def fixture(policy: dict, k: dict) -> dict:
    n = int(k["fixture_vectors"])
    N = int(k["window_samples"])
    x, y, _ = synth.make_windows(n, int(k["seed"]) + SEED_FIX, k)
    xq = synth.quantize_input(x, k)
    rng = np.random.default_rng(int(k["seed"]) + SEED_FIX)
    edge = np.stack([np.zeros(N), np.full(N, golden.INT8_MAX), np.full(N, golden.INT8_MIN),
                     np.where(np.arange(N) % 2, golden.INT8_MAX, golden.INT8_MIN),
                     rng.integers(golden.INT8_MIN, golden.INT8_MAX + 1, N),
                     rng.integers(golden.INT8_MIN, golden.INT8_MAX + 1, N)]).astype(np.int8)
    tags = ["synth"] * n + ["zero", "max", "min", "alternating", "random", "random"]
    xq = np.concatenate([xq, edge])
    labels = [int(v) for v in y] + [None] * len(edge)
    logit, hidden = golden.forward_int8(policy, xq)
    weights = json.dumps(policy["layers"], sort_keys=True).encode()
    return {"source": {"generator": "python -m ml.knock.train (ml.knock.golden.forward_int8)",
                       "policy_sha256": hashlib.sha256(weights).hexdigest()},
            "vectors": [{"tag": tags[i], "label": labels[i], "x_q": xq[i].astype(int).tolist(),
                         "hidden": [h[i].astype(int).tolist() for h in hidden], "logit": int(logit[i])}
                        for i in range(len(xq))]}


def float_layers(policy: dict):
    return [(np.asarray(L["W"]), np.asarray(L["b"])) for L in policy["float_layers"]]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)
    k = ml_spec()["knock"]
    policy, _ = build(k)
    golden.POLICY_PATH.parent.mkdir(parents=True, exist_ok=True)
    golden.POLICY_PATH.write_text(json.dumps(policy, separators=(",", ":")) + "\n")
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    FIXTURE_PATH.write_text(json.dumps(fixture(policy, k), separators=(",", ":")) + "\n")
    print("knock detector (MODELLED, synthetic windows), test split, threshold at "
          f"{100 * float(k['false_alarm_rate']):.1f} % false alarms on validation:")
    for name, r in policy["metrics"].items():
        print(f"  {name:16s} AUC {r['auc']:.4f}  detection {r['detection']:.3f} (>= {float(k['strong_excess_deg']):.1f} deg: {r['detection_strong']:.3f})  false alarm {r['false_alarm']:.4f}")
    print(f"wrote {golden.POLICY_PATH} and {FIXTURE_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
