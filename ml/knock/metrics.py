"""Detector metrics without extra dependencies: ROC AUC (rank statistic) and detection rate at a
fixed false-alarm rate."""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def roc_auc(score, label) -> float:
    score = np.asarray(score, dtype=np.float64)
    label = np.asarray(label).astype(bool)
    n_pos, n_neg = int(label.sum()), int((~label).sum())
    r = rankdata(score)                         # average ranks: ties count one half
    return float((r[label].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def threshold_at_far(score_neg, far: float):
    """Smallest threshold t such that P(score > t | no knock) <= far on the given negatives."""
    s = np.sort(np.asarray(score_neg))
    k = int(np.ceil((1.0 - far) * len(s))) - 1
    return s[min(max(k, 0), len(s) - 1)]


def rates(score, label, threshold):
    score = np.asarray(score)
    label = np.asarray(label).astype(bool)
    hit = score > threshold
    return float(hit[label].mean()), float(hit[~label].mean())      # detection, false alarm
