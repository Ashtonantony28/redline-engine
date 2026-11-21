"""Synthetic block-accelerometer windows for the knock detector (spec/ml.yaml knock). MODELLED.

One window per combustion event, opening at CA50. Knock is a damped burst at the cylinder's first
circumferential acoustic mode (Draper: f = 1.841 c / (pi B), plus a weaker second mode) starting a few
degrees after CA50; its amplitude is proportional to the spark advance beyond the knock limit, the
same law the MVEM knock model uses for the window integral. The background is broadband noise that
grows with the reciprocating inertia force (rpm^2), the low-frequency block response to the
pressure rise, and valve-seating impacts. Everything is drawn from one seeded generator."""
from __future__ import annotations

import math

import numpy as np

from ml import ml_spec
from mvem import params

INT8_MAX = 127
INT8_MIN = -128
DEG_PER_REV = 360.0


def first_mode_hz(t_burned: float, k: dict | None = None) -> float:
    k = k or ml_spec()["knock"]
    bore = params.geometry().bore
    c = math.sqrt(float(k["burned_gamma"]) * float(k["burned_gas_constant"]) * t_burned)
    return float(k["mode_constant_1"]) * c / (math.pi * bore)


def input_scale(k: dict | None = None) -> float:
    """Acceleration [g] per int8 code."""
    k = k or ml_spec()["knock"]
    return float(k["adc_full_scale_g"]) / INT8_MAX


def quantize_input(x_g: np.ndarray, k: dict | None = None) -> np.ndarray:
    return np.clip(np.rint(x_g / input_scale(k)), INT8_MIN, INT8_MAX).astype(np.int8)


def _burst(t, t0, amp, f, tau, phase):
    dt = t - t0
    on = dt >= 0.0
    return np.where(on, amp * np.exp(-np.where(on, dt, 0.0) / tau) * np.sin(2.0 * np.pi * f * dt + phase), 0.0)


def make_windows(n: int, seed: int, k: dict | None = None):
    """Returns (x_g [n, N] float64, label [n] 0/1, meta dict of per-window arrays)."""
    k = k or ml_spec()["knock"]
    rng = np.random.default_rng(seed)
    es = params.engine_spec()
    redline = float(es["engine"]["redline_rpm"])
    ncyl = params.geometry().cylinders
    fs = float(k["sample_rate_hz"])
    N = int(k["window_samples"])
    t = np.arange(N) / fs
    u = lambda lohi, size=None: rng.uniform(float(lohi[0]), float(lohi[1]), size)
    x = np.zeros((n, N))
    label = (rng.random(n) < float(k["knock_probability"])).astype(np.int64)
    meta = {key: np.zeros(n) for key in ("rpm", "map", "excess_deg", "knock_g", "f1_hz", "cyl")}
    ratio = float(k["mode_constant_2"]) / float(k["mode_constant_1"])
    for i in range(n):
        rpm = u(k["rpm_range"])
        p_map = u(k["map_range_pa"])
        load = p_map / float(k["map_range_pa"][1])
        cyl = int(rng.integers(ncyl))
        path = float(k["cylinder_path_gain"][cyl])
        sec_per_deg = params.SEC_PER_MIN / (rpm * DEG_PER_REV)
        sigma = (float(k["noise_rms_idle_g"]) + float(k["noise_rms_inertia_g"]) * (rpm / redline) ** 2)
        sigma *= 1.0 + float(k["noise_gain_spread"]) * rng.uniform(-1.0, 1.0)
        w = sigma * rng.standard_normal(N)
        w += _burst(t, 0.0, u(k["combustion_thump_g"]) * load * path, u(k["combustion_thump_hz"]),
                    float(k["combustion_thump_decay_s"]), rng.uniform(0.0, 2.0 * np.pi))
        if rng.random() < float(k["valve_impact_probability"]):
            w += _burst(t, rng.uniform(0.0, t[-1]), u(k["valve_impact_g"]), u(k["valve_impact_hz"]),
                        u(k["valve_impact_decay_s"]), rng.uniform(0.0, 2.0 * np.pi))
        f1 = first_mode_hz(u(k["burned_temperature_k"]), k)
        if label[i]:
            excess = u(k["spark_excess_deg"])
            amp = float(k["knock_g_per_deg"]) * excess * (0.5 + 0.5 * load) * path
            t0 = u(k["knock_onset_after_ca50_deg"]) * sec_per_deg
            tau = u(k["knock_decay_s"])
            w += _burst(t, t0, amp, f1, tau, rng.uniform(0.0, 2.0 * np.pi))
            w += _burst(t, t0, amp * u(k["mode_2_amplitude_ratio"]), f1 * ratio, tau, rng.uniform(0.0, 2.0 * np.pi))
            meta["excess_deg"][i], meta["knock_g"][i] = excess, amp
        x[i] = w
        meta["rpm"][i], meta["map"][i], meta["f1_hz"][i], meta["cyl"][i] = rpm, p_map, f1, cyl
    return x, label, meta


def bandpass_energy(x_g: np.ndarray, k: dict | None = None) -> np.ndarray:
    """Classical reference detector: spectral energy in the knock band of each window."""
    k = k or ml_spec()["knock"]
    fs = float(k["sample_rate_hz"])
    N = x_g.shape[-1]
    spec = np.abs(np.fft.rfft(x_g * np.hanning(N), axis=-1)) ** 2
    f = np.fft.rfftfreq(N, 1.0 / fs)
    lo, hi = (float(v) for v in k["bandpass_hz"])
    return spec[..., (f >= lo) & (f <= hi)].sum(axis=-1)
