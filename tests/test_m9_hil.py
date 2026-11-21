"""M9 exit criteria: the ECU core on the Vortex-V RTL in lockstep with the MVEM plant, the int8
knock sensor on Xmlacc, and the learned transient-fueling controller under its supervisor.
All results MODELLED (MVEM plant, synthetic accelerometer; there is no physical engine).

Verify:  make -C hil cosim && pytest -q tests/test_m9_hil.py tests/test_m9_1_rt.py tests/test_m8_ecu.py && python -m ml.eval
(M9.1: the timing gaps this file reports in frozen time are measured in real time in tests/test_m9_1_rt.py.)

Anchors first, then the exit criteria. Stated tolerances (never loosened silently):
  requantize rule vs Vortex-V's mlacc.rqs model and NeuroFlight's cases .. exact
  portable expf / powf (ecu_libm.c) vs glibc ......................... <= 1 ulp (binary32)
  int8 knock MLP vs its float model ................................. inside the propagated worst-case
                                                                     quantization bound, every window
  RV64 sync-path units .............................................. no soft-float routine referenced
  step-wise session vs hil.runner.run ............................... identical combustion events
  events (tick, type, cyl, arg) and slow-loop outputs, RTL vs host ... bit-exact
  crank-tooth capture at 7000 rpm, 50 MHz core clock ................ <= 50 % of the tooth period
  slow loop ......................................................... <= its 1 ms period
      (two gaps are reported, not hidden: see test_cycle_budget_at_redline)
  int8 sensor on the RTL vs golden fixture .......................... bit-exact, every vector
  knock detector (spec/ml.yaml) ..................................... AUC >= 0.95, detection >= 0.75
                                                                     (>= 0.90 at >= 2 deg) at 1 % FAR
  supervisor ........................................................ per-event lambda inside the safe
                                                                     band in every evaluation cycle
Tests that need riscv64-unknown-elf-gcc / verilator / the Vortex-V checkout skip without them.
"""
from __future__ import annotations

import ctypes as C
import importlib.util
import json
import math
import subprocess

import numpy as np
import pytest

from hil import bridge, runner, rv64
from hil.session import Session
from ml import ml_spec
from ml.knock import export as knock_export
from ml.knock import golden, metrics, quant, synth, train
from mvem import params
from redline.spec import ROOT

ULP_TOL = 1                      # portable libm vs glibc, units in the last place of binary32
N_BOUND_WINDOWS = 400            # windows for the float-vs-int8 bound anchor
N_RQS_RANDOM = 4000
RV64_DIR = ROOT / "ecu" / "targets" / "rv64"

needs_rv = pytest.mark.skipif(bool(rv64.toolchain_missing()) or not rv64.VORTEX_V_ROOT.exists(),
                              reason="needs riscv64-unknown-elf-gcc, verilator and $VORTEX_V_ROOT")
needs_vortex = pytest.mark.skipif(not rv64.VORTEX_V_ROOT.exists(), reason="needs $VORTEX_V_ROOT")


# --------------------------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def knock_cfg():
    return ml_spec()["knock"]


@pytest.fixture(scope="module")
def knock_policy():
    return golden.load_policy()


@pytest.fixture(scope="module")
def knock_fixture():
    return json.loads(train.FIXTURE_PATH.read_text())


@pytest.fixture(scope="module")
def rtl():
    """The firmware on the Vortex-V SoC model (builds model, firmware and host reference)."""
    from hil import cosim
    ecu = rv64.Rv64Ecu()
    runs = {}

    def window(name):
        if name not in runs:
            runs[name] = cosim.run_window(name, ecu)
        return runs[name]
    return ecu, window


@pytest.fixture(scope="module")
def fuel_cfg():
    return ml_spec()["fuel_ppo"]


@pytest.fixture(scope="module")
def fuel_table(fuel_cfg):
    from ml import eval as ml_eval
    return ml_eval.table(fuel_cfg)


def _vortex_module(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, rv64.VORTEX_V_ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _f32_bits(x) -> np.ndarray:
    return np.asarray(x, dtype=np.float32).view(np.int32).astype(np.int64)


# --------------------------------------------------------------------------------------------
# anchors
# --------------------------------------------------------------------------------------------
@needs_vortex
def test_requant_matches_vortex_v_rqs_model():
    """golden.requant(acc, M, S) == mlacc.rqs(acc * M, S) from Vortex-V's own golden model, on random
    operands and on NeuroFlight's rshift_round cases that Vortex-V carries as a fixture."""
    xm = _vortex_module("model/xmlacc_model.py", "vortex_xmlacc_model")
    rng = np.random.default_rng(int(ml_spec()["knock"]["seed"]))
    acc = rng.integers(golden.INT32_MIN, golden.INT32_MAX, N_RQS_RANDOM, endpoint=True)
    m = rng.integers(1, 1 << int(ml_spec()["knock"]["multiplier_bits"]), N_RQS_RANDOM)
    s = rng.integers(golden.SHIFT_MIN, golden.SHIFT_MAX, N_RQS_RANDOM, endpoint=True)
    ours = golden.requant(acc, m, s)
    theirs = [xm.rqs(int(a) * int(mm), int(ss))[0] for a, mm, ss in zip(acc, m, s)]
    assert ours.tolist() == theirs
    cases = json.loads((rv64.VORTEX_V_ROOT / "verif/workloads/neuroflight/fixtures/rshift_cases.json").read_text())
    for v, sh, expect in cases["cases"]:
        if golden.SHIFT_MIN <= sh <= golden.SHIFT_MAX and abs(v) < (1 << 62):
            assert int(golden.requant(v, 1, sh)) == expect, (v, sh)


def test_golden_fixture_is_self_consistent(knock_policy, knock_fixture):
    """The committed vectors are what the committed policy's golden model produces."""
    x = np.array([v["x_q"] for v in knock_fixture["vectors"]], dtype=np.int8)
    logit, hidden = golden.forward_int8(knock_policy, x)
    assert logit.tolist() == [v["logit"] for v in knock_fixture["vectors"]]
    for i, h in enumerate(hidden):
        assert h.tolist() == [v["hidden"][i] for v in knock_fixture["vectors"]]
    assert all(0 <= a <= golden.INT8_MAX for h in hidden for a in h.ravel())
    assert knock_export.main(["--check"]) == 0, "gen/knock_weights.h is stale"


def test_int8_mlp_within_quantization_bound_of_float(knock_policy, knock_cfg):
    """|int8 logit - float logit| <= the worst-case bound propagated through the layers (weight,
    bias and requantize rounding, multiplier error, activation saturation), for every window."""
    layers = train.float_layers(knock_policy)
    x, _, _ = synth.make_windows(N_BOUND_WINDOWS, int(knock_cfg["seed"]) + train.SEED_FIX + 1, knock_cfg)
    xq = synth.quantize_input(x, knock_cfg)
    bound, logit_f = quant.error_bound(layers, knock_policy, xq)
    logit_q = golden.forward_int8(knock_policy, xq)[0] * float(knock_policy["logit_scale"])
    err = np.abs(logit_q - logit_f)
    assert np.all(err <= bound), f"worst excess {np.max(err - bound):.3g}"
    # the bound is not vacuous: it is small against the spread of the logits it protects
    assert np.median(bound) < np.std(logit_f)
    assert np.median(err) < 0.1 * np.std(logit_f)


def test_portable_libm_within_one_ulp_of_glibc():
    """ecu_libm.c (used by the RV64 build and the host reference) against glibc, binary32."""
    subprocess.run(["make", "-s", "-C", str(RV64_DIR), "ref"], check=True)
    lib = C.CDLL(str(rv64.REF_LIB))
    for f, n in ((lib.ecu_expf, 1), (lib.ecu_powf, 2), (lib.ecu_fabsf, 1)):
        f.restype, f.argtypes = C.c_float, [C.c_float] * n
    libm = C.CDLL("libm.so.6")
    libm.expf.restype, libm.expf.argtypes = C.c_float, [C.c_float]
    libm.powf.restype, libm.powf.argtypes = C.c_float, [C.c_float, C.c_float]
    rng = np.random.default_rng(0)
    xs = np.concatenate([rng.uniform(-20.0, 20.0, 3000), rng.uniform(-1e-3, 1e-3, 500), [0.0, -0.0, 88.0, -100.0]])
    worst = 0
    for x in xs.astype(np.float32):
        worst = max(worst, abs(int(_f32_bits(lib.ecu_expf(x)) - _f32_bits(libm.expf(x)))))
    for rpm_ratio, y in zip(rng.uniform(0.05, 12.0, 3000).astype(np.float32), rng.uniform(-2.0, 2.0, 3000).astype(np.float32)):
        worst = max(worst, abs(int(_f32_bits(lib.ecu_powf(rpm_ratio, y)) - _f32_bits(libm.powf(rpm_ratio, y)))))
    assert worst <= ULP_TOL, f"{worst} ulp"
    assert lib.ecu_fabsf(C.c_float(-2.5)) == 2.5 and math.copysign(1.0, lib.ecu_fabsf(C.c_float(-0.0))) == 1.0
    assert lib.ecu_expf(C.c_float(1000.0)) == math.inf and lib.ecu_expf(C.c_float(-1000.0)) == 0.0


@needs_rv
def test_rv64_sync_units_have_no_float_and_core_is_unchanged():
    """rv64im / lp64 turns any float operation into a soft-float call, so the sync-path
    objects (decoder, scheduler, knock control) and the int8 kernel must reference none; the
    slow-loop object must (it is the proof that the check can fail). The portable core sources are
    the committed ones."""
    subprocess.run(["make", "-s", "-C", str(RV64_DIR), "all"], check=True,
                   env={**__import__("os").environ, "VORTEX_V_ROOT": str(rv64.VORTEX_V_ROOT)})

    def undefined(obj):
        out = subprocess.run(["riscv64-unknown-elf-nm", "-u", str(RV64_DIR / "build" / "rv" / obj)],
                             check=True, capture_output=True, text=True).stdout
        return {ln.split()[-1] for ln in out.splitlines() if ln.strip()}

    soft = lambda syms: {s for s in syms if s.startswith("__") and any(t in s for t in ("sf", "df", "tf"))}
    # M9.1 (ADR 0005): the three core units are compiled as one (sync_unit.c includes them unchanged),
    # with rv_div.c for the divisions
    for obj in ("sync_unit.o", "knock_kernel.o", "rv_div.o"):
        assert not soft(undefined(obj)), f"{obj} references soft-float: {soft(undefined(obj))}"
    assert soft(undefined("ecu.o")), "slow loop expected to use soft-float (check is vacuous otherwise)"
    dis = subprocess.run(["riscv64-unknown-elf-objdump", "-d", str(rv64.ELF)], check=True, capture_output=True,
                         text=True).stdout
    # no F/D instruction and no float CSR; the integer CSRs of the timer interrupt (mstatus, mie,
    # mtvec, mcause) are expected since M9.1
    assert not any(tok in dis for tok in ("\tfadd", "\tfmul", "\tfmv", "\tfld", "\tflw", "\tfcvt", "fcsr", "\tfrcsr",
                                          "\tfscsr", "\tfrrm", "\tfsrm", "fflags"))
    if (ROOT / ".git").exists():
        st = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "ecu/src", "ecu/include", "ecu/host"],
                            capture_output=True, text=True).stdout
        assert st.strip() == "", f"portable core / host HAL modified:\n{st}"


def test_session_reproduces_runner():
    """hil.session.Session with zero trim is hil.runner.run step for step (the PPO environment's base)."""
    scn = runner.Scenario.from_cycle("tipin", duration=2.6)
    ref = runner.run(scn, bridge.HostEcu())
    s = Session(scn, bridge.HostEcu())
    while not s.done:
        s.step()
    assert s.eng.events == ref.events
    assert s.ecu.status() == ref.status


# --------------------------------------------------------------------------------------------
# exit criteria: lockstep co-simulation on the Vortex-V RTL
# --------------------------------------------------------------------------------------------
@needs_rv
def test_rtl_event_stream_and_slow_outputs_bit_exact(rtl):
    """Short window (sync, first injections and sparks, start of a tip-in): every drained event
    (tick, type, cyl, arg) and every slow-loop output struct equals the host reference build."""
    _, window = rtl
    r = window("smoke")
    assert r["events"] >= 12 and r["slow_steps"] >= 50 and r["status"]["sync"] == 2
    assert r["diff"]["host_ref"] == {"events": 0, "slow": 0, "first": None}, r["diff"]["host_ref"]["first"]
    assert r["status"]["overflow"] == 0
    # informational: the stock M8 host library uses glibc's expf / powf instead of ecu_libm.c
    d = r["diff"]["host_glibc"]
    print(f"\nvs glibc host build: {d['events']} event and {d['slow']} slow-loop differences in {r['calls']} calls")


@needs_rv
def test_cycle_budget_at_redline(rtl):
    """Cycles per crank-tooth capture and per slow-loop iteration on the RTL, against the tooth
    period at the redline speed for the assumed core clock (hil/cosim.yaml: Vortex-V's Arty target).

    Met: worst crank-tooth capture <= budget_margin of the tooth period; slow loop <= its period.
    Gaps, reported by hil.cosim and asserted to be reported here rather than hidden:
      * slow loop (soft-float) + teeth need more than one core at redline (cpu_utilisation > 1 is
        possible at 50 MHz), and
      * Vortex-V has no interrupts yet, so a tooth arriving during a slow-loop iteration waits for
        it (polled_worst_tooth_latency_us > tooth period)."""
    from hil import cosim
    _, window = rtl
    runs = [window("smoke"), window("redline")]
    assert window("redline")["diff"]["host_ref"]["events"] == 0 and window("redline")["diff"]["host_ref"]["slow"] == 0
    b = cosim.budget(runs)
    print("\n" + cosim.report(runs, b, None))
    cfg = rv64.cosim_spec()
    assert b["tooth_period_us"] == pytest.approx(142.857, abs=0.01)
    assert b["crank_post_to_ack_cycles_max"] >= b["crank_handler_cycles_max"] > 0
    assert b["crank_fraction_of_tooth_period"] <= float(cfg["budget_margin"])
    assert b["slow_fraction_of_period"] <= 1.0
    for key in ("cpu_utilisation_at_redline", "polled_worst_tooth_latency_us", "polled_latency_within_tooth_period"):
        assert key in b
    assert math.isfinite(b["cpu_utilisation_at_redline"])


@needs_rv
def test_int8_knock_sensor_bit_exact_on_rtl(rtl, knock_fixture):
    """Every fixture vector through the Xmlacc kernel (macc / mul / rqs) on the RTL == golden."""
    from hil import cosim
    ecu, _ = rtl
    k = cosim.knock_on_rtl(ecu)
    print(f"\nknock MLP on Xmlacc: {k['cycles_max']} cycles per inference (max), {k['cycles_mean']:.0f} mean")
    assert k["vectors"] == len(knock_fixture["vectors"]) and k["mismatches"] == 0
    tooth_cycles = cosim.tooth_period_s(float(params.engine_spec()["engine"]["redline_rpm"])) * float(
        rv64.cosim_spec()["core_clock_hz"])
    assert k["cycles_max"] < tooth_cycles          # one inference fits between two teeth at redline


@needs_rv
@pytest.mark.slow
def test_rtl_full_tipin_window_bit_exact(rtl):
    """The `make -C hil cosim` window: the M8 tip-in (2000 rpm, pedal 0.10 -> 0.40)."""
    _, window = rtl
    r = window("tipin")
    assert r["events"] > 40
    assert r["diff"]["host_ref"]["events"] == 0 and r["diff"]["host_ref"]["slow"] == 0


# --------------------------------------------------------------------------------------------
# exit criteria: knock detector accuracy
# --------------------------------------------------------------------------------------------
def test_knock_detector_accuracy(knock_policy, knock_cfg):
    """Recomputed from the committed int8 policy on a regenerated test split (seeded)."""
    n = int(knock_cfg["n_test"])
    x, y, meta = synth.make_windows(n, int(knock_cfg["seed"]) + train.SEED_TEST, knock_cfg)
    xq = synth.quantize_input(x, knock_cfg)
    score = golden.forward_int8(knock_policy, xq)[0]
    auc = metrics.roc_auc(score, y)
    det, far = metrics.rates(score, y, int(knock_policy["threshold_q"]))
    strong = meta["excess_deg"] >= float(knock_cfg["strong_excess_deg"])
    m = knock_policy["metrics"]
    assert auc == pytest.approx(m["int8"]["auc"], abs=1e-12)
    assert auc >= float(knock_cfg["min_auc"])
    assert det >= float(knock_cfg["min_detection_at_far"])
    assert float((score[strong] > int(knock_policy["threshold_q"])).mean()) >= float(knock_cfg["min_detection_strong"])
    assert far <= 2.0 * float(knock_cfg["false_alarm_rate"])
    assert m["float"]["auc"] - auc <= float(knock_cfg["max_auc_drop_int8"])
    assert auc > m["bandpass_energy"]["auc"]         # beats the classical band-energy detector


# --------------------------------------------------------------------------------------------
# exit criteria: learned fueling under the supervisor
# --------------------------------------------------------------------------------------------
def test_supervisor_keeps_lambda_in_safe_band_in_eval(fuel_table, fuel_cfg):
    """PPO + supervisor on every evaluation cycle: every combustion event from the steady window on
    is inside the safe band, and the applied trim never exceeds its authority."""
    lo, hi = (float(v) for v in fuel_cfg["safe_band"])
    for r in fuel_table:
        p = r["ppo"]
        assert lo <= p["lambda_min"] and p["lambda_max"] <= hi, r
        assert p["max_abs_trim"] <= float(fuel_cfg["trim_limit"]) + 1e-12


def test_supervisor_contains_adversarial_policies(fuel_cfg):
    """The safety claim does not rest on the policy: constant, bang-bang and transient-seeking
    adversaries (lean on tip-in, rich on tip-out) stay inside the safe band too."""
    from ml import eval as ml_eval
    lo, hi = (float(v) for v in fuel_cfg["safe_band"])
    count = [0]

    def bang(_obs):
        count[0] += 1
        return 1.0 if (count[0] // 25) % 2 else -1.0

    def worst(obs):                       # obs[1] = pedal rate, obs[3] = MAP rate (ml.fuel.env.OBS_NAMES)
        return -1.0 if obs[1] > 0.0 or obs[3] > 0.0 else 1.0

    for policy in (lambda _o: 1.0, lambda _o: -1.0, bang, worst):
        for cycle in fuel_cfg["eval_cycles"]:
            for rpm in (fuel_cfg["eval_rpm"][0], fuel_cfg["eval_rpm"][-1]):
                r = ml_eval.run_cycle(cycle, float(rpm), policy, True, fuel_cfg)
                assert lo <= r["lambda_min"] and r["lambda_max"] <= hi, (cycle, rpm, r)


def test_ppo_does_not_worsen_classical_excursion(fuel_table):
    """Comparison table sanity: averaged over the evaluation cycles the supervised PPO controller's
    rms lambda error in the transient window is no worse than the classical ECU's."""
    mean = lambda who, key: float(np.mean([r[who][key] for r in fuel_table]))
    print(f"\nrms: classical {mean('classical', 'rms_dev'):.4f} -> PPO+supervisor {mean('ppo', 'rms_dev'):.4f}; "
          f"max: classical {mean('classical', 'max_dev'):.4f} -> PPO+supervisor {mean('ppo', 'max_dev'):.4f}")
    assert mean("ppo", "rms_dev") <= mean("classical", "rms_dev")
