"""M9.1 exit criteria: the ECU on Vortex-V in REAL time (plant time = core cycles, ADR 0005).
All results MODELLED (MVEM plant, Verilator model of the Vortex-V SoC; there is no physical engine).

Verify:  make -C hil cosim && pytest -q tests/test_m9_hil.py tests/test_m9_1_rt.py tests/test_m8_ecu.py

Anchors first, then the exit criteria. Stated tolerances (never loosened silently):
  rv_softfp.c / rv_div.c vs the host FPU and integer divide .......... bit-identical, every vector
                                                                      (special-value pairs + seeded random)
  frozen-time mode, events and slow-loop outputs, RTL vs host ........ bit-exact (unchanged from M9)
  slow loop on the RTL (frozen-time handler cost) .................... <= 50 % of its 1 ms period at 50 MHz
  real time, 7000 rpm and 2000 rpm tip-in, 50 MHz (the board clock):
      captures dropped / loss of sync / dropouts / slow-loop overruns  0
      CPU load ...................................................... <= 60 %
      spark angle error vs the host build in ideal time ............. <= 0.5 deg (the M8 scheduler spec)
      every event, and the worst tooth-to-event latency as an angle .. <= 0.5 deg at 2000 rpm
  real time, 7000 rpm, 100 MHz (the board oscillator undivided) ...... all of the above, every event and
                                                                      the worst-case latency <= 0.5 deg
  real time, 7000 rpm, 50 MHz, since M9.3 (ADR 0008 scheduler arm-ahead): every event and the worst
  case (worst tooth-to-event latency minus the scheduler's lead) <= 0.5 deg. Until M9.3 this was a
  KNOWN GAP (0.70 deg observed, 0.75 deg worst case): an event right behind its reference tooth was
  armed late.
  the M9 firmware (polled, libgcc soft-float) in real time ........... fails: spark error > 5 deg (slow-marked)
Tests that need riscv64-unknown-elf-gcc / verilator / the Vortex-V checkout skip without them.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from hil import rv64
from redline.spec import ROOT

RV64_DIR = ROOT / "ecu" / "targets" / "rv64"
ADR = ROOT / "docs" / "adr" / "0005-rv64-real-time-hal.md"
SPEC_DEG = float(rv64.cosim_spec()["realtime"]["angle_error_max_deg"])
BOARD_HZ = float(rv64.cosim_spec()["core_clock_hz"])
OSC_HZ = 100.0e6                 # Arty A7 oscillator (Vortex-V fpga/arty_a7_100t.xdc), MMCM divider 10
LOAD_MAX = 0.60
SLOW_FRACTION_MAX = 0.50
M9_SPARK_ERROR_MIN_DEG = 5.0     # the M9 firmware is at least this wrong in real time
START_INJ_DIFF_MAX = 4           # one start-up injection per cylinder
SOFTFP_VECTORS = 300000

needs_rv = pytest.mark.skipif(bool(rv64.toolchain_missing()) or not rv64.VORTEX_V_ROOT.exists(),
                              reason="needs riscv64-unknown-elf-gcc, verilator and $VORTEX_V_ROOT")


@pytest.fixture(scope="module")
def bench():
    from hil import realtime
    return realtime.Bench()


@pytest.fixture(scope="module")
def redline_50(bench):
    return bench.run("redline_rt", "irq", BOARD_HZ)


# --------------------------------------------------------------------------------------------
# anchors
# --------------------------------------------------------------------------------------------
def test_adr_0005_present_and_accepted():
    text = ADR.read_text()
    assert text.startswith("# ADR 0005")
    assert "Status: Accepted" in text
    for word in ("MTIMECMP", "snapshot", "frozen", "bit-exact", "Consequences"):
        assert word in text, word


def test_softfp_and_divide_bit_identical_to_host():
    """The target's soft-float and reciprocal divide, compiled natively, against the host's FPU and
    integer divide: add / sub / mul / div in binary32 and binary64, comparisons, conversions,
    64-bit divide, on special-value pairs and seeded random operands."""
    out = subprocess.run(["make", "-s", "-C", str(RV64_DIR), "softfp-check", f"SOFTFP_VECTORS={SOFTFP_VECTORS}"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "PASS" in out.stdout
    rows = [ln.split() for ln in out.stdout.splitlines() if "mismatches" in ln]
    assert len(rows) == 11 and all(int(r[1]) > SOFTFP_VECTORS and int(r[3]) == 0 for r in rows)


def test_portable_core_changed_only_by_adr_0008():
    """ADR 0003 / 0005: the RV64 work is HAL, build and target arithmetic only. The core itself
    changed once since, on purpose, in M9.3 (ADR 0008: scheduler arm-ahead): relative to the M9.2
    commit only the ADR 0008 files may differ. (Until M9.3 this required a clean `git status`.)"""
    changed = rv64.core_files_changed_since(rv64.CORE_BASE_M92)
    if changed is None:
        pytest.skip("not a git checkout, or the M9.2 commit is unknown")
    assert changed <= rv64.ADR_0008_CORE_FILES, f"core changed outside ADR 0008: {sorted(changed - rv64.ADR_0008_CORE_FILES)}"


@needs_rv
def test_fast_build_links_no_libgcc_arithmetic_and_no_sync_divide():
    """The M9.1 build takes its soft-float from rv_softfp.c and has no `div` in the sync path."""
    subprocess.run(["make", "-s", "-C", str(RV64_DIR), "all"], check=True,
                   env={**os.environ, "VORTEX_V_ROOT": str(rv64.VORTEX_V_ROOT)})
    syms = rv64.elf_symbols()
    assert "__clz_tab" not in syms and "__clzdi2" not in syms          # libgcc soft-fp drags these in
    for name in ("__divdf3", "__mulsf3", "__udivdi3"):
        assert name in syms
    dis = subprocess.run(["riscv64-unknown-elf-objdump", "-d", str(RV64_DIR / "build" / "rv" / "sync_unit.o")],
                         check=True, capture_output=True, text=True).stdout
    assert not any(f"\t{op}" in dis for op in ("div\t", "divu\t", "rem\t", "remu\t", "divw\t", "divuw\t", "remw\t", "remuw\t"))


# --------------------------------------------------------------------------------------------
# frozen time: bit-exactness is kept
# --------------------------------------------------------------------------------------------
@needs_rv
def test_frozen_time_still_bit_exact_and_slow_loop_fits(bench):
    """The frozen-time lockstep mode on the M9.1 firmware: events and slow-loop outputs bit-exact
    against the host reference build, and the slow loop now fits half its period."""
    from hil import cosim
    ecu = bench.ecu("fast")
    runs = [cosim.run_window(w, ecu, with_glibc=False) for w in ("smoke", "redline")]
    for r in runs:
        assert r["events"] > 0 and r["slow_steps"] > 0
        assert r["diff"]["host_ref"] == {"events": 0, "slow": 0, "first": None}
    b = cosim.budget(runs)
    print(f"\nfrozen time, M9.1 build: crank handler mean {b['crank_handler_cycles_mean']:.0f} / max "
          f"{b['crank_handler_cycles_max']} cycles, slow loop mean {b['slow_handler_cycles_mean']:.0f} / max "
          f"{b['slow_handler_cycles_max']} cycles (MODELLED)")
    assert b["slow_fraction_of_period"] <= SLOW_FRACTION_MAX
    assert b["cpu_utilisation_at_redline"] < 1.0


# --------------------------------------------------------------------------------------------
# real time
# --------------------------------------------------------------------------------------------
def _held(m):
    assert m["captures_dropped"] == 0 and m["captures_unserviced"] <= 1
    assert m["n_loss"] == 0 and m["n_dropout"] == 0 and m["n_noise"] == 0 and m["n_cam_err"] == 0
    assert m["sync"] == 2 and m["n_sync"] == m["host_n_sync"] == 1
    assert m["slow_overruns"] == 0 and abs(m["slow_steps"] - m["slow_steps_expected"]) <= 1
    assert m["spark_count_difference"] == 0 and m["sparks_compared"] > 0
    # the first injection of a cylinder can fall one cycle earlier or later than on the host: after
    # sync the slow loop's first pulse widths reach the scheduler when its step is committed
    assert abs(m["events"] - m["host_events"]) <= START_INJ_DIFF_MAX


@needs_rv
def test_real_time_redline_at_board_clock(redline_50):
    """7000 rpm, 50 MHz, time advancing: no missed tooth, no loss of sync, the slow loop keeps its
    1 ms period, the CPU has headroom, and every spark is within 0.5 deg of the host's."""
    from hil import realtime
    r, m = redline_50
    print("\n" + realtime.line(m))
    _held(m)
    assert m["rpm"] == pytest.approx(7000.0)
    assert m["cpu_load"] <= LOAD_MAX
    assert m["spark_angle_error_max_deg"] <= SPEC_DEG
    assert m["tooth_to_serviced_max_us"] < m["tooth_period_us"]         # every tooth done before the next
    assert m["irqs"] > 0 and m["spin_timeouts"] <= 5                    # only while the decoder finds the gap


@needs_rv
def test_real_time_redline_meets_spec_at_board_clock(redline_50):
    """M9.3 (ADR 0008): the scheduler arms an event that lies less than the arm-ahead guard behind a
    tooth from the tooth before, so the handler's latency (about 19 us, 0.8 deg at 7000 rpm) no longer
    makes it late: every event and the worst case are within 0.5 deg at 50 MHz in the M9.1 mailbox
    mode too. (Until M9.3 this test asserted the known gap: worst case above 0.5 deg, below 1 deg.)"""
    _, m = redline_50
    assert m["events_late"] == 0
    assert m["lead_min_deg"] > m["tooth_to_event_max_deg"]
    assert m["event_angle_error_max_deg"] <= SPEC_DEG and m["worst_case_angle_error_deg"] <= SPEC_DEG
    assert m["meets_spec"] and m["rt_ok"]


@needs_rv
def test_real_time_redline_meets_spec_at_100mhz(bench):
    """With the board oscillator undivided the whole spec holds at 7000 rpm: every event and the
    worst tooth-to-event latency are within 0.5 deg of crank travel."""
    from hil import realtime
    _, m = bench.run("redline_rt", "irq", OSC_HZ)
    print("\n" + realtime.line(m))
    _held(m)
    assert m["event_angle_error_max_deg"] <= SPEC_DEG
    assert m["worst_case_angle_error_deg"] <= SPEC_DEG
    assert m["meets_spec"]


@needs_rv
@pytest.mark.slow
def test_real_time_tipin_at_board_clock(bench):
    """2000 rpm tip-in, 50 MHz: the whole spec holds (0.5 deg is 42 us there)."""
    from hil import realtime
    _, m = bench.run("tipin", "irq", BOARD_HZ)
    print("\n" + realtime.line(m))
    _held(m)
    assert m["cpu_load"] <= LOAD_MAX
    assert m["event_angle_error_max_deg"] <= SPEC_DEG and m["worst_case_angle_error_deg"] <= SPEC_DEG
    assert m["meets_spec"]


@needs_rv
@pytest.mark.slow
def test_m9_firmware_fails_in_real_time(bench):
    """The time model is honest: the M9 structure (polled captures, libgcc soft-float, hardware
    divide) passes the frozen-time comparison but, with time advancing, fires sparks degrees late
    and cannot keep the slow loop's period at 7000 rpm."""
    from hil import realtime
    _, m = bench.run("redline_rt", "m9", BOARD_HZ)
    print("\n" + realtime.line(m))
    assert m["spark_angle_error_max_deg"] > M9_SPARK_ERROR_MIN_DEG
    assert m["tooth_to_serviced_max_us"] > m["tooth_period_us"]
    assert m["slow_overruns"] > 0
    assert not m["meets_spec"]
