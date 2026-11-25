"""M9.2 exit criteria: the ECU on Vortex-V with the SoC's own capture/compare timer (ADR 0007;
Vortex-V ADR 0017, CAPCMP=1). All results MODELLED (MVEM plant, Verilator model of the Vortex-V
SoC; there is no physical engine).

Verify:  make -C hil cosim && pytest -q tests/test_m9_hil.py tests/test_m9_1_rt.py tests/test_m9_2_capcmp.py tests/test_m8_ecu.py

Anchors first, then the exit criteria. Stated tolerances (never loosened silently):
  portable core (ecu/src, ecu/include, ecu/host) .................... unchanged
  frozen-time mode on the CAPCMP=1 SoC, events and slow-loop outputs . bit-exact vs the host build
  M9 and M9.1 firmware builds ........................................ still build (comparison modes)
  real time, capture/compare mode, 7000 rpm, 50 MHz (the board clock):
      every crank / cam edge in the RTL FIFO with the host's tick ..... all, CAP_OVF 0
      loss of sync / dropouts / noise / slow-loop overruns ........... 0
      compare slots: pending ........................................ <= 8 at every cycle, 0 pushes
                                                                     dropped by the unit, 0 held back
      every armed compare accounted for on cmp_out .................. no missing edge, every
                                                                     injector pulse exactly `arg` ticks
      an event armed before its tick ................................ its pin edge lands on that tick
      spark angle error vs the host build in ideal time ............. <= 0.5 deg
      CPU load ...................................................... <= 60 % (the M9.1 bound)
  KNOWN GAP, asserted so that it is reported and not hidden: the unit fires on time whatever is
  armed on time, but the portable scheduler still arms an event from the tooth just before it, so an
  event inside the handler's latency behind a tooth (the injection "late catch-up" starts, right at
  the tooth) is pushed after its tick and fires one clock later. At 7000 rpm and 50 MHz that is more
  than 0.5 deg (and less than 1 deg) for the worst event and for the worst-case latency bound: the
  "every event within 0.5 deg at 50 MHz" criterion is NOT met, see ADR 0007 for the fix (scheduler
  arm-ahead, a core change) and the clock that meets it.
  real time, capture/compare mode, 7000 rpm, at the minimum clock .... every event and the worst case
                                                                     <= 0.5 deg (slow-marked)
  real time, capture/compare mode, 2000 rpm tip-in, 50 MHz ........... the whole spec (slow-marked)
Tests that need riscv64-unknown-elf-gcc / verilator / the Vortex-V checkout skip without them.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from hil import rv64
from redline.spec import ROOT

RV64_DIR = ROOT / "ecu" / "targets" / "rv64"
ADR = ROOT / "docs" / "adr" / "0007-capture-compare-hal.md"
CORE_PATHS = ("ecu/src", "ecu/include", "ecu/host")
SPEC_DEG = float(rv64.cosim_spec()["realtime"]["angle_error_max_deg"])
BOARD_HZ = float(rv64.cosim_spec()["core_clock_hz"])
CC_SLOTS = 8                     # vortex_capcmp N_CMP (Vortex-V ADR 0017)
LOAD_MAX = 0.60
GAP_DEG_MAX = 1.0                # upper bound of the documented 50 MHz worst case
MIN_CLOCK_HZ = 80.0e6            # lowest whole-MHz clock meeting the spec at 7000 rpm (ADR 0007 sweep)
START_INJ_DIFF_MAX = 4           # one start-up injection per cylinder (as in M9.1)
TOOTH_TO_IRQ_MAX_CYCLES = 64     # pin edge -> interrupt taken: synchroniser + pipeline drain

needs_rv = pytest.mark.skipif(bool(rv64.toolchain_missing()) or not rv64.VORTEX_V_ROOT.exists(),
                              reason="needs riscv64-unknown-elf-gcc, verilator and $VORTEX_V_ROOT")


@pytest.fixture(scope="module")
def bench():
    from hil import realtime
    return realtime.Bench()


@pytest.fixture(scope="module")
def redline_50(bench):
    return bench.run("redline_rt", "capcmp", BOARD_HZ)


# --------------------------------------------------------------------------------------------
# anchors
# --------------------------------------------------------------------------------------------
def test_adr_0007_present_and_accepted():
    text = ADR.read_text()
    assert text.startswith("# ADR 0007")
    assert "Status: Accepted" in text
    for word in ("CAPCMP", "external interrupt", "compare slot", "knock", "bit-exact", "arm-ahead", "Consequences"):
        assert word in text, word


def test_portable_core_unchanged():
    """ADR 0003 / 0007: the capture/compare work is HAL, co-sim and build only. Checks the working
    tree and, when `main` is known, every commit of this branch."""
    if not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    st = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", *CORE_PATHS],
                        capture_output=True, text=True).stdout
    assert st.strip() == "", f"portable core / host HAL modified:\n{st}"
    base = subprocess.run(["git", "-C", str(ROOT), "merge-base", "HEAD", "main"], capture_output=True, text=True)
    if base.returncode == 0:
        diff = subprocess.run(["git", "-C", str(ROOT), "diff", "--stat", base.stdout.strip(), "HEAD", "--", *CORE_PATHS],
                              capture_output=True, text=True).stdout
        assert diff.strip() == "", f"portable core / host HAL changed on this branch:\n{diff}"


def test_model_built_with_the_capture_compare_unit():
    """One SoC model serves every mode, and it has the unit (CAPCMP=1, Vortex-V ADR 0017)."""
    soc = (RV64_DIR / "soc.mk").read_text()
    assert "CAPCMP=1" in soc
    vlt = (ROOT / "hil" / "cosim" / "public.vlt").read_text()
    for sig in ("wr_push", "pending", "cmp_drop", "cap_ovf", "fifo_ts"):
        assert f'-var "{sig}"' in vlt, sig


@needs_rv
def test_m9_and_m91_builds_still_build():
    """The comparison firmware stays buildable: the M9 code generation and the M9.1 mailbox mode
    (same image as M9.2, selected at run time by HIL_RT_F_CAPCMP)."""
    env = {**os.environ, "VORTEX_V_ROOT": str(rv64.VORTEX_V_ROOT)}
    subprocess.run(["make", "-s", "-C", str(RV64_DIR), "all"], check=True, env=env)
    subprocess.run(["make", "-s", "-C", str(RV64_DIR), "m9"], check=True, env=env)
    syms = rv64.elf_symbols()
    for name in ("irq_handler", "g_mbox", "ecu_crank_edge"):
        assert name in syms
    assert rv64.M9_ELF.exists()
    assert rv64.mbox_enums()["HIL_RT_F_CAPCMP"] not in (rv64.mbox_enums()["HIL_RT_F_IRQ"],
                                                        rv64.mbox_enums()["HIL_RT_F_PREDICT"])


# --------------------------------------------------------------------------------------------
# frozen time: bit-exactness is kept on the CAPCMP=1 SoC
# --------------------------------------------------------------------------------------------
@needs_rv
def test_frozen_time_still_bit_exact(bench):
    """The frozen-time lockstep mode on the M9.2 model and firmware: events and slow-loop outputs
    bit-exact against the host reference build."""
    from hil import cosim
    ecu = bench.ecu("fast")
    r = cosim.run_window("smoke", ecu, with_glibc=False)
    assert r["events"] > 0 and r["slow_steps"] > 0 and r["status"]["sync"] == 2
    assert r["diff"]["host_ref"] == {"events": 0, "slow": 0, "first": None}, r["diff"]["host_ref"]["first"]


# --------------------------------------------------------------------------------------------
# real time on the capture/compare unit
# --------------------------------------------------------------------------------------------
def _held(m):
    assert m["captures_dropped"] == 0 and m["cc_captures_lost"] == 0 and m["cap_ovf"] == 0
    assert m["cc_capture_tick_mismatches"] == 0                       # the RTL latched the host's tick
    assert m["captures_unserviced"] <= 1
    assert m["n_loss"] == 0 and m["n_dropout"] == 0 and m["n_noise"] == 0 and m["n_cam_err"] == 0
    assert m["sync"] == 2 and m["n_sync"] == m["host_n_sync"] == 1
    assert m["slow_overruns"] == 0 and abs(m["slow_steps"] - m["slow_steps_expected"]) <= 1
    assert m["spark_count_difference"] == 0 and m["sparks_compared"] > 0
    assert abs(m["events"] - m["host_events"]) <= START_INJ_DIFF_MAX
    # slots: never more than the unit has, nothing dropped by it, nothing held back by the HAL
    assert m["cc_pending_max"] <= CC_SLOTS and m["cc_cmp_drop"] == 0 and m["cc_deferred"] == 0
    assert m["slots_ok"] and m["cc_pushes"] == m["cc_pushes_fw"] == m["events"]
    # every compare accounted for on the pins
    assert m["cc_edge_missing"] == 0 and m["cc_pulse_mismatches"] == 0 and m["cc_pulses_checked"] > 0
    assert m["tooth_to_irq_max_cycles"] <= TOOTH_TO_IRQ_MAX_CYCLES


@needs_rv
def test_capcmp_redline_at_board_clock(redline_50):
    """7000 rpm, 50 MHz: every edge captured by the RTL with the host's tick, no loss of sync, the
    slow loop keeps its period, at most 8 compares pending, every compare's pin edge accounted for,
    every spark within 0.5 deg of the host's."""
    from hil import realtime
    r, m = redline_50
    print("\n" + realtime.line(m))
    print(f"max pending compares {m['cc_pending_max']} of {m['cc_slots']}; cycles with n pending: "
          f"{m['cc_pending_hist']}")
    _held(m)
    assert m["rpm"] == pytest.approx(7000.0)
    assert m["rt_ok"]
    assert m["cpu_load"] <= LOAD_MAX
    assert m["spark_angle_error_max_deg"] <= SPEC_DEG
    assert m["tooth_to_serviced_max_us"] < m["tooth_period_us"]


@needs_rv
def test_capcmp_on_time_compares_fire_on_their_tick(redline_50):
    """The unit removes output jitter: every event pushed before its tick changes its pin in the
    clock after TIME reaches that tick (the edge is found by exactly that tick), and every event
    pushed late is one the scheduler armed after its tick."""
    r, m = redline_50
    on_time = [e for e in r.events if not e["pushed_late"] and not e["no_edge"]]
    assert on_time and all(e["edge_tick"] == (e["tick"] + m["cc_t0"]) & 0xFFFFFFFF for e in on_time)
    late = [e for e in r.events if e["pushed_late"]]
    assert all(e["arm_cycle"] / r.clk_hz > (e["tick"] - 0.5) / 1.0e6 for e in late)
    assert m["events_late"] <= m["cc_pushed_late"]


@needs_rv
def test_capcmp_redline_known_gap_at_board_clock(redline_50):
    """KNOWN GAP (ADR 0007): at 7000 rpm and 50 MHz the hardware compare does not close the M9.1
    gap. Events right behind their tooth are still armed after their tick by the handler's
    latency; the worst one and the worst-case bound are above 0.5 deg and below 1 deg."""
    _, m = redline_50
    worst = max(m["worst_case_angle_error_deg"], m["event_angle_error_max_deg"])
    assert m["cc_pushed_late"] > 0
    assert SPEC_DEG < worst < GAP_DEG_MAX
    assert not m["meets_spec"] and m["rt_ok"]


@needs_rv
@pytest.mark.slow
def test_capcmp_redline_meets_spec_at_min_clock(bench):
    """At the minimum clock of the ADR 0007 sweep the whole spec holds at 7000 rpm: every event
    and the worst-case latency within 0.5 deg; one MHz below, it does not."""
    from hil import realtime
    _, m = bench.run("redline_rt", "capcmp", MIN_CLOCK_HZ)
    print("\n" + realtime.line(m))
    _held(m)
    assert m["event_angle_error_max_deg"] <= SPEC_DEG and m["worst_case_angle_error_deg"] <= SPEC_DEG
    assert m["meets_spec"]
    _, below = bench.run("redline_rt", "capcmp", MIN_CLOCK_HZ - 1.0e6)
    assert not below["meets_spec"]


@needs_rv
@pytest.mark.slow
def test_capcmp_tipin_at_board_clock(bench):
    """2000 rpm tip-in, 50 MHz: the whole spec holds on the capture/compare unit."""
    from hil import realtime
    _, m = bench.run("tipin", "capcmp", BOARD_HZ)
    print("\n" + realtime.line(m))
    _held(m)
    assert m["cpu_load"] <= LOAD_MAX
    assert m["event_angle_error_max_deg"] <= SPEC_DEG and m["worst_case_angle_error_deg"] <= SPEC_DEG
    assert m["meets_spec"]
