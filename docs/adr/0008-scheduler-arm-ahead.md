# ADR 0008: Scheduler arm-ahead: arm an event from the tooth before when it is close behind a tooth

Status: Accepted
Date: 2026-11-25
Supersedes: the event arming contract of ADR 0003 (the window an event is armed from). The rest of
ADR 0003 (HAL boundary, angle units, event struct, generated calibration) stands.

All figures below are MODELLED (MVEM plant, Verilator model of the Vortex-V SoC at 21fe1b1, no
physical engine). They are reproduced by `make -C ecu test`, `make -C hil cosim`, `make -C hil cc-sweep`
and checked in `tests/test_m8_ecu.py`, `tests/test_m9_1_rt.py` and `tests/test_m9_2_capcmp.py`.

## Context
ADR 0003 left the arming rule to `ecu_sched_tooth`: at each tooth, arm the events whose target lies
in `[this tooth, next real tooth)`. An event therefore had zero lead over the tooth that armed it
when it fell right at that tooth. On a target, the edge-to-armed latency of the handler (about 19 us,
0.80 deg at 7000 rpm on the 50 MHz Vortex-V) made such events late: 0.70 deg on the M9.1 mailbox
peripheral (ADR 0005) and 0.72 deg on the RTL capture/compare unit (ADR 0007), worst case 0.80 deg,
against the 0.5 deg spec. ADR 0007 showed the gap is arming latency, not output timing, and rejected
HAL look-ahead pre-arming (a predicted tooth tick moves an injection start by about 0.95 deg per
tick of misprediction, and the unit has no per-slot cancel). The fix belongs in the core, which must
stay deterministic and integer-only so the host and RV64 builds stay bit-exact.

## Decision
1. **Window.** At tooth N (cycle angle `A_N`, span `S_N` to the next real tooth: one pitch, three
   across the missing-tooth gap) the scheduler arms every event whose target lies in
   `[A_N, A_N + S_N + G)`. The part beyond the next tooth, `[A_{N+1}, A_{N+1} + G)`, is extrapolated
   with tooth N's measured period (the existing per-pitch period averaged over `ECU_PERIOD_HIST`
   intervals, Q8 sub-tick). The windows of consecutive teeth overlap by G; the existing repeat guard
   (360 deg) keeps each event to one arming per cycle, so an event lying less than G behind tooth N+1
   is armed once, from tooth N, at least G ahead of its target. The late catch-up (a target that slipped
   up to one pitch behind the current tooth fires at once) is unchanged. No state and no float are
   added; the guard costs one more integer division by the tooth period, the divisor the scheduler
   already uses for the dwell (on RV64 a cached reciprocal, ADR 0005), so the host and RV64 builds
   compute the same events and frozen-time co-sim stays bit-exact.
2. **G is a time, generated from the spec.** The handler latency it must cover is a time, and so is
   most of the target motion (point below), so the guard is `ECU_ARM_AHEAD_TICKS` converted to angle
   with the tooth period (as the dwell is), capped at `ECU_ARM_AHEAD_MAX_U` = one pitch - 1 unit, so a
   window never reaches two teeth beyond the arming one. Both constants are generated into
   `ecu_calib_gen.h` by `ecu/tools/gen_calib.py` from `spec/ecu.yaml` `decoder.arm_ahead_s`.
   **arm_ahead_s = 100 us** (100 ticks): 4.2 deg at 7000 rpm, 1.2 deg at 2000 rpm, 0.9 deg at
   1500 rpm. Sizing, from measurement:
   - edge-to-armed latency of the RV64 handler at 7000 rpm, 50 MHz: up to 20.5 us (M9.3 capture/compare
     run; 19.1 us in M9.2, the guard adds about 30 cycles). The unit stamps a capture 2 clocks after
     the pin edge and a compare output changes 1 clock after TIME reaches its tick (Vortex-V ADR 0017;
     the co-simulation models both): together 60 ns at 50 MHz, 0.003 deg at 7000 rpm, absorbed by G
     with nothing to change;
   - the tooth-to-tooth motion of injection-start targets. The target is `EOI - pulse`, the pulse
     converted to angle with the period that includes the new tooth: 1 us timer quantisation of a
     four-interval period moves a 12 ms pulse by up to ~1 deg per tooth at 7000 rpm (same pulse), and
     the 1 kHz slow loop changes the pulse between teeth. An injection whose target moves earlier by m
     between tooth N and N+1 is armed with lead G - m, so G must exceed latency + m. Host runs of the
     same core (redline WOT, the M8 tip-in, a free-rev into the rev limiter, 4500 rpm WOT; about 1300
     injections) with fixed-angle candidates: 1.5 deg left 24 of 229 redline injections with less than
     1.5 deg lead and 11 fired by catch-up; 3 deg left 2 below 0.8 deg; 4 deg none below 2.8 deg except
     the first injection after sync (its pulse appears when the target is already behind the tooth,
     which no window can arm earlier), i.e. up to ~2.3 deg (55 us) of motion at 7000 rpm. With the
     100 us guard the smallest injection lead in the same runs is 35 us (free-rev into the limiter;
     55-59 us elsewhere), sparks and dwell starts at least 97 us, against 20.5 us of latency.
   - **Why a time and not a fixed angle:** a fixed 4 deg was tried first and met the event spec, but it
     moved injection decisions up to 0.44 ms earlier at 1500 rpm, where a 4 deg target motion would take
     4x longer than at 7000 rpm, and the M9 supervisor test with the transient-seeking adversary went
     from lambda 1.1964 to 1.2122 on the 1500 rpm tip-in, outside its 0.8-1.2 safe band
     (`tests/test_m9_hil.py::test_supervisor_contains_adversarial_policies`). With the 100 us guard it
     is back to 1.1964, identical to M9.2.
   - G costs extrapolation (point 3), compare-slot occupancy (Consequences) and sync-loss exposure
     (point 4), all small; beyond one pitch it would not help (a target that moves more than a pitch
     between teeth is a pulse-width step, not jitter).
3. **Acceleration.** An event d deg past its arming tooth, timed with the mean speed over the last
   `hist` deg of teeth, is off by `alpha d (hist + d) / (2 omega^2)` under a constant acceleration alpha
   (the mean speed lags the instantaneous one by `alpha t_hist / 2`), plus timer quantisation (half a
   tick on the event and on the arming capture, one tick on the averaged span scaled by d / hist). The
   largest d is across the gap: 3 pitches + G from tooth 57, whose history is 4 normal pitches
   (24 deg). At the M8 acceleration case, 12000 rpm/s at 5800 rpm (G = 3.5 deg, d = 21.5 deg), the
   extrapolation term is 0.029 deg (0.022 deg without the window: arm-ahead adds 0.007 deg); at
   7000 rpm (G = 4.2 deg) 0.021 deg. With quantisation the bound is under 0.1 deg.
   `ecu/tests/test_ecu.c` (`test_arm_ahead_window`) sweeps spark targets through every 1/16 deg of a
   pitch, also right after the gap, at 7000 rpm steady, 7000 rpm with -12000 rpm/s and 5800 rpm with
   +12000 rpm/s: every spark is armed from a tooth at least the guard before it (within 0.05 deg: the
   guard is converted with the averaged speed) and at most span + guard before it, one spark per
   cylinder per cycle, every error within the analytic bound and within 0.5 deg (worst 0.078 deg, bound
   there 0.099 deg).
4. **Missing-tooth gap, loss of sync, rev limiter.**
   - Gap: the window of tooth 57 is `[A_57, A_0 + G)`, 18 deg + G; events inside the gap have always
     been armed from tooth 57, since there is no tooth to arm them from. Tooth 0's window starts at
     `A_0`; its first G are tooth 57's. The windows tile the cycle across the gap.
   - Loss of sync: nothing is armed from a tooth without full sync (unchanged), so nothing is armed
     from the tooth that detects a loss or until the cam re-establishes the cycle angle. Events armed
     before the loss was detected are not cancelled: the HAL has no cancel and the unit has no per-slot
     cancel (ADR 0007). They are bounded: at most `S_N + G` past the last good tooth, i.e. at most G
     beyond what ADR 0003 already armed (the gap's events were always armed 18 deg across the tooth at
     which a missing gap is detected). `test_arm_ahead_sync_loss` checks both (an accepted spurious
     edge, resync on the gap at the wrong count).
   - Rev limiter: the cut flag is evaluated at the arming tooth and decides each event once (a cut
     event is marked fired, so the next tooth does not arm it). An event armed ahead therefore follows
     the cut state of the tooth before its own: cut onset and release move by up to G of crank angle,
     against a cut decided on a full revolution's period with 150 rpm of hysteresis. `test_rev_limiter` still holds.

## Consequences
- **Frozen time:** events and slow-loop outputs bit-exact against the host reference build and the
  glibc host build on both frozen-time windows (`make -C hil cosim`).
- **Fixtures:** no committed host event fixture exists (the co-simulation compares against the host
  build of the same sources at run time). The one regenerated committed artefact is
  `ecu/include/ecu/ecu_calib_gen.h` (adds `ECU_ARM_AHEAD_TICKS` and `ECU_ARM_AHEAD_MAX_U`, from the new
  `spec/ecu.yaml` key `decoder.arm_ahead_s`).
- **Real time at 50 MHz, M9.2 -> M9.3 (capture/compare unit; the M9.1 mailbox mode, same core,
  in brackets):**

  | | 7000 rpm, M9.2 | 7000 rpm, M9.3 | 2000 rpm tip-in, M9.2 | 2000 rpm tip-in, M9.3 |
  |---|---|---|---|---|
  | any event angle error, max | 0.72 deg [0.70] | 0.27 deg [0.25] | 0.22 deg [0.21] | 0.005 deg [0.000] |
  | events late (compare pushed after its tick) | 9 of 74 [9] | 0 of 74 [0] | 4 of 60 [4] | 0 of 60 [0] |
  | worst tooth-to-event latency | 19.1 us | 20.5 us (0.86 deg) | 20.2 us | 21.1 us |
  | guard G at that speed | - | 4.2 deg | - | 1.2 deg |
  | smallest lead of an event over its arming tooth | 0 | 3.78 deg | 0 | 1.42 deg |
  | worst case (latency - smallest lead, >= 0) | 0.80 deg [0.75] | 0 (margin 2.9 deg) | 0.24 deg [0.23] | 0 |
  | worst event vs host: dwell / spark | 0.024 / 0.024 deg | 0.270 / 0.060 deg | 0.005 / 0.005 deg | 0.005 / 0.005 deg |
  | max pending compares (of 8) | 2 | 2 | 2 | 2 |
  | cycles with 2 compares pending | 3214 | 17759 | 20842 | 61158 |
  | CPU load | 48.0 % [47.9] | 48.5 % [48.4] | 38.4 % [37.7] | 38.5 % [37.8] |
  | captures lost / CAP_OVF / loss of sync / slow overruns | 0/0/0/0 | 0/0/0/0 | 0/0/0/0 | 0/0/0/0 |

  "Worst case" changes meaning with this ADR: it was the worst latency as an angle (an event could sit
  right at its arming tooth); it is now the worst latency minus the smallest lead an event can get
  (G, or the lead an event actually got when that was less: catch-up, first injection after sync).
  The 0.27 deg figure at 7000 rpm is not lateness: it is one dwell start shortly after sync that the
  target armed 6 us (6 ticks) earlier than the host. The slow loop's outputs
  reach the scheduler at a different tooth on target than on the host (ADR 0005's snapshot commit),
  and with arm-ahead that dwell is decided one tooth earlier, from the older advance. The 0.060 deg
  spark is the same effect (1 tick). Both are inside the spec; neither is a late compare.
- **Every event within 0.5 deg at 7000 rpm and 50 MHz: met** in both the capture/compare and the
  mailbox modes. The minimum clock is no longer set by the angle spec but by the CPU load (real time
  is lost when slow-loop steps overrun). **Clock sweep at 7000 rpm (M9.3, capture/compare):** 10, 20,
  21, 24 and 25 MHz lose real time (100 / 53 / 52 / 48 / 42 slow steps skipped); 26 MHz holds it and
  meets the spec (every event 0.27 deg, worst case 0, CPU load 92.6 %, slow-loop step 995 us of its
  1000 us period: almost no margin); every whole MHz from 27 to 32, and 50, 55 and 100 MHz meet it too.
  **Minimum whole-MHz clock: 26 MHz** (M9.2: 80 MHz). With the M9.1 load bound of 60 % as margin the
  practical floor is about 40 MHz. The 0.27 deg event error at 26-55 MHz is the slow-loop timing effect
  above, not arming latency. **2000 rpm tip-in:** meets the spec at 26 MHz (worst case 0.49 deg) and
  50 MHz; at 25 MHz it holds real time but misses the spec by 0.007 deg (worst case 0.507 deg: one
  injection start fired by late catch-up during the pedal step, observed 0.28 deg late, bounded by the
  42 us latency; M9.2: 0.49 deg there, the guard's division made the handler about 3 % longer); 20 MHz
  loses real time. A pulse-width step larger than the guard is not something any arm-ahead window can
  pre-arm; the catch-up remains for it.
- **Slots:** arm-ahead keeps an event pending up to G longer, so two compares are pending more often
  (5.5x the cycles at 7000 rpm, 2.9x in the tip-in), but never more than 2 of 8; nothing dropped or
  held back.
- The worst-case metric and the arming-tooth attribution in `hil/realtime.py` changed with this ADR
  (an event's arming tooth is the crank capture in service when it was pushed, no longer the last
  tooth before its tick).
- Late catch-up stays for targets that jump behind the tooth by more than the guard (the first
  injection after sync, a pulse-width step in a fast transient), late by the handler latency. At 50 MHz
  it did not occur in the measured windows; at 25-26 MHz it occurred once in the tip-in. Removing it
  would change start-up and transient fuelling (M8), so it is left as is.
- **Host-side regression check:** with the time guard every M8 / M9 / M10 test passes unchanged,
  including the M9 supervisor's adversarial policies (lambda 1.1964 on the 1500 rpm tip-in, as in
  M9.2). The M8 scheduler figures in `make -C ecu test` move from 0.186 / 0.280 deg (steady 7000 rpm /
  5800 rpm + 12000 rpm/s) to 0.192 / 0.288 deg, against 0.5 deg.
- Vortex-V has not verified FPGA timing at any of these clocks; the capture/compare pins are not on
  Arty pins yet (ADR 0007).
