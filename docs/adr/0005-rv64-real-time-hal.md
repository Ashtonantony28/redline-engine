# ADR 0005: Real-time RV64 HAL: timer-interrupt sync path, snapshot slow loop, exact fast arithmetic

Status: Accepted
Date: 2026-11-21

All figures below are MODELLED (MVEM plant, Verilator model of the Vortex-V SoC, no physical engine).
They are reproduced by `python -m hil.realtime --sweep` and checked in `tests/test_m9_1_rt.py`.

## Context
M9 ran the ECU core (ADR 0003) on the Vortex-V RTL in lockstep: every entry point was a polled
mailbox command and plant time stood still while the core computed. That proves bit-exactness and
hides time. The M9 report listed two gaps at 7000 rpm on the 50 MHz Arty clock (tooth period
142.9 us = 7143 cycles): a tooth arriving during a slow-loop step waits for it (up to 973 us), and
the slow loop on libgcc soft-float (46.7k cycles per 1 ms) plus the teeth needs 107 % of the core.

Vortex-V now has a machine timer interrupt (its ADR 0015: `MTIMECMP` in the SoC, `irq_handler()`
from `crt0`, taken after the out-of-order pipeline drains). It is the only interrupt: there is no
external interrupt, no capture/compare timer and no FPU. Vortex-V is consumed read-only.

A retire-PC profile of the M9 firmware (`python -m hil.profile`) showed where the cycles went:
- crank-tooth handler, 1430 cycles: about 1300 in 20 integer divisions. Vortex-V divides one bit per
  clock (64 cycles) and GCC emitted `remw` for `% ECU_CYCLE_U` in `ecu_angle_wrap`, 12 times a tooth.
- slow loop, 44k cycles: 15.7k in `__divdf3` (47 binary64 divisions per step in the expf / powf
  series of `ecu_libm.c`), most of the rest in libgcc's generic soft-fp and its `__clzdi2`.
- every function return is a JALR that stops dispatch until it resolves, and the core speculates
  past one branch only, so small helper calls cost more than their work.

## Decision
1. **Two time models, one firmware.** The frozen-time command mode of M9 stays, unchanged in
   behaviour, for the bit-exactness tests. A real-time mode is added (`HIL_CMD_RT`,
   `hil/realtime.py`): plant time is core cycles. The Verilator model plays the capture/compare timer
   the SoC lacks, cycle by cycle, through the mailbox (`hil_mbox.h`): a free-running timer register,
   an input-capture FIFO that the plant fills at the cycle an edge occurs (a full FIFO drops the
   edge), and output compares that take effect at `max(event tick, cycle the event reached RAM)`.
   Sensor samples and slow-loop outputs are double-buffered. Lateness and overload therefore appear
   as timing errors in the events the plant acts on.

2. **Sync path in the machine timer interrupt, armed for the next predicted tooth.** `irq_handler()`
   drains the capture FIFO in order and calls `ecu_crank_edge` / `ecu_cam_edge` /
   `ecu_knock_window`. It then writes `MTIMECMP` for the predicted time of the next crank edge
   (last accepted tick + averaged period, three pitches at the gap) minus a lead of 150 cycles that
   covers the pipeline drain and the register save. If it wakes before the edge it waits for it
   inside the handler, at most 400 cycles; an edge that does not come (search, stall, a wrong
   prediction) falls back to a fixed 2500-cycle poll. The values are in `hil/cosim.yaml`.
   - Rejected: a fixed periodic tick (`--modes irq_periodic`). Measured at 7000 rpm, 50 MHz with a
     2500-cycle period: tooth serviced up to 68.9 us late, an event up to 2.7 deg late. A period short enough for 0.5 deg
     (under 600 cycles including the handler) would spend most of the core on interrupt entries.
   - Prediction error only moves latency, never correctness: the handler acts on captured ticks.
     Its size is the speed change over one tooth, far below one timer tick at any engine speed
     where latency matters.

3. **Slow loop in the background, on a snapshot.** `ecu_slow_step` cannot be made atomic by masking
   interrupts (it is longer than a tooth period) and must not race the handler: a tooth newer than
   `now_tick` reads as a stall and drops sync. The HAL keeps a second `ecu_t`. Each step it copies
   the handler-owned state (decoder, scheduler, knock retard) with a retry-if-preempted read, runs
   the unchanged `ecu_slow_step` on the copy with interrupts enabled, then commits advance, pulse
   widths and the limp flag with single-word stores and a stall verdict under a short mask, only if
   no edge arrived since the snapshot. `now_tick` is the firmware's own timer read after the snapshot.
   The step stays at `ECU_DT_SLOW` (1 kHz) for all four loops (fuel, lambda, throttle, idle): the
   period is compiled into the unchanged core, so a lower rate would change every integrator, and
   after item 4 it is not needed. A step that overruns its period skips to the next grid point and
   is counted.

4. **Cut the load by exact arithmetic, not by approximation.** RV64 target only, `FAST=1` (default):
   - `rv_softfp.c` replaces libgcc's soft-float: IEEE-754 binary32 / binary64, round to nearest
     even, written for a 64-bit core with a one-cycle multiplier, helpers inlined, divisions by a
     Newton reciprocal with an exact remainder correction (`rv_recip.h`).
   - `rv_div.c` + `-mno-div` on the sync path: constant divisions become multiplies, variable ones a
     reciprocal multiply with exact correction; the reciprocal of the last divisor is cached (the
     scheduler divides by the same tooth period five times a tooth).
   - `sync_unit.c` compiles decoder, scheduler and knock control (unchanged) as one translation
     unit at `-O3 -mbranch-cost=6`, so `ecu_angle_wrap` and friends inline.
   - `ecu_libm.c` evaluates its series with constant coefficients instead of dividing by the term
     index at run time. This file is shared by the target and the host reference build.
   Fixed point was not used: every result above is the exact one, so **host/target bit-exactness
   holds for every output** and no tolerance had to be introduced. `tests/softfp_check.c` compares
   the target arithmetic with the host FPU and integer divide bit for bit.
   `make -C ecu/targets/rv64 m9` keeps the M9 code generation as the "before" build.

## Consequences
- Frozen time: events and slow-loop outputs remain bit-exact against the host reference build (and,
  on the committed windows, against the glibc host build). Crank handler 1430 -> 705 cycles mean,
  slow loop 44.2k -> 16.4k cycles mean (363 us of its 1000 us period at 50 MHz).
- Real time at 50 MHz, M9 firmware -> M9.1 firmware:

  | | 7000 rpm, M9 | 7000 rpm, M9.1 | 2000 rpm tip-in, M9 | 2000 rpm tip-in, M9.1 |
  |---|---|---|---|---|
  | worst tooth-to-serviced latency | 965 us | 46.3 us (start-up poll; 19 us once predicting) | 956 us | 54.6 us |
  | worst tooth-to-event latency | 946 us | 18.1 us | 946 us | 18.9 us |
  | spark angle error vs host, max | 34.9 deg | 0.042 deg | 10.6 deg | 0.000 deg |
  | any event late by, max | 38.0 deg | 0.70 deg | 11.2 deg | 0.20 deg |
  | CPU load | 85.7 % with 35 of 120 slow steps skipped | 48.2 % | 95.1 % | 37.6 % |
  | slow-loop step, wall time max | 933 us | 430 us | 920 us | 360 us |
  | captures dropped / loss of sync | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

  (The M9 firmware loses no tooth only because the modelled capture FIFO holds 32 edges; a
  single capture register would have lost several per slow-loop step.)
- **Not met at 50 MHz, 7000 rpm: the 0.5 deg bound for every event.** The portable scheduler arms
  an event from the tooth immediately before it. The handler needs up to 907 cycles from the edge to
  the last event's store (18.1 us = 0.76 deg at 7000 rpm), so an event that falls within that
  distance behind a tooth is armed after its tick and fires late: observed 0.70 deg on injection and
  dwell starts in the measured window; the sparks there landed mid-window and were within 0.05 deg.
  Sparks have no protection from this: at another advance they would be late by the same amount.
  The latency is a cycle count, so it scales with the clock. Sweep at 7000 rpm: 25 MHz loses real
  time (46 of 120 slow steps skipped, worst case 1.50 deg); 50 MHz holds real time, worst case
  0.76 deg; 75 MHz 0.51 deg; **77.5 MHz is the lowest clock that meets the bound** (0.49 deg, 2.5 MHz
  search step); at the board oscillator's 100 MHz the worst case is 0.38 deg. Vortex-V states
  50 MHz as its Arty target and has not yet verified timing at any frequency. The 2000 rpm tip-in
  meets the bound at 50 MHz (0.23 deg) and at 25 MHz (0.46 deg).
- Closing the gap at 50 MHz needs one of: the scheduler arming an event from the previous tooth
  when it is closer than the handler latency (a change to the portable core, so a new ADR and a new
  host fixture); a hardware capture/compare unit with its own interrupt in the SoC (a Vortex-V
  change); or wider branch speculation in the core. None is in M9.1's scope.
- The HAL now reads decoder fields (`last_tick`, `period_q8`, `tooth`, `sync`, `n_edges`,
  `have_last`) and relies on `stft_i` being the first slow-loop field of `ecu_t`. A layout change
  in `ecu.h` must revisit `hal_rv64.c`.
- `rv_div.c` keeps a one-entry cache and is not reentrant: it is linked for the capture context only.
- NaN payloads and out-of-range float-to-integer conversions follow RISC-V / libgcc conventions
  and are the only results not compared with the host; the ECU produces neither.
- The co-simulation model reads the out-of-order core's retire and interrupt debug signals
  (`hil/cosim/public.vlt`), so the real-time model is tied to `CORE_OOO=1`.
