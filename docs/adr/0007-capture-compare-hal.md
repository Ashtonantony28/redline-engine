# ADR 0007: RV64 HAL on the Vortex-V capture/compare timer

Status: Accepted
Date: 2026-11-25

All figures below are MODELLED (MVEM plant, Verilator model of the Vortex-V SoC, no physical engine).
They are reproduced by `python -m hil.realtime` (`make -C hil cosim`) and
`python -m hil.realtime --cc-sweep` (`make -C hil cc-sweep`), and checked in `tests/test_m9_2_capcmp.py`.

## Context
M9.1 (ADR 0005) ran the unchanged ECU core on Vortex-V in real time. The SoC had no capture/compare
hardware, so the Verilator wrapper played one through a RAM mailbox, and tooth handling was woken by
`MTIMECMP` for the predicted next edge, with a bounded spin. At 7000 rpm and 50 MHz the worst event
was 0.76 deg late against the 0.5 deg spec: the portable scheduler arms an event from the tooth just
before it, and the handler needs up to about 18 us from the edge to the last event store.

Vortex-V now has a capture/compare timer in RTL (its ADR 0017, merged in PR #17; built with
`vortex_soc #(.CAPCMP(1))`): a 32-bit timer at PRESC+1 clocks per tick, a 32-entry capture FIFO on
`cap_in[1:0]` (crank, cam) with a machine external interrupt while it is not empty, and 8 one-shot
compare slots on `cmp_out[7:0]` with pulse / set / clear actions. It has no knock input. ADR 0017 also
changed the software contract: `crt0` now calls `irq_handler(uint64_t mcause)`, which the M9.1 HAL
did not compile against.

## Decision
1. **One firmware image, a third run-time mode.** `HIL_RT_F_CAPCMP` (with `HIL_RT_F_IRQ`) selects the
   capture/compare mode inside the real-time command (`hil_mbox.h`). The M9 frozen-time command mode
   and the M9.1 mailbox mode stay in the same image and are selected the same way; `make -C
   ecu/targets/rv64 m9` still builds the M9 code generation. The portable core (`ecu/src`,
   `ecu/include`, `ecu/host`) is unchanged (ADR 0003). One SoC model, built with `CAPCMP=1`, serves
   every mode: the M9 / M9.1 modes leave the unit idle.
2. **Captures from the RTL FIFO under the external interrupt.** `irq_handler` (`mcause` = 2^63 + 11)
   drains the FIFO in order, calling `ecu_crank_edge` / `ecu_cam_edge` with the latched tick, and pops
   each entry after servicing it (a pop is a `CAP_RD` write; reads have no side effects). `MTIMECMP`,
   the next-edge prediction and the spin are gone from tooth handling. The timer runs at the ECU's 1 MHz:
   the firmware writes `PRESC = f_clk / ECU_TIMER_HZ - 1`, so clocks are whole MHz.
3. **One compare slot per event.** `hal_emit_event` pushes: INJ = pulse of `arg` ticks on output `cyl`,
   DWELL = set and SPARK = clear on output `ECU_NCYL + cyl` (the mapping of the Vortex-V register map).
   **The 8-slot limit is enforced in software, not by the unit's drop counter:** the HAL keeps a count
   of slots known to be free (decremented per push; refreshed from `CMP_STATUS` only when it reaches 0;
   slots only free on their own, so it is a lower bound). An event that finds no free slot goes to a
   backlog, is pushed on the next capture, and is counted (`HIL_RT_S_CC_DEFERRED`); the unit's drop
   count is read back at the end of every run. Both must be 0.
4. **Knock stays on the mailbox.** The unit has no knock input; knock windows stay in the mailbox
   FIFO with the timer tick of the window's end. The handler merges them with the hardware captures in
   tick order (a knock window strictly older than the oldest hardware capture goes first), so the core
   sees the order the host build sees. They are serviced in the next capture interrupt.
5. **Slow-loop commit is atomic against captures.** The slow loop still runs on a snapshot (ADR 0005),
   but its commit of advance, pulse width and limp flag runs with `mstatus.MIE` clear, so a tooth sees
   either all old or all new outputs, never one cylinder's new advance with its old pulse width.
6. **Co-simulation on the pins.** The plant drives `cap_in` (an 8-cycle pulse) at the cycle a crank or
   cam edge occurs and acts on `cmp_out` edges (injector open while high, fuel credited on closing and
   at the cylinder's IVC if still open; spark on the coil's falling edge). The model observes the RTL
   only: every `cmp_out` edge (cycle, TIME), every `CMP_PUSH` as the unit takes it, every capture and
   pop, slot occupancy per cycle, `CAP_OVF` and the push-drop count. Plant time 0 is placed on the
   timer's phase so that a captured tick equals `round(t * 1 MHz) + T0`, the tick the host build is
   given for the same edge (checked for every edge). Each push is matched to the pin edge it caused;
   one that changes nothing (a clear on a coil that never charged) is kept with the unit's fire time
   and checked against the pin level.
7. **Rejected: HAL look-ahead pre-arming.** Prototyped and measured, not kept. After each tooth the
   handler ran the core once more on a copy of its state with the predicted next tooth tick and
   pre-armed the events that fall within 1500 cycles behind that tooth (the compare slots make this
   possible). At 7000 rpm, 50 MHz no event was pushed late and the largest deviation from the real
   pass was 10 ticks (0.42 deg), but: the scheduler converts the injector pulse to angle with the
   period that includes the new tooth, so a 1-tick misprediction (the timer's quantisation, not
   avoidable with 1 us ticks) moves an injection start by about pw / (4 pitches) ~ 23 ticks, 0.95 deg
   at 7000 rpm full load: worse than the latency it removes. A pre-armed compare that turns out wrong
   usually fires before the real pass can correct it, and the unit has no per-slot cancel. It also
   cost about 1010 cycles per tooth (CPU load 48 % -> 61 %) and changed when slow-loop outputs reach the
   scheduler. Prediction-based arming is not exact; the fix belongs in the core (below).

## Consequences
- **Frozen time:** events and slow-loop outputs remain bit-exact against the host reference build on
  the `CAPCMP=1` SoC.
- **Real time at 50 MHz, M9.1 (mailbox peripheral, re-measured on Vortex-V main) -> M9.2 (RTL unit):**

  | | 7000 rpm, M9.1 | 7000 rpm, M9.2 | 2000 rpm tip-in, M9.1 | 2000 rpm tip-in, M9.2 |
  |---|---|---|---|---|
  | worst tooth-to-serviced latency | 46.8 us (start-up poll) | 20.2 us | 55.5 us (start-up poll) | 20.8 us |
  | worst tooth-to-event latency | 18.0 us | 19.1 us | 19.1 us | 20.2 us |
  | pin edge to interrupt taken, max | (woken ahead, spins) | 32 cycles | (woken ahead, spins) | 21 cycles |
  | spark angle error vs host, max | 0.042 deg | 0.024 deg | 0.000 deg | 0.005 deg |
  | any event angle error, max | 0.70 deg | 0.72 deg | 0.21 deg | 0.22 deg |
  | worst case (latency as crank angle) | 0.75 deg | 0.80 deg | 0.23 deg | 0.24 deg |
  | CPU load | 47.9 % | 48.0 % | 37.7 % | 38.4 % |
  | captures lost / CAP_OVF / loss of sync | 0 / - / 0 | 0 / 0 / 0 | 0 / - / 0 | 0 / 0 / 0 |
  | max pending compares (of 8) | - | 2 | - | 2 |
  | compare pushes dropped / held back | - | 0 / 0 | - | 0 / 0 |

- **Not met at 50 MHz, 7000 rpm: 0.5 deg for every event.** The unit fires every compare that is armed
  before its tick exactly on that tick, and removes the prediction, the spin and the start-up poll. It
  does not remove the M9.1 gap, because the gap was never output jitter (M9.1's played compare already
  fired at `max(tick, arm time)`, which is what the RTL does): it is the scheduler arming an event
  from the tooth immediately before it. The late events are injection starts the scheduler fires
  "now" at the tooth (target slipped behind the tooth, ECU late catch-up) or within the handler's
  latency after it. The handler is the same core code, plus the interrupt entry the predicted wake
  used to hide (up to 32 cycles to the trap, then the register save), so the worst case is slightly
  worse than M9.1 (0.80 vs 0.75 deg).
- **Clock sweep at 7000 rpm (M9.2):** 20 MHz loses real time (52 slow steps skipped, worst case 1.95 deg);
  25 MHz loses real time (35 skipped, 1.61 deg); 50 MHz holds real time, 0.80 deg; 55 MHz 0.72 deg;
  78 MHz 0.51 deg; 79 MHz 0.501 deg; **80 MHz is the lowest whole-MHz clock that meets the spec**
  (0.495 deg; M9.1: 77.5 MHz on a 2.5 MHz grid); 100 MHz 0.40 deg. The CPU load at 25 MHz is lower than
  M9.1's (35 instead of 46 slow steps skipped: no spin), but not low enough to hold real time. The 2000 rpm tip-in meets the spec at 50 MHz
  (0.24 deg) and at 25 MHz (0.49 deg), and loses real time at 20 MHz.
- **Slots:** at most 2 of 8 are ever pending in either window (an event is armed at most one tooth
  window, 3 pitches at the gap, ahead of its tick, and an injector pulse frees its slot when it starts).
  8 is enough at 7000 rpm with a factor of 4 to spare; the HAL backlog makes an overflow visible
  instead of silent.
- **Proposed fix for 50 MHz (not in M9.2's scope): scheduler arm-ahead.** Let `ecu_sched_tooth` arm,
  from tooth N, the events whose target lies within a guard angle `G` behind tooth N+1 (window
  `[tooth N, tooth N+1 + G)`), extrapolated with tooth N's period; the existing repeat guard stops tooth
  N+1 arming them again. It is deterministic, so the host build computes the same events and bit-exact
  co-simulation holds after a reviewed host-fixture regeneration. With `G` = 1 deg (the worst handler
  latency, 19.1 us, is 0.80 deg at 7000 rpm) no event is armed after its tick at 50 MHz. It needs a new
  ADR (it changes the core's arming contract of ADR 0003), a new host fixture, and about one more
  pending slot. Alternatively run the core at or above the minimum clock above.
- `cap_in` / `cmp_out` are not on Arty pins yet (Vortex-V ADR 0017 follow-on), and Vortex-V has not
  verified FPGA timing at any of these clocks.
- The co-simulation model reads the unit's internal signals (`hil/cosim/public.vlt`): it is tied to
  `vortex_capcmp`'s signal names as well as to `CORE_OOO=1`.
