# ADR 0003: ECU HAL boundary, fixed-point angle units and event contract

Status: Accepted
Date: 2026-11-21

## Context
The same ECU core must run on the host (ctypes HIL, C unit tests) and, in M9, as RV64 code on
Vortex-V in Verilator co-sim. Crank-synchronous code must be deterministic and float-free so
it behaves identically on a core without an FPU path in the HIL budget.

## Decision
- The portable core (`ecu/src`) calls only `ecu/include/ecu/ecu_hal.h` (`hal_emit_event`). Targets
  provide the HAL (`ecu/host/hal_host.c` today; an RV64 HAL in M9) and call the core entry points in
  `ecu.h` from crank/cam capture, knock-window and slow-loop contexts.
- Angles in the crank-synchronous path are integers: 64 units per crank degree
  (`ECU_ANGLE_U_PER_DEG`), one cycle = 46080 units (`ECU_CYCLE_U`), 0 = firing TDC of cylinder 1
  (ADR 0001). Decoder, scheduler and knock units are built with the float-forbidding flag set by
  the Makefile `NOFLOAT` variable; the slow loop (fuel, lambda, throttle, idle) may use float.
- Events are `ecu_event_t {tick, arg, type, cyl}` with absolute wrap-safe 32-bit timer ticks:
  INJ opens at `tick` for `arg` ticks; DWELL starts charging; SPARK cuts coil current.
- Calibration tables are generated from spec/thermo outputs into `ecu_calib_gen.h` by
  `ecu/tools/gen_calib.py`; the header is committed and checked for reproducibility.

## Consequences
M9 only adds a HAL and overrides `NOFLOAT` for RV64; core sources stay unchanged.
