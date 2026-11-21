# ADR 0004: Turbocharger coupling and synthetic turbomachinery maps

Status: Accepted
Date: 2026-11-21

## Context
The M10 stretch (PLAN option B) adds a turbocharged variant, Redline I4-2.0T. It needs compressor and
turbine maps, which are hard to source openly for a machine of the right size, and it has to sit on the M4
1D gas-dynamics kernel (gasdyn.fv, numba) without changing the NA results that M5, M7, M10 and the report
are built on.

## Decision
- **Coupling outside the kernel.** Compressor, turbine, intercooler and wastegate are 0D elements applied
  between kernel calls (operator split every `steady.substep_deg`, gasdyn.elements), through optional hooks
  on gasdyn.engine.EngineModel (`boost` layout, `substep_fn`, `cycle_fn`). With the hooks unset the M4
  engine is built and marched exactly as before (tests/test_turbo.py checks the M4 fixture point).
- **Turbine as the kernel orifice.** The turbine-housing pipe ends in an orifice of the turbine's
  nozzle-equivalent area, so the kernel's isentropic orifice boundary is the turbine flow map; the shaft
  work (isentropic work x blade-speed-ratio efficiency) is withdrawn from the post-turbine volume.
- **Compressor as a quasi-steady source** into a 0D charge-air plenum, from a map inverted against the
  plenum pressure plus intercooler and air-filter losses.
- **Synthetic maps** from documented non-dimensional models (turbo/maps.py, constants in spec/turbo.yaml),
  labelled SYNTHETIC wherever results are shown.
- **Boost control in Python** (turbo/control.py): the C11 ECU has no wastegate output, its MAP range ends at
  1.15 bar, and its output struct is shared bit-exactly with the RV64 co-simulation.
- **Variant as an override layer** (spec/turbo.yaml over engine.yaml / gasdyn.yaml); engine.yaml untouched.

## Consequences
- Turbine and compressor dynamics are resolved at the sub-step (5 deg) rather than the kernel time step;
  the turbine work uses the end-cell state of the housing pipe, so the split has a small time-integration
  error but conserves mass exactly.
- Results depend on synthetic maps: absolute boost, lag and surge-margin numbers are indicative only.
- Moving the elements into the numba kernel later would be a new ADR superseding this one.
