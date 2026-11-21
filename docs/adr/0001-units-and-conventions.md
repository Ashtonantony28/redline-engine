# ADR 0001: Units, crank-angle and sign conventions

Status: Accepted
Date: 2026-11-21

## Context
Thermodynamics, gas dynamics, dynamics, FEA and the ECU all exchange crank-angle-indexed
data. A mismatch of origin or sign silently corrupts every downstream result.

## Decision
- SI units everywhere in code and data (m, kg, s, Pa, K, J, rad). Degrees only at I/O
  boundaries and in `spec/engine.yaml` keys suffixed `_deg`.
- Crank angle 0 = firing TDC of cylinder 1. One cycle spans 0..720 deg; 360 = gas-exchange TDC.
  Cylinder k (firing order 1-3-4-2) fires at phase offset 180 deg * (position in firing order).
- Piston displacement `s` is measured downward from TDC; positive crank torque drives the
  crank forward (engine rotation direction).
- Body axes for shaking forces: x along the crankshaft (front to rear), z along the bore
  axis (up), y completing a right-handed set.
- `spec/engine.yaml` is the single source of truth; no magic numbers in model code.

## Consequences
Every module converts degrees at its edges. Tests pin the origin (volume minimum at 0 deg).
