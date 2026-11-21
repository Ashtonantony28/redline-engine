# ADR 0002: Crank rotation direction in body axes

Status: Accepted
Date: 2026-11-21

## Context
ADR 0001 fixes the crank-angle origin and the body axes but not which way the shaft turns
about +x. CAD placement, shaking moments and bearing-load polar plots all depend on it.

## Decision
Increasing crank angle theta turns the crankshaft by -theta about +x (clockwise seen from the
front, looking rearward along +x). The crankpin of cylinder k sits at
(y, z) = r * (sin(theta - phase_k), cos(theta - phase_k)), with phase_k from the firing order.

## Consequences
cad/assembly.py and cad/crankshaft.py implement this; dyn/ and later FEA load cases must use
the same sign when mapping crankpin forces into body axes.
