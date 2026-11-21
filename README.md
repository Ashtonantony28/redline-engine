# Redline

A code-first digital twin of a 2.0 L naturally aspirated inline-4: parametric CAD (build123d),
crank-angle thermodynamics and combustion (Cantera), 1D gas dynamics, kinematics, dynamics and NVH,
structural FEA, and an engine control unit that runs hardware-in-the-loop on the Vortex-V RISC-V core.

All numbers are **modelled**. Plan: see `docs/PLAN.md`.

## Quick start
```
pip install -r requirements.txt
pytest -q
```
