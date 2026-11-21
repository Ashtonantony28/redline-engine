# Redline: inline-4 engine digital twin (CAD, thermodynamics, dynamics, ECU)

Status: DRAFT plan, 2026-11-19. Working name "Redline", suggested repo `Ashtonantony28/redline-engine`.
Everything runs in the cloud container. No hardware, no dyno, nothing on the user's computer.

## 1. Options considered

| | Option | What it shows | Fit with existing work | Risk |
|---|---|---|---|---|
| **A (recommended)** | **Full digital twin of a 2.0 L NA inline-4: parametric CAD, crank-angle thermodynamics, 1D gas dynamics, kinematics and dynamics, structural FEA, and an ECU running in hardware-in-the-loop** | Mechanical design, thermo/combustion, vibration, FEA, controls, embedded firmware | Reuses drone-design's build123d + FEA flow and Zener E-CAD, NeuroFlight's PPO scaffold, and closes the loop on **Vortex-V** (ECU firmware and an int8 virtual sensor run as RV64 code in co-sim, same pattern as Talon M6) | Medium: physics is textbook (Heywood), tools are pip/apt installable |
| B | Turbocharged 1.5 L I4 (A plus compressor/turbine maps, wastegate, boost control) | Everything in A plus turbomachinery and boost control | Same as A | Medium-high: needs compressor and turbine maps that are hard to source openly. Kept as stretch milestone M10 |
| C | Wankel rotary or opposed-piston two-stroke | Unusual kinematics, strong "wow" factor | Weaker: less published validation data | High: few reference cases to check against |
| D | Hybrid powertrain (small ICE + e-motor + battery, energy management) | Systems/controls breadth | Overlaps Sentinel/NeuroFlight battery modelling | Medium, but the engine itself becomes a lookup table |

**Recommendation: A**, with B's turbo as a later stretch. A naturally aspirated engine has the most published reference data to validate against, every subsystem has an analytic check, and it gives a complete story from a CAD part to a running ECU.

## 2. Reference engine

A generic modern DOHC 16-valve NA inline-4 (not a copy of any production engine):

| Parameter | Value |
|---|---|
| Bore x stroke | 86 x 86 mm (square), 1998 cc |
| Compression ratio | 10.5 : 1 (decision 2 can change this) |
| Conrod length | 145 mm (rod/stroke 1.69) |
| Valvetrain | DOHC, 4 valves/cyl, direct-acting bucket tappets, intake cam phasing (VVT) |
| Firing order | 1-3-4-2, flat-plane crank, 5 main bearings |
| Fuel | Port injection, gasoline (iso-octane / PRF surrogate) |
| Targets (modelled) | ~110 kW at 6500 rpm, ~200 Nm at 4500 rpm, redline 7000 rpm, min BSFC ~245 g/kWh |

All parameters live in `spec/engine.yaml`, the single source of truth (same rule as Talon's `vehicle.yaml`). Every number the project reports is labelled "modelled" because there is no physical engine.

## 3. Headline deliverables

1. **Parametric CAD (build123d)**: block with bores, water jacket and main bulkheads; head with ports, combustion chamber (pent-roof) and cam bearings; crankshaft with counterweights; conrods; pistons with rings and pin; cams and valves; intake manifold and exhaust header. Mass properties come straight from CAD into the dynamics model. STEP, STL and glTF exports.
2. **Kinematics**: exact slider-crank (position, velocity, acceleration to 4th harmonic), cam profile design (polynomial / polydyne lift with jerk limits), valve lift-duration curves, piston-to-valve clearance check across the VVT range.
3. **Thermodynamics and combustion** (crank-angle resolved, 0.1 deg steps):
   - Two-zone model (burned / unburned) with a Wiebe burn rate, then a quasi-dimensional turbulent flame model as an upgrade.
   - Real gas properties and equilibrium burned gas via **Cantera** (GRI-Mech / a gasoline surrogate mechanism).
   - Woschni heat transfer, ring blow-by, residual gas fraction.
   - Knock: Livengood-Wu integral with a Douaud-Eyzat autoignition correlation, giving knock-limited spark advance.
   - Emissions: extended Zeldovich NOx and equilibrium CO in the burned zone.
4. **1D gas dynamics**: finite-volume (HLLC / MUSCL) solver for intake runners, plenum, throttle, exhaust header and collector, coupled to the cylinders through valve flow coefficients. This captures intake tuning and ram effects that give the torque curve its shape. Cross-check against OpenWAM (open-source 1D engine code) if it builds in the container.
5. **Dynamics and NVH**:
   - Gas plus inertia forces, crank torque curves, flywheel sizing for cyclic speed variation.
   - Primary and secondary shaking forces and moments; the I4's 2nd-order imbalance, and an optional Lanchester balance shaft pair as a design study.
   - Main and big-end bearing load polar diagrams; hydrodynamic journal bearing check (mobility method), minimum oil film thickness.
   - Crankshaft torsional vibration: lumped multi-mass model, natural frequencies (Holzer plus eigenvalue check), forced response by order, torsional damper tuning.
   - Valvetrain dynamics: multi-mass valve spring model, surge and float speed, cam/follower Hertz contact stress.
6. **Structures (FEA)**: gmsh meshing of the CAD plus CalculiX (or scikit-fem as fallback). Conrod (buckling, fatigue at peak firing and at redline inertia), crankshaft fillets (fatigue with Goodman), piston crown (thermal plus mechanical with Woschni heat flux), main bearing caps. Each part iterates until its safety factor target holds, same loop as drone-design's frame.
7. **Engine control unit**: C11 firmware, deterministic, no heap.
   - 60-2 crank wheel plus cam sensor decoding, sync and angle-based event scheduling.
   - Speed-density fueling with a VE table generated by the simulator, sequential injection timing, wall-wetting compensation (Aquino / x-tau).
   - Spark map with knock control (retard and recover), closed-loop lambda (PI with short and long term fuel trims).
   - Electronic throttle (DC motor plus return spring plant, PID with friction compensation), idle speed control, VVT phaser control, rev limiter, sensor fault handling with limp mode.
   - A **mean-value engine model** (Hendricks style) for fast control development, with the crank-angle model for detailed validation.
8. **Hardware-in-the-loop on Vortex-V**: the ECU firmware compiled with riscv64 gcc and run in Vortex-V Verilator co-sim in lockstep with the Python plant. Report the cycle budget per engine event at 7000 rpm. An **int8 neural virtual sensor** (knock detection from simulated block accelerometer signal, or cylinder air charge estimation) runs on Xmlacc with the shared `rshift_round` requant rules, bit-exact against its numpy golden model.
9. **Learned control comparison**: a PPO transient fueling or idle controller (reusing NeuroFlight's training scaffold) compared against the classical ECU on lambda excursion during tip-in and tip-out, with a supervisor fallback the same way NeuroFlight uses Simplex.
10. **Vehicle level**: longitudinal vehicle model (gearbox, tyres, aero drag), BSFC map from the simulator, WLTP drive cycle fuel consumption and acceleration (0-100 km/h).
11. **Visualisation and report**: a browser three.js replay of the engine from the glTF export, with pistons, crank and cams animated from the kinematics, and a live P-V diagram, torque trace and valve lift beside it. Published as an Artifact, plus a written technical report.

## 4. Architecture

```
redline-engine/
  spec/engine.yaml        single source of truth (geometry, materials, cams, maps)
  cad/                    build123d parts and assembly -> STEP/STL/glTF + mass properties
  kin/                    slider-crank, cam design, clearance checks
  thermo/                 two-zone cycle, Wiebe + flame model, Cantera props, heat transfer, knock, NOx
  gasdyn/                 1D finite-volume pipes, junctions, valves, throttle
  dyn/                    forces, balance, bearings, torsional, valvetrain dynamics
  fea/                    gmsh + CalculiX runners, load cases, fatigue
  mvem/                   mean-value model for control work
  ecu/                    C11 firmware (portable core + host and RV64 targets)
  hil/                    Python plant <-> ECU bridge, Vortex-V co-sim
  ml/                     int8 virtual sensor, PPO controller
  vehicle/                drive cycle and performance
  viz/                    three.js replay
  tests/                  pytest + C unit tests; analytic and golden fixtures
  docs/adr/               decisions (ADR 0001: units, crank-angle convention, sign conventions)
```

Conventions: SI units, crank angle 0 deg = TDC firing of cylinder 1, positive torque drives the crank forward. Golden fixtures are committed and regenerated only as a reviewed change, as in the other repos.

## 5. Milestones

| M | Scope | Verify command | Exit criteria |
|---|---|---|---|
| M0 | Repo scaffold, spec, ADRs, CI; slider-crank kinematics; ideal Otto cycle | `pytest -q tests/` | Kinematics match closed form to 1e-12; Otto efficiency matches 1 - r^(1-gamma) exactly |
| M1 | Parametric CAD of the full rotating and reciprocating assembly, block and head; mass properties | `python -m cad.build && pytest -q tests/test_cad.py` | Assembly builds with no interference; reciprocating mass and crank inertia within 3% of spec budget |
| M2 | Two-zone thermo with Wiebe, Cantera properties, Woschni, single cylinder | `python -m thermo.run --rpm 3000 --load full` | Matches a published Heywood worked example within stated tolerance; energy balance closes to 0.5% |
| M3 | Cam design + valvetrain kinematics and dynamics; clearance across VVT range | `pytest -q tests/test_m3_valvetrain.py` | No piston-valve contact; no float below 7500 rpm; Hertz stress under material limit |
| M4 | 1D gas dynamics coupled to 4 cylinders; full-load torque curve | `python -m sim.wot --sweep 1000:7000:250` | Shock tube and acoustic resonance tests pass; torque curve shape and peak within 10% of targets; tuned-runner peak where the formula predicts |
| M5 | Dynamics and NVH: balance, bearing loads, torsional, flywheel | `python -m dyn.report` | Secondary force matches analytic 2nd-order formula; torsional modes agree between Holzer and eigen solve within 1% |
| M6 | FEA of conrod, crank, piston, main cap; fatigue | `python -m fea.run --all` | Every part meets its SF target (conrod buckling >= 3, fatigue >= 1.5); mesh convergence shown |
| M7 | Knock, quasi-dimensional flame model, emissions; full speed-load maps (BSFC, VE, MBT/KLSA spark, NOx) | `python -m sim.maps` | Maps generated and committed as fixtures; knock limit moves the right way with CR and intake temperature |
| M8 | ECU firmware on host against the mean-value model, then the crank-angle model | `make -C ecu test && python -m hil.run --target host --cycle tipin` | Lambda within +/-3% in steady state, idle held at 800 +/- 25 rpm, knock control recovers |
| M9 | HIL on Vortex-V + int8 virtual sensor on Xmlacc; PPO controller comparison | `make -C hil cosim` and `python -m ml.eval` | ECU outputs bit-exact vs host build; cycle budget at 7000 rpm reported; int8 sensor bit-exact vs golden |
| M10 | Vehicle drive cycle, three.js replay, report. Stretch: turbo variant (option B), ECU PCB in Zener, intake port CFD | `python -m vehicle.wltp` | WLTP consumption and 0-100 reported; replay published |

Milestones run in this order because each one feeds the next: CAD mass properties feed dynamics, the thermo model feeds the gas dynamics, the maps feed the ECU.

## 6. Links to other projects

| Project | Link |
|---|---|
| **Vortex-V** | ECU firmware runs as RV64 code in Verilator co-sim; int8 virtual sensor uses Xmlacc `dot8`/`macc`/`rqs`. Same HIL pattern as Talon M6, so the bridge code can be shared. |
| **drone-design** | Same build123d CAD workflow and FEA iteration loop; optional ECU board (injector and ignition drivers, VR sensor conditioner) designed in Zener with the same `pcb` toolchain and KiCad 10 Docker image. |
| **NeuroFlight** | PPO training scaffold and Simplex-style supervisor reused for the learned transient fueling controller; same int8 requant rules. |
| **Talon** | Shared conventions (spec YAML as single source of truth, analytic tests first, three.js replay). |
| **MLACC-UVM** | Optional: the virtual sensor's int8 kernel adds a directed test to the Xmlacc bench. |

## 7. Tooling (all cloud)

Python 3.11, numpy, scipy, numba (1D solver speed), Cantera, build123d, gmsh (pip), CalculiX (`apt install calculix-ccx`), matplotlib, gymnasium, torch (CPU). C11 with gcc and riscv64-unknown-elf-gcc; Verilator (both already used by Vortex-V). three.js for the replay. Fable 5.2 only for the hardest pieces (1D solver stability, torsional and FEA debugging, the HIL timing), per the project model policy.

## 8. Risks

- **No real engine to calibrate against.** Mitigation: validate each model against textbook worked examples (Heywood, Stone) and analytic cases, compare headline numbers to published spec-sheet figures for engines of similar size, and label all results as modelled.
- **1D gas dynamics stability** at valve opening and closing (strong gradients, tiny areas). Mitigation: HLLC with positivity limiter, CFL control, validate on Sod shock tube and an organ-pipe resonance case first.
- **Cantera cost** inside a 0.1 deg loop. Mitigation: precompute property and equilibrium tables, interpolate in the loop, check against direct Cantera calls in tests.
- **FEA meshing of CAD** (thin fillets, ports). Mitigation: defeature for meshing, mesh convergence study, scikit-fem fallback if CalculiX is unavailable.
- **HIL cycle budget at 7000 rpm**: a 60-2 wheel gives a tooth every ~143 us. Mitigation: angle-domain scheduling with precomputed tables; fixed-point controller first, neural sensor as the stretch part.
- **Scope**: this is the largest plan so far. Turbo, CFD and the ECU PCB are explicitly outside the core until M10.

## 9. Decisions for the user

Defaults are in bold; work can start on the defaults.

1. **Direction**: **A (NA inline-4 digital twin)**, B (turbo), C (rotary), or D (hybrid).
2. **Reference engine**: **2.0 L NA, 86 x 86 mm, CR 10.5**, or another displacement / layout (for example a 1.6 L or a V6).
3. **ECU HIL target**: **Vortex-V co-sim** (ties into your core), an emulated STM32 via Renode, or both.
4. **Stretch items to include in M10**: turbo variant, ECU PCB in Zener, intake port CFD with OpenFOAM (CPU heavy, slowest in the cloud). Default: **turbo and ECU PCB, no CFD**.
5. **Repo**: create it from your terminal with `gh repo create Ashtonantony28/redline-engine --public` (no README), or pick another name. Until then work lives in the container with a bundle backup in `/mnt/project-files/project-plans/car-engine/`.
