# Redline: technical report (MODELLED)

> **Every number in this report is MODELLED.** Redline is a code-first digital twin of a generic 2.0 L naturally aspirated inline-4; there is no physical engine, dyno or car. Numbers are computed from the code, `spec/*.yaml` and the committed fixtures by `python -m report.build`, which also regenerates the figures under `out/` (git-ignored, so the image links resolve after a build). Built 2026-11-21.

## Headline numbers

| Quantity | Modelled | Target / reference | Source |
|---|---|---|---|
| Peak brake torque | 200.2 N m at 4500 rpm | 200 N m at 4500 rpm | M4 WOT fixture |
| Peak brake power | 108.8 kW at 6500 rpm | 110 kW at 6500 rpm | M4 WOT fixture |
| Minimum BSFC (WOT sweep / M7 maps) | 237 / 244 g/kWh | 245 g/kWh | M4 / M7 fixtures |
| WLTP Class 3b fuel | 7.03 L/100 km iso-octane (6.74 L/100 km E5-equivalent) | 6.0-8.5 L/100 km band for 2.0 NA compacts | M10 vehicle |
| WLTP CO2 | 150 g/km | | M10 vehicle |
| 0-100 km/h | 9.61 s | 7.5-10.5 s band | M10 vehicle |
| 80-120 km/h in 5th | 13.39 s | | M10 vehicle |
| Secondary shaking force at 7000 rpm | 13.29 kN | analytic 12.99 kN (4 m r w^2 lam) | M5 |
| First torsional mode | 560 Hz | Holzer within 7e-16 | M5 |
| ECU crank-tooth handler at 7000 rpm | 19.0 us (13.3 % of a tooth) | 142.9 us tooth period | M9.1 |

## M0: scaffold, spec, kinematics, Otto anchor

Reference engine from `spec/engine.yaml`: 86 x 86 mm, 1998 cc, CR 10.5, rod/stroke 1.69, firing order 1-3-4-2. Exact slider-crank kinematics are checked against closed forms (to 1e-12) and finite differences; the ideal Otto cycle on the real volume trace reproduces 1 - r^(1-gamma) = 0.6096 (gamma 1.4). At the 7000 rpm redline the mean piston speed is 20.1 m/s (peak 32.9 m/s) and piston acceleration spans -1675 to +3055 g.

*Known gaps:* zero wrist-pin offset only.

## M1: parametric CAD and mass properties

| Item | CAD | Budget | Error |
|---|---|---|---|
| Piston assembly | 0.3322 kg | 0.33 kg | 0.68 % |
| Conrod | 0.4836 kg | 0.48 kg | 0.75 % |
| Reciprocating mass | 0.4773 kg | 0.474 kg | 0.70 % |
| Crank inertia | 0.02003 kg m^2 | 0.02 kg m^2 | 0.14 % |
| Clearance volume (CR from CAD) | 52.60 cm^3 (CR 10.497) | 52.58 cm^3 | 0.028 % |

Exports: STEP / STL / glTF in `out/cad/` (`python -m cad.build`); the glTF parts drive the M10 replay.

*Known gaps:* the block is monolithic (no separate main caps), crank fillets are not in the solid, the head is a simplified slab with chambers, seats and port throats, and the camshafts, valves, intake manifold and exhaust header of the plan are not CAD parts (cams exist as kin.cam profiles only).

## M2: two-zone crank-angle thermodynamics (3000 rpm, full load)

Two-zone Wiebe cycle with Cantera property tables, Woschni heat transfer and filling/emptying gas exchange: net IMEP 11.70 bar (gross 11.90), net indicated efficiency 0.404, peak pressure 52.4 bar at 17.5 deg ATDC, CA50 10.9 deg, peak burned-gas temperature 2685 K, residual fraction 0.039, volumetric efficiency 0.908. Closed-cycle energy balance residual -3.0e-05 (exit criterion 0.5 %).

![P-V diagram](../out/thermo/pv_3000_full.png)

*Known gaps:* the stand-alone M2 cycle still uses a sine-squared valve-lift stand-in (M4 uses kin.cam); no ring blow-by model; no wave dynamics (M4 adds them).

## M3: cams, valvetrain kinematics and dynamics

| Valve | Lift | Seat-to-seat | Float (multi-mass) | Min piston-valve clearance | Peak Hertz | Seating speed at redline |
|---|---|---|---|---|---|---|
| intake | 10.50 mm | 329 deg | 8800 rpm | 1.23 mm (VVT 40 deg) | 400 MPa | 0.265 m/s |
| exhaust | 9.50 mm | 331 deg | 10190 rpm | 3.54 mm (VVT 0 deg) | 356 MPa | 0.263 m/s |

Limits: no float below 7500 rpm, Hertz 1200 MPa, clearance 1.00 mm. All met.

![cam_profiles](../out/valvetrain/cam_profiles.png)
![float_sweep](../out/valvetrain/float_sweep.png)

*Known gaps:* the first ramp design bounced at lash take-up (near 4400 rpm per the M3 commit); ramps were slowed and contact damping raised. Spec timings are interpreted at the 1.0 mm check lift, so seat-to-seat events are 89-91 deg longer than the spec durations.

## M4: 1D gas dynamics and full-load curve

HLLC/MUSCL finite-volume pipes (intake duct, throttle, plenum, 4 runners; 4-2-1 exhaust, catalyst, silencer) coupled to four two-zone cylinders. Committed WOT fixture (`fixtures/wot/wot_coarse.csv`, 13 speeds): peak 200.2 N m at 4500 rpm and 108.8 kW at 6500 rpm (targets 200 N m / 110 kW, within 0.1 % / 1.1 %), trapped VE 0.66-1.09, minimum BSFC 237 g/kWh, all points converged: True. Validated on Sod, orifice and quarter-wave anchors; mass conserved to round-off.

![WOT curve](../out/report/wot_coarse.png)

*Known gaps:* catalyst and silencer internals are lumped as friction multipliers in the exhaust pipes; the torque curve has wave-tuning dips at 1500 rpm (122 N m), 2500 rpm (125 N m), 5000 rpm (178 N m) on this grid; spark is the fixed spec Wiebe timing (no knock limit) in the WOT sweep; the planned OpenWAM cross-check has not been done.

## M5: dynamics and NVH (now on REAL pressure traces)

Pressure source: M4 1D gas-dynamics engine (gasdyn.engine, WOT, two-zone Wiebe closed part) at full load; M2 quasi-steady two-zone cycle (thermo.cycle), intake pressure solved for zero brake torque at idle (`python -m dyn.report --pressure m4`; the M5 placeholder stays the default for the fast tests).

* Balance at 7000 rpm: primary forces and all pitch/yaw moments cancel to round-off (max 5e-13 N m); secondary vertical force 13290 N vs 12992 N from 4 m r w^2 lam (exact kinematics add the higher lam terms); a Lanchester pair needs m e = 3.09 kg mm per shaft.
* Torsional: modes 560, 1515, 2364 Hz (Holzer agrees to 7e-16); dominant resonance order 6 at 5605 rpm, twist 0.271 deg bare -> 0.037 deg with the tuned rubber damper.
* Flywheel (spec 0.09 kg m^2): cyclic irregularity C_s 0.039 at idle, 0.081 at 1500 rpm full load (target 0.05); meeting the target at 1500 rpm full load needs 0.161 kg m^2.
* Bearings at 4500 rpm full load: big-end peak 31.9 kN, minimum oil film 2.04 um big end, 3.84 um main 4 (limit ~1 um).
* Bearings at 7000 rpm full load: big-end peak 22.2 kN, minimum oil film 1.71 um big end, 2.68 um main 3 (limit ~1 um).

Change of key numbers when the placeholder Wiebe trace is replaced by the M4 / thermo traces:

| Quantity | Placeholder | M4 / thermo | Change |
|---|---|---|---|
| peak pressure, full load @ torque rpm [bar] | 54.52 | 70.07 | +28.5 % |
| peak pressure, full load @ redline [bar] | 54.52 | 60.06 | +10.2 % |
| peak pressure, idle [bar] | 12.96 | 7.964 | -38.6 % |
| IMEP, full load @ torque rpm [bar] | 12.11 | 15.15 | +25.1 % |
| IMEP, full load @ redline [bar] | 12.11 | 12.64 | +4.3 % |
| IMEP, idle [bar] | 3.045 | 0.5483 | -82.0 % |
| peak crank torque @ torque rpm [N m] | 367.6 | 520.6 | +41.6 % |
| peak crank torque @ redline [N m] | 914.3 | 872.8 | -4.5 % |
| peak rod force @ torque rpm [N] | 2.58e+04 | 3.484e+04 | +35.0 % |
| flywheel dE, idle [J] | 79.64 | 31.31 | -60.7 % |
| C_s idle [-] | 0.09934 | 0.03905 | -60.7 % |
| C_s full load 1500 [-] | 0.1029 | 0.0811 | -21.2 % |
| flywheel needed, idle [kg m^2] | 0.2027 | 0.06499 | -67.9 % |
| flywheel needed, full load 1500 [kg m^2] | 0.2108 | 0.161 | -23.6 % |
| big-end max load 4500 [N] | 2.29e+04 | 3.191e+04 | +39.4 % |
| big-end max load 7000 [N] | 2.194e+04 | 2.216e+04 | +1.0 % |
| big-end h_min 4500 [um] | 2.129 | 2.036 | -4.4 % |
| big-end h_min 7000 [um] | 1.785 | 1.713 | -4.0 % |
| main h_min 4500 (worst main) [um] | 4.466 | 3.842 | -14.0 % |
| main h_min 7000 (worst main) [um] | 2.753 | 2.681 | -2.6 % |
| dominant torsional order twist, bare [deg] | 0.2314 | 0.2712 | +17.2 % |
| dominant torsional order twist, damped [deg] | 0.03152 | 0.03694 | +17.2 % |
| all-order twist peak, bare [deg] | 0.5478 | 0.6421 | +17.2 % |
| all-order twist peak, damped [deg] | 0.3839 | 0.4406 | +14.8 % |
| dominant torsional order [-] | 6 | 6 |  |

![pressure_m4](../out/dyn/m4/pressure_m4.png)
![crank_torque](../out/dyn/m4/crank_torque.png)
![bearing_polar_4500](../out/dyn/m4/bearing_polar_4500.png)
![torsional_forced_response](../out/dyn/m4/torsional_forced_response.png)

*Known gaps:* the M4 full-load traces use the spec Wiebe burn at fixed spark, not the M7 flame model or the knock-limited spark, so peak pressure at the torque peak is likely overstated; the idle trace comes from the quasi-steady M2 cycle; torsional gas-torque orders use one trace (torque-peak speed) for the whole speed sweep; the short-bearing mobility solution assumes rigid shells.

## M6: structural FEA

**In progress** in a parallel work stream (`fea/`, `spec/fea.yaml`: gmsh + CalculiX for conrod, crank, piston and main cap with Goodman fatigue and mesh convergence). Interim results in `out/fea/` are not final and are not summarised here. Exit criterion: conrod buckling SF >= 3, fatigue SF >= 1.5, convergence shown.

## M7: knock, flame model, emissions and speed-load maps

Quasi-dimensional turbulent flame, Livengood-Wu / Douaud-Eyzat knock and extended Zeldovich NOx on the M4 engine. Fixture `fixtures/maps/maps_m4.csv`: 33 of 42 grid points (9 skipped beyond the 0.25 residual dilution limit). Minimum BSFC 244 g/kWh at 3000 rpm (wot, BMEP 8.7 bar); knock limits spark at 6 points (1000/wot, 2000/wot, 2000/0.90bar, 2000/0.80bar, 4000/wot, 4000/0.90bar); NOx 4.6-72.7 g/kWh.

![bsfc_g_per_kwh](../out/report/maps/bsfc_g_per_kwh.png)
![spark_deg](../out/report/maps/spark_deg.png)

*Known gaps:* MBT CA50 lands at 3.4-6.2 deg ATDC, earlier than the usual ~8-10 deg (short burn-up tail of the flame model); the M4 network's high residuals at low speed and part load leave 9 points of the grid empty (covered only by `--quasi-steady`), which is why the M10 drive cycle extrapolates fuel flow below the map.

## M8: C11 ECU on host against the mean-value model

* Idle: mean 800.7 rpm after the accessory step, max deviation 3.2 rpm (exit: 800 +/- 25), minimum 656 rpm during the step.
* Tip-in at 2000 rpm: steady lambda mean 0.9926 (max deviation 0.017; exit +/- 0.03), transient lambda 0.943-1.093, final trims STFT -0.020 / LTFT -0.019.

![tip-in](../out/report/hil_tipin.png)

*Known gaps:* the closed loop runs against the mean-value plant; running it against the crank-angle model is still open. The M8 commit message quotes lambda within 0.3 % steady; the tip-in steady window recomputed here (with the 4 % injector flow error the plant carries) gives the deviation above, inside the exit band.

## M9: ECU on Vortex-V RTL, int8 knock sensor, PPO fuel trim

* Lockstep co-sim: ECU events and slow-loop outputs bit-exact vs the host build (0 event mismatches over 2 windows).
* Frozen-time budget at 7000 rpm on a 50 MHz core (M9.1 firmware; M9 figure in brackets): tooth handler worst 19.0 us [39.3] of a 142.9 us tooth (13.3 %); soft-float slow loop 363 us [934] of 1000 us; CPU utilisation 42.0 % [106.7 %].
* int8 knock MLP on Xmlacc: 54 vectors, 0 mismatches vs the numpy golden model, 3010 cycles per inference; ROC AUC 0.955 (float 0.955, band-pass energy baseline 0.888).
* PPO + Simplex supervisor vs classical over 6 tip-in/tip-out cycles: mean rms lambda deviation 0.0157 -> 0.0139 (-11 %), mean peak deviation 0.0637 -> 0.0666 (+4 %).

*Known gaps:* frozen-time co-sim hides latency and overload: the M9 firmware needed more than one 50 MHz core at 7000 rpm and polled its teeth (M9.1 below measures and fixes this). PPO only beats the classical ECU on rms lambda, not on peak excursion.

## M9.1: the ECU on Vortex-V in real time (ADR 0005)

Plant time advances with core cycles; the frozen-time mode is kept for bit-exactness. All values MODELLED (`make -C hil cosim`, `make -C hil rt-sweep`).

* Design: crank / cam / knock captures are serviced from the Vortex-V machine timer interrupt, armed for the next predicted tooth; the slow loop runs in the background on a snapshot of the decoder and scheduler state; soft-float and sync-path division are replaced by exact, faster target routines. The portable core `ecu/src` is unchanged and every output stays bit-exact against the host build in frozen time.
* Frozen time: crank handler 1430 -> 705 cycles mean; slow loop 44.1k -> 16.4k cycles mean.
* Real time at 50 MHz, M9 firmware vs M9.1 firmware:

| | 7000 rpm (redline_rt), M9 | 7000 rpm (redline_rt), M9.1 | 2000 rpm (tipin), M9 | 2000 rpm (tipin), M9.1 |
|---|---|---|---|---|
| worst tooth-to-event latency | 946.0 us | 18.1 us | 945.9 us | 18.9 us |
| worst tooth-to-serviced latency | 965.1 us | 46.3 us | 958.0 us | 54.3 us |
| spark angle error vs host, max | 34.936 deg | 0.042 deg | 10.580 deg | 0.000 deg |
| any event angle error, max | 38.02 deg | 0.69 deg | 11.19 deg | 0.21 deg |
| worst case (latency as crank angle) | 39.73 deg | 0.76 deg | 11.35 deg | 0.23 deg |
| CPU load | 85.7 % | 48.2 % | 95.1 % | 37.6 % |
| slow-loop steps skipped (overrun) | 35 of 120 | 0 of 120 | 0 of 350 | 0 of 350 |
| captures dropped / loss of sync | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

* Clock sweep on redline_rt (M9.1 firmware), worst case angle error: 20 MHz 1.87 deg (real time lost, 52 slow steps skipped); 60 MHz 0.63 deg; 70 MHz 0.54 deg; 75 MHz 0.51 deg; 77.5 MHz 0.49 deg; 80 MHz 0.47 deg; 100 MHz 0.38 deg. Minimum clock meeting the spec for every event: 77.5 MHz.

*Known gaps:* at 7000 rpm and 50 MHz the 0.5 deg bound is not met for an event that falls within 18 us behind its reference tooth (worst case 0.76 deg, observed 0.69 deg; sparks in the measured window 0.042 deg because they landed mid-window). Closing it at this clock needs the scheduler to arm from the previous tooth (a core change) or a capture/compare unit in the SoC. The M9 firmware drops no tooth only because the modelled capture FIFO holds 32 edges. The worst tooth-to-serviced figures are the start-up poll before the decoder has a period. Vortex-V has not verified FPGA timing at any of these clocks.

## M10: vehicle, drive cycle, replay

Generic compact sedan (`spec/vehicle.yaml`, sources cited there): WLTP test mass 1442 kg, c_d A 0.638 m^2, c_r 0.01, 6-speed manual [3.62, 2.12, 1.45, 1.1, 0.89, 0.74] x final drive 4.06, driveline efficiency 0.92. Fuel flow from a Willans fit to the M7 maps (throttled points reproduced to 0.96 % rms, max 2.4 %), WOT envelope from the M4 fixture, overrun fuel cut, idle at standstill.

| WLTP Class 3b | Distance | L/100 km | g CO2/km |
|---|---|---|---|
| low | 3.095 km | 8.33 | 178 |
| medium | 4.756 km | 6.49 | 138 |
| high | 7.162 km | 6.30 | 134 |
| extra_high | 8.254 km | 7.48 | 160 |
| **combined** | 23.266 km | **7.03** | **150** |

Cycle-mean BSFC 323 g/kWh, positive wheel energy 3.10 kWh, idle share of fuel 2.8 %, trace followed with 0 s of torque shortfall. Acceleration (mass in running order 1350 kg, 81 W/kg): 0-100 km/h 9.61 s, 80-120 km/h in 5th 13.39 s.

![WLTP](../out/vehicle/wltp.png)
![acceleration](../out/vehicle/accel.png)

Replay: `viz/replay/index.html` (self-contained page, three.js 0.160.0 pinned on jsDelivr, built by `python -m viz.build`): the CAD glTF crank, rods and pistons animated from the exact slider-crank at 3000 rpm WOT with live P-V, pressure, crank-torque and valve-lift panels from the M4 engine (brake torque 137.1 N m, peak pressure 49.1 bar).


![replay](../out/viz/replay.png)

*Known gaps:* only 48 % of the WLTP fuel is burned inside the speed/BMEP range the M7 maps sampled; the rest relies on the Willans extrapolation to light load and idle. Hot start only (no cold-start enrichment or catalyst heating), economy shift strategy rather than the GTR 15 Annex 2 algorithm, constant driveline efficiency, no road grade. The replay is not yet published as an Artifact.

## What remains

* **M6 FEA** (in progress): conrod, crankshaft fillets, piston crown and main cap to their safety-factor targets with mesh convergence.
* **M10 stretch** (plan decision 4 default): turbocharged variant (option B: compressor/turbine maps, wastegate, boost control) and the ECU PCB in Zener. Intake-port CFD is not planned.
* Open items carried from earlier milestones: ECU loop against the crank-angle model (M8), event scheduling from the previous tooth or a capture/compare unit in the SoC to reach 0.5 deg for every event at 7000 rpm on a 50 MHz core (M9.1), maps on a knock-limited M4 at low-speed part load (M7), M4 traces with the flame model for the dynamics (M5/M10), publishing the replay, and the GitHub repository.

## Regenerate

```
python -m report.build            # this report and its figures (~2 min)
python -m vehicle.wltp && python -m vehicle.accel
python -m dyn.report --pressure m4
python -m viz.build --check
pytest -q tests/
```
