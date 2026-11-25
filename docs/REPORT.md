# Redline: technical report (MODELLED)

> **Every number in this report is MODELLED.** Redline is a code-first digital twin of a generic 2.0 L naturally aspirated inline-4; there is no physical engine, dyno or car. Numbers are computed from the code, `spec/*.yaml` and the committed fixtures by `python -m report.build`, which also regenerates the figures under `out/` (git-ignored, so the image links resolve after a build). Built 2026-11-25.

## Headline numbers

| Quantity | Modelled | Target / reference | Source |
|---|---|---|---|
| Peak brake torque | 199.7 N m at 4500 rpm | 200 N m at 4500 rpm | M4 WOT fixture |
| Peak brake power | 108.7 kW at 6500 rpm | 110 kW at 6500 rpm | M4 WOT fixture |
| Minimum BSFC (WOT sweep / M7 maps) | 239 / 245 g/kWh | 245 g/kWh | M4 / M7 fixtures |
| WLTP Class 3b fuel | 7.14 L/100 km iso-octane (6.85 L/100 km E5-equivalent) | 6.0-8.5 L/100 km band for 2.0 NA compacts | M10 vehicle |
| WLTP CO2 | 152 g/km | | M10 vehicle |
| 0-100 km/h | 9.62 s | 7.5-10.5 s band | M10 vehicle |
| 80-120 km/h in 5th | 13.47 s | | M10 vehicle |
| Secondary shaking force at 7000 rpm | 14.61 kN | analytic 14.28 kN (4 m r w^2 lam) | M5 |
| First torsional mode | 482 Hz | Holzer within 1e-15 | M5 |
| ECU crank-tooth handler at 7000 rpm | 19.0 us (13.3 % of a tooth) | 142.9 us tooth period | M9.1 |

## M0: scaffold, spec, kinematics, Otto anchor

Reference engine from `spec/engine.yaml`: 86 x 86 mm, 1998 cc, CR 10.5, rod/stroke 1.69, firing order 1-3-4-2. Exact slider-crank kinematics are checked against closed forms (to 1e-12) and finite differences; the ideal Otto cycle on the real volume trace reproduces 1 - r^(1-gamma) = 0.6096 (gamma 1.4). At the 7000 rpm redline the mean piston speed is 20.1 m/s (peak 32.9 m/s) and piston acceleration spans -1675 to +3055 g.

*Known gaps:* zero wrist-pin offset only.

## M1: parametric CAD and mass properties

| Item | CAD | Budget | Error |
|---|---|---|---|
| Piston assembly | 0.3776 kg | 0.377 kg | 0.16 % |
| Conrod | 0.4873 kg | 0.48 kg | 1.53 % |
| Reciprocating mass | 0.5238 kg | 0.521 kg | 0.54 % |
| Crank inertia | 0.02001 kg m^2 | 0.02 kg m^2 | 0.07 % |
| Clearance volume (CR from CAD) | 52.52 cm^3 (CR 10.512) | 52.58 cm^3 | -0.131 % |

| Flywheel + ring gear + bolts inertia | 0.0901 kg m^2 | 0.09 kg m^2 | 0.09 % |

Exports: STEP / STL / glTF in `out/cad/` (`python -m cad.build`, absolute-deflection tessellation from cad/tess.py). Clash, clearance, support and tessellation audit: `python -m cad.audit --tag final` (docs/CAD-AUDIT.md, renders in docs/renders/).

*Known gaps:* crank fillets are not in the CAD solid (FEA only); threads are line-to-line cylinders; the head has no water jacket, oil drains or cam cover; springs are placed at installed length and checked at full lift separately. The M5 bearing model now uses the CAD journal diameters (main 52 mm, big end 45 mm; they were 56/48 mm), checked by tests/test_m5_dynamics.py.

## M2: two-zone crank-angle thermodynamics (3000 rpm, full load)

Two-zone Wiebe cycle with Cantera property tables, Woschni heat transfer and filling/emptying gas exchange: net IMEP 11.72 bar (gross 11.92), net indicated efficiency 0.404, peak pressure 52.5 bar at 17.5 deg ATDC, CA50 10.9 deg, peak burned-gas temperature 2684 K, residual fraction 0.039, volumetric efficiency 0.909. Closed-cycle energy balance residual -3.0e-05 (exit criterion 0.5 %).

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

HLLC/MUSCL finite-volume pipes (intake duct, throttle, plenum, 4 runners; 4-2-1 exhaust, catalyst, silencer) coupled to four two-zone cylinders. Committed WOT fixture (`fixtures/wot/wot_coarse.csv`, 13 speeds): peak 199.7 N m at 4500 rpm and 108.7 kW at 6500 rpm (targets 200 N m / 110 kW, within 0.1 % / 1.2 %), trapped VE 0.66-1.09, minimum BSFC 239 g/kWh, all points converged: True. Validated on Sod, orifice and quarter-wave anchors; mass conserved to round-off.

![WOT curve](../out/report/wot_coarse.png)

*Known gaps:* catalyst and silencer internals are lumped as friction multipliers in the exhaust pipes; the torque curve has wave-tuning dips at 1500 rpm (122 N m), 2500 rpm (125 N m), 5000 rpm (178 N m) on this grid; spark is the fixed spec Wiebe timing (no knock limit) in the WOT sweep; the planned OpenWAM cross-check has not been done.

## M5: dynamics and NVH (now on REAL pressure traces)

Pressure source: M4 1D gas-dynamics engine (gasdyn.engine, WOT, two-zone Wiebe closed part) at full load; M2 quasi-steady two-zone cycle (thermo.cycle), intake pressure solved for zero brake torque at idle (`python -m dyn.report --pressure m4`; the M5 placeholder stays the default for the fast tests).

* Balance at 7000 rpm: primary forces and all pitch/yaw moments cancel to round-off (max 5e-13 N m); secondary vertical force 14607 N vs 14280 N from 4 m r w^2 lam (exact kinematics add the higher lam terms); a Lanchester pair needs m e = 3.40 kg mm per shaft.
* Torsional: modes 482, 1303, 2032 Hz (Holzer agrees to 1e-15); dominant resonance order 6 at 4822 rpm, twist 0.355 deg bare -> 0.049 deg with the tuned rubber damper.
* Flywheel (spec 0.09 kg m^2): cyclic irregularity C_s 0.038 at idle, 0.080 at 1500 rpm full load (target 0.05); meeting the target at 1500 rpm full load needs 0.158 kg m^2.
* Bearings at 4500 rpm full load: big-end peak 31.4 kN, minimum oil film 1.93 um big end, 3.68 um main 2 (limit ~1 um).
* Bearings at 7000 rpm full load: big-end peak 23.6 kN, minimum oil film 1.61 um big end, 2.48 um main 3 (limit ~1 um).

Change of key numbers when the placeholder Wiebe trace is replaced by the M4 / thermo traces:

| Quantity | Placeholder | M4 / thermo | Change |
|---|---|---|---|
| peak pressure, full load @ torque rpm [bar] | 54.52 | 70.1 | +28.6 % |
| peak pressure, full load @ redline [bar] | 54.52 | 60.07 | +10.2 % |
| peak pressure, idle [bar] | 12.96 | 7.963 | -38.6 % |
| IMEP, full load @ torque rpm [bar] | 12.11 | 15.15 | +25.1 % |
| IMEP, full load @ redline [bar] | 12.11 | 12.64 | +4.4 % |
| IMEP, idle [bar] | 3.045 | 0.5483 | -82.0 % |
| peak crank torque @ torque rpm [N m] | 404.9 | 489.1 | +20.8 % |
| peak crank torque @ redline [N m] | 1008 | 966 | -4.1 % |
| peak rod force @ torque rpm [N] | 2.526e+04 | 3.431e+04 | +35.8 % |
| flywheel dE, idle [J] | 78.5 | 30.31 | -61.4 % |
| C_s idle [-] | 0.09776 | 0.03774 | -61.4 % |
| C_s full load 1500 [-] | 0.1012 | 0.07956 | -21.4 % |
| flywheel needed, idle [kg m^2] | 0.1993 | 0.06195 | -68.9 % |
| flywheel needed, full load 1500 [kg m^2] | 0.2073 | 0.1576 | -23.9 % |
| big-end max load 4500 [N] | 2.237e+04 | 3.139e+04 | +40.4 % |
| big-end max load 7000 [N] | 2.335e+04 | 2.357e+04 | +0.9 % |
| big-end h_min 4500 [um] | 2.012 | 1.932 | -4.0 % |
| big-end h_min 7000 [um] | 1.675 | 1.611 | -3.9 % |
| main h_min 4500 (worst main) [um] | 4.275 | 3.685 | -13.8 % |
| main h_min 7000 (worst main) [um] | 2.54 | 2.482 | -2.3 % |
| dominant torsional order twist, bare [deg] | 0.3023 | 0.3547 | +17.3 % |
| dominant torsional order twist, damped [deg] | 0.04216 | 0.04946 | +17.3 % |
| all-order twist peak, bare [deg] | 0.8684 | 1.049 | +20.8 % |
| all-order twist peak, damped [deg] | 0.5672 | 0.6437 | +13.5 % |
| dominant torsional order [-] | 6 | 6 |  |

![pressure_m4](../out/dyn/m4/pressure_m4.png)
![crank_torque](../out/dyn/m4/crank_torque.png)
![bearing_polar_4500](../out/dyn/m4/bearing_polar_4500.png)
![torsional_forced_response](../out/dyn/m4/torsional_forced_response.png)

*Known gaps:* the M4 full-load traces use the spec Wiebe burn at fixed spark, not the M7 flame model or the knock-limited spark, so peak pressure at the torque peak is likely overstated; the idle trace comes from the quasi-steady M2 cycle; torsional gas-torque orders use one trace (torque-peak speed) for the whole speed sweep; the short-bearing mobility solution assumes rigid shells.

## M6: structural FEA

gmsh second-order tets + CalculiX, Goodman fatigue, three mesh levels per part (ADR 0006 for the piston). Loads come from the M4/M2 pressure traces and the M5 inertia model. Finest-mesh results:

| Part | Safety factors | Targets | Peak-stress change, finest two meshes | Result |
|---|---|---|---|---|
| conrod | buckling 5.93, fatigue 1.53 | fatigue 1.5, buckling 3.0 | 0.72 % | PASS |
| crankshaft | fatigue 2.58 | fatigue 1.5 | 0.48 % | PASS |
| piston | fatigue 1.57, static 2.41 | fatigue 1.5, static 1.5 | 4.19 % | PASS |
| main_cap | fatigue 5.55, bolt_fatigue 4.98 | fatigue 1.5, bolt_fatigue 1.5 | 2.89 % | PASS |

The M1 piston (pin boss hung from a thin rib) could not reach fatigue SF 1.5 in any mass budget. It became a forged 2618-T61 pillar-boss design at 0.377 kg (budget raised from 0.330 kg) with under-crown oil jets (Stevens & Webb impingement correlation); the converged 505 K piston temperature is fed back into the cycle model and the jet pump cost into FMEP. Conrod flange widened to keep fatigue SF above 1.5 under the higher reciprocating mass.

*Known gaps:* the piston peak stress still rises with refinement (no Richardson extrapolation), so a fourth mesh level should confirm SF 1.57; skirt side thrust, piston inertia and pin ovalisation are not modelled; 2618 fatigue above 230 C is extrapolated from tensile-strength ratios; a blocked oil jet would put the crown back above the cited data (oil-pressure diagnosis needed).

## M7: knock, flame model, emissions and speed-load maps

Quasi-dimensional turbulent flame, Livengood-Wu / Douaud-Eyzat knock and extended Zeldovich NOx on the M4 engine. Fixture `fixtures/maps/maps_m4.csv`: 33 of 42 grid points (9 skipped beyond the 0.25 residual dilution limit). Minimum BSFC 245 g/kWh at 3000 rpm (wot, BMEP 8.6 bar); knock limits spark at 6 points (1000/wot, 2000/wot, 2000/0.90bar, 2000/0.80bar, 4000/wot, 4000/0.90bar); NOx 4.8-73.8 g/kWh.

![bsfc_g_per_kwh](../out/report/maps/bsfc_g_per_kwh.png)
![spark_deg](../out/report/maps/spark_deg.png)

*Known gaps:* MBT CA50 lands at 3.4-6.3 deg ATDC, earlier than the usual ~8-10 deg (short burn-up tail of the flame model); the M4 network's high residuals at low speed and part load leave 9 points of the grid empty (covered only by `--quasi-steady`), which is why the M10 drive cycle extrapolates fuel flow below the map.

## M8: C11 ECU on host against the mean-value model

* Idle: mean 800.7 rpm after the accessory step, max deviation 3.7 rpm (exit: 800 +/- 25), minimum 648 rpm during the step.
* Tip-in at 2000 rpm: steady lambda mean 0.9928 (max deviation 0.017; exit +/- 0.03), transient lambda 0.943-1.091, final trims STFT -0.020 / LTFT -0.019.

![tip-in](../out/report/hil_tipin.png)

*Known gaps:* the closed loop runs against the mean-value plant; running it against the crank-angle model is still open. The M8 commit message quotes lambda within 0.3 % steady; the tip-in steady window recomputed here (with the 4 % injector flow error the plant carries) gives the deviation above, inside the exit band.

## M9: ECU on Vortex-V RTL, int8 knock sensor, PPO fuel trim

* Lockstep co-sim: ECU events and slow-loop outputs bit-exact vs the host build (0 event mismatches over 2 windows).
* Frozen-time budget at 7000 rpm on a 50 MHz core (M9.1 firmware; M9 figure in brackets): tooth handler worst 19.0 us [39.3] of a 142.9 us tooth (13.3 %); soft-float slow loop 363 us [934] of 1000 us; CPU utilisation 42.0 % [106.7 %].
* int8 knock MLP on Xmlacc: 54 vectors, 0 mismatches vs the numpy golden model, 3007 cycles per inference; ROC AUC 0.955 (float 0.955, band-pass energy baseline 0.888).
* PPO + Simplex supervisor vs classical over 6 tip-in/tip-out cycles: mean rms lambda deviation 0.0157 -> 0.0139 (-11 %), mean peak deviation 0.0634 -> 0.0667 (+5 %).

*Known gaps:* frozen-time co-sim hides latency and overload: the M9 firmware needed more than one 50 MHz core at 7000 rpm and polled its teeth (M9.1 below measures and fixes this). PPO only beats the classical ECU on rms lambda, not on peak excursion.

## M9.1: the ECU on Vortex-V in real time (ADR 0005)

Plant time advances with core cycles; the frozen-time mode is kept for bit-exactness. All values MODELLED (`make -C hil cosim`, `make -C hil rt-sweep`).

* Design: crank / cam / knock captures are serviced from the Vortex-V machine timer interrupt, armed for the next predicted tooth; the slow loop runs in the background on a snapshot of the decoder and scheduler state; soft-float and sync-path division are replaced by exact, faster target routines. The portable core `ecu/src` is unchanged and every output stays bit-exact against the host build in frozen time.
* Frozen time: crank handler 1430 -> 705 cycles mean; slow loop 44.1k -> 16.4k cycles mean.
* Real time at 50 MHz, M9 firmware vs M9.1 firmware:

| | 7000 rpm (redline_rt), M9 | 7000 rpm (redline_rt), M9.1 | 2000 rpm (tipin), M9 | 2000 rpm (tipin), M9.1 |
|---|---|---|---|---|
| worst tooth-to-event latency | 940.3 us | 18.0 us | 948.6 us | 19.1 us |
| worst tooth-to-serviced latency | 967.6 us | 46.8 us | 961.3 us | 55.5 us |
| spark angle error vs host, max | 35.121 deg | 0.042 deg | 10.606 deg | 0.000 deg |
| any event angle error, max | 38.16 deg | 0.70 deg | 11.13 deg | 0.21 deg |
| worst case (latency as crank angle) | 39.49 deg | 0.75 deg | 11.38 deg | 0.23 deg |
| CPU load | 85.7 % | 47.9 % | 95.3 % | 37.7 % |
| slow-loop steps skipped (overrun) | 35 of 120 | 0 of 120 | 0 of 350 | 0 of 350 |
| captures dropped / loss of sync | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

* Clock sweep on redline_rt (M9.1 firmware), worst case angle error: 20 MHz 1.87 deg (real time lost, 52 slow steps skipped); 60 MHz 0.63 deg; 70 MHz 0.54 deg; 75 MHz 0.51 deg; 77.5 MHz 0.49 deg; 80 MHz 0.47 deg; 100 MHz 0.38 deg. Minimum clock meeting the spec for every event: 77.5 MHz.

*Known gaps:* at 7000 rpm and 50 MHz the 0.5 deg bound is not met for an event that falls within 18 us behind its reference tooth (worst case 0.75 deg, observed 0.70 deg; sparks in the measured window 0.042 deg because they landed mid-window). Closing it at this clock needs the scheduler to arm from the previous tooth (a core change) or a capture/compare unit in the SoC. The M9 firmware drops no tooth only because the modelled capture FIFO holds 32 edges. The worst tooth-to-serviced figures are the start-up poll before the decoder has a period. Vortex-V has not verified FPGA timing at any of these clocks.

## M9.2: the ECU on the Vortex-V capture/compare timer (ADR 0007)

The SoC is built with Vortex-V's capture/compare unit (its ADR 0017, `CAPCMP=1`): the plant drives the crank and cam pins and acts on the eight compare outputs. All values MODELLED (`make -C hil cosim`, `make -C hil cc-sweep`).

* Design: crank / cam edges are timestamped by the RTL and serviced from the machine external interrupt (no timer wake, prediction or spin); every event is one compare slot (injector pulse, coil set, coil clear); the HAL counts free slots and holds an event back rather than let the unit drop it; knock windows stay in the mailbox and are merged by tick; the slow-loop commit is atomic against captures. The portable core is unchanged and frozen-time outputs stay bit-exact against the host build.
* Real time at 50 MHz, M9.1 (played mailbox peripheral) vs M9.2 (RTL unit):

| | 7000 rpm (redline_rt), M9.1 | 7000 rpm (redline_rt), M9.2 | 2000 rpm (tipin), M9.1 | 2000 rpm (tipin), M9.2 |
|---|---|---|---|---|
| worst tooth-to-serviced latency | 46.8 us | 20.2 us | 55.5 us | 20.8 us |
| worst tooth-to-event latency | 18.0 us | 19.1 us | 19.1 us | 20.2 us |
| pin edge to interrupt taken, max | - | 32 cycles | - | 21 cycles |
| spark angle error vs host, max | 0.042 deg | 0.024 deg | 0.000 deg | 0.005 deg |
| any event angle error, max | 0.70 deg | 0.72 deg | 0.21 deg | 0.22 deg |
| worst case (latency as crank angle) | 0.75 deg | 0.80 deg | 0.23 deg | 0.24 deg |
| events pushed after their tick | - | 9 of 74 | - | 4 of 60 |
| CPU load | 47.9 % | 48.0 % | 37.7 % | 38.4 % |
| captures lost / CAP_OVF / loss of sync | 0 / - / 0 | 0 / 0 / 0 | 0 / - / 0 | 0 / 0 / 0 |
| max pending compares | - | 2 of 8 | - | 2 of 8 |
| compare pushes dropped / held back | - | 0 / 0 | - | 0 / 0 |

* Clock sweep on redline_rt (M9.2), worst case angle error: 10 MHz 4.158 deg (real time lost, 102 slow steps skipped); 20 MHz 1.953 deg (real time lost, 52 slow steps skipped); 25 MHz 1.614 deg (real time lost, 35 slow steps skipped); 50 MHz 0.802 deg; 55 MHz 0.716 deg; 78 MHz 0.507 deg; 79 MHz 0.501 deg; 80 MHz 0.495 deg; 81 MHz 0.490 deg; 84 MHz 0.472 deg; 89 MHz 0.444 deg; 100 MHz 0.395 deg. Minimum whole-MHz clock meeting the spec for every event: 80 MHz.
* tipin (2000 rpm), worst case: 20 MHz 0.61 deg (real time lost); 25 MHz 0.49 deg; 50 MHz 0.24 deg.

*Known gaps:* the unit fires every compare armed before its tick exactly on that tick, but at 7000 rpm and 50 MHz the scheduler still arms events right behind a tooth after their tick (9 of 74 pushed late; worst 0.72 deg, worst case 0.80 deg against 0.5 deg). The M9.1 gap was arming latency, not output jitter, so a hardware compare alone does not close it; ADR 0007 proposes a scheduler arm-ahead window (a core change with a new host fixture) and rejects HAL look-ahead pre-arming, which was measured. Knock windows still use the mailbox. The capture/compare pins are not on the Arty board yet, and Vortex-V has not verified FPGA timing at any of these clocks.

## M9.3: scheduler arm-ahead (ADR 0008)

*No M9.3 runs in out/hil/realtime.json: run `make -C hil cosim` (needs Verilator and the Vortex-V tree).*

## M10: vehicle, drive cycle, replay

Generic compact sedan (`spec/vehicle.yaml`, sources cited there): WLTP test mass 1442 kg, c_d A 0.638 m^2, c_r 0.01, 6-speed manual [3.62, 2.12, 1.45, 1.1, 0.89, 0.74] x final drive 4.06, driveline efficiency 0.92. Fuel flow from a Willans fit to the M7 maps (throttled points reproduced to 0.95 % rms, max 2.4 %), WOT envelope from the M4 fixture, overrun fuel cut, idle at standstill.

| WLTP Class 3b | Distance | L/100 km | g CO2/km |
|---|---|---|---|
| low | 3.095 km | 8.55 | 182 |
| medium | 4.756 km | 6.61 | 141 |
| high | 7.162 km | 6.40 | 137 |
| extra_high | 8.254 km | 7.55 | 161 |
| **combined** | 23.266 km | **7.14** | **152** |

Cycle-mean BSFC 328 g/kWh, positive wheel energy 3.10 kWh, idle share of fuel 3.0 %, trace followed with 0 s of torque shortfall. Acceleration (mass in running order 1350 kg, 80 W/kg): 0-100 km/h 9.62 s, 80-120 km/h in 5th 13.47 s.

![WLTP](../out/vehicle/wltp.png)
![acceleration](../out/vehicle/accel.png)

Replay: `viz/replay/index.html` (self-contained page, three.js 0.160.0 pinned on jsDelivr, built by `python -m viz.build`): the CAD glTF crank, rods and pistons animated from the exact slider-crank at 3000 rpm WOT with live P-V, pressure, crank-torque and valve-lift panels from the M4 engine (brake torque 136.4 N m, peak pressure 49.1 bar).


![replay](../out/viz/replay.png)

*Known gaps:* only 49 % of the WLTP fuel is burned inside the speed/BMEP range the M7 maps sampled; the rest relies on the Willans extrapolation to light load and idle. Hot start only (no cold-start enrichment or catalyst heating), economy shift strategy rather than the GTR 15 Annex 2 algorithm, constant driveline efficiency, no road grade.

## What remains

* **Stretch, done:** turbocharged I4-2.0T (`turbo/`, ADR 0004; modelled 275 Nm peak, 165 kW, 2.1 s to 90 % boost at 2000 rpm, synthetic compressor/turbine maps, turbine inlet ~1320 K without enrichment) and the STM32H743 ECU board in Zener (`ecad/`, 230 parts, KiCad DRC 0/0/0; conformal-coating clearance, ETB copper and datasheet VERIFY items open, see docs/ecad.md). Intake-port CFD is not planned.
* Open items: ECU loop against the crank-angle model (M8); the M9.3 scheduler arm-ahead (ADR 0008) meets the 0.5 deg event spec at 7000 rpm and 50 MHz on Vortex-V (minimum clock now 26 MHz, set by CPU load), but the capture/compare pins are not on the Arty board and FPGA timing is unverified (M9.2 / M9.3); late catch-up remains for targets that jump behind a tooth by more than the guard (M9.3); maps at low-speed part load where M4 residuals are high (M7); M4 traces with the flame model for the dynamics (M5/M10); a fourth piston mesh level (M6); and the GitHub repository.

## Regenerate

```
python -m report.build            # this report and its figures (~2 min)
python -m vehicle.wltp && python -m vehicle.accel
python -m dyn.report --pressure m4
python -m viz.build --check
pytest -q tests/
```
