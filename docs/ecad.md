# Redline ECU board (stretch goal, PLAN.md 3.7 / M10)

A 4-layer, 120 x 90 mm engine controller for the modelled 2.0 L inline-4. The schematic is in Zener
(`ecad/`), the layout in KiCad 10, and both are checked by scripts. **Nothing here has been built or
measured.** The board is a design that the firmware in `ecu/` could run on. Like every other result in
this repo, its numbers are modelled or come from datasheets.

The **Vortex-V core on an FPGA is still the research target** (ADR 0005). This board is the practical
target: a commercial MCU with automotive-grade drivers that would run the same portable core
(`ecu/include/ecu/*.h`, ADR 0003).

## 1. Verify

```sh
cd ecad
python3 -m pytest -q check_design.py           # Zener build, firmware I/O map, ERC, datasheet limits
layout/run.sh all fab                          # place, pre-route, autoroute, pours, DRC, plots, fab (docker)
python3 -m pytest -q check_layout.py           # outline, DRC clean, ground plane, decoupling, power copper, fab
cd .. && python3 -m pytest -q tests/test_ecad.py   # both of the above; skips without pcb / docker
```

| Check | Result |
|---|---|
| `check_design.py` | 13 passed: build has no warnings, all firmware I/O mapped, ERC clean, every datasheet limit holds |
| KiCad DRC (`layout/run.sh drc`) | **0 errors, 0 warnings, 0 unconnected** (kicad-cli 10.0.8, `--severity-all`) |
| `check_layout.py` | 10 passed |
| Part count | 230 footprints: 86 R, 83 C, 22 U, 13 TP, 8 Q, 6 D, 4 J, 4 H, 1 each of FB, L, SW, Y |

## 2. Toolchain

The toolchain is rebuilt the same way as drone-design:

- `pcb` 0.4.103 is built from source with cargo into `~/.local/bin/pcb`. Always run it as
  `pcb +local ...`, which selects the locally built toolchain.
- KiCad symbols come from a sparse clone of gitlab `kicad-symbols`. `ecad/lib/build_lib.py` copies or
  authors every symbol and footprint into `ecad/lib`, so the build never reaches the network.
- KiCad 10 runs from the Docker image `kicad/kicad:10.0.8`.
  - `pcb layout` needs pcbnew, so it runs inside the image:
    `docker run --rm -u 0 -e HOME=/root -v /root:/root -v /home/user:/home/user -w $PWD kicad/kicad:10.0.8 /root/.local/bin/pcb +local layout --no-open --offline boards/redline_ecu/redline_ecu.zen`
  - The `layout/run.sh` steps also run in the image.

## 3. MCU and pin plan

**STM32H743VIT6 (LQFP-100).** The reasons:

- It is the drone-design MCU, so its symbol, footprint, decoupling and crystal design are reused as-is.
- The 480 MHz Cortex-M7 has a double-precision FPU, which removes the float slow-loop overrun that
  Vortex-V showed in M9.
- **TIM2 is 32-bit.** It is the 1 MHz wrap-safe event timebase of ADR 0003, and crank and cam capture
  on it, so captured ticks and event ticks share one clock.
- TIM1 and TIM8 give 4 + 4 output compares for injection and ignition.
- TIM4 gives the 4 PWM outputs: ETB IN1, ETB IN2, VVT and the lambda heater.
- 14 ADC channels are used, plus FDCAN and USB FS.

The full map from firmware I/O to MCU pin to net to harness pin is `ecad/io_map.yaml`.
`check_design.py` proves every entry on the netlist.

| Function | Pins |
|---|---|
| ADC | PA0 TPS1, PA1 TPS2, PA2 PPS1, PA3 PPS2, PA4 MAP, PA5 5V_S2 mon, PA6 IAT, PA7 CLT, PB0 VBAT, PB1 lambda UA, PC0 lambda UR, PC1 knock, PC4 ETB current, PC5 5V_S1 mon |
| Crank / cam | PA15 TIM2_CH1, PB10 TIM2_CH3 |
| Injectors | PE9 / PE11 / PE13 / PE14, TIM1_CH1-4 |
| Ignition | PC6-PC9, TIM8_CH1-4 |
| ETB / VVT / heater | PD12, PD13 (TIM4_CH1/2); PD14 VVT; PD15 heater |
| ETB control | PD8 nFAULT, PD9 DISABLE, PD11 nSLEEP |
| Fuel pump | PD10 |
| SPI2 (knock, lambda) | PB13 SCK, PB14 MISO, PB15 MOSI; PB12 knock CS, PE10 lambda CS, PE12 knock INT/HOLD; PA8 MCO1 = 8 MHz to the TPIC8101 |
| CAN, USB, debug | PD0 / PD1 FDCAN1; PA11 / PA12 USB FS; PA13 / PA14 SWD, PB3 SWO; PE3 / PE4 LEDs |

The alternate-function numbers are from ST DS12110 Tables 9/10, recalled from memory, because the ST
site was not reachable. **VERIFY** them before fab.

## 4. Circuit blocks (ecad/modules)

### Power (`power.zen`)

The supply path from the harness:

1. `VBAT_IN` arrives on J1.1 and J1.2.
2. SMCJ24CA bidirectional TVS (load dump, 24 V jump-start stand-off).
3. LM74700-Q1 ideal diode with a DMTH6016LK3Q FET (reverse battery to -14 V), giving `VBAT`.
4. Two 47 µF / 50 V electrolytics.

The rails made from `VBAT`:

- **5V:** an LMR36015-Q1 buck at 2.1 MHz with 3.3 µH.
- **3V3:** an LDL1117 LDO from 5V.
- **5V_S1 and 5V_S2:** two TPS7B4250-Q1 trackers for the off-board sensor supplies.
- **USB:** VBUS feeds 5V through a B5819W diode, so the logic powers up from USB for flashing.

The `check()` calls compare the TVS clamp against every part on VBAT and check the buck ripple and
output.

### Actuators (`outputs.zen`, `etb.zen`)

- **Injectors, VVT, fuel-pump relay:** ZXMS6006DGQ self-protected low-side switches. Each has a
  60-70 V active clamp, which sits above the TVS clamp, so the switch does not conduct during a
  clamped load dump. VVT also has an SBR3U60P1 freewheel diode.
- **Ignition:** UCC27524A-Q1 at 5 V through 100 Ω 1206 resistors, for logic-level coil-on-plug coils
  with the igniter in the coil. **VERIFY:** this assumes smart coils.
- **ETB:** DRV8873H-Q1 H-bridge in PWM mode, with ITRIP enabled and IPROPI into 390 Ω read on PC4.

### Sensors (`sensors.zen`, `knock_lambda.zen`)

- **Crank:** LM1815 VR conditioner.
- **Cam:** hall sensor with a 3.3 V pull-up and an RC filter.
- **0-5 V inputs:** 13.3k / 20k dividers.
- **NTC inputs:** 2.49k pull-ups to VDDA.
- **ESD:** TPD2E2U06-Q1 on the ADC side of every series resistor.
- **Knock:** TPIC8101.
- **Wideband lambda:** CJ125 with the rusEFI network.
- **Lambda heater:** low-side DMTH6016 driven by a UCC27524A.
- **Level shifting:** SN74AHCT125-Q1 shifters to 5 V SPI, because the TPIC8101 needs VIH ≥ 0.7·VDD.

### Comms (`comms.zen`)

- TCAN1044AV-Q1 CAN transceiver with split termination and ESD2CAN24-Q1 protection.
- USB-C with 5.1k CC resistors and TPD2E2U06-Q1 protection.
- Tag-Connect TC2030 for SWD.

### Harness

Two Molex Micro-Fit 3.0 2x12 right-angle headers:

- **J1:** battery and ground, actuators, ETB motor, VBAT out to the heater, CAN.
- **J2:** sensor supplies and returns, analog inputs, crank, cam, knock, LSU.

## 5. Layout

KiCad 10, 4 layers on the JLC04161H-7628 stack:

- **L1 (F.Cu):** parts and routing.
- **L2 (In1.Cu):** one solid ground plane.
- **L3 (In2.Cu):** routing.
- **L4 (B.Cu):** routing.

Ground pours fill L1, L3 and L4. All parts are on the top side.

| Rule | Value |
|---|---|
| Clearance (all classes) | 0.15 mm |
| Signal track | 0.2 mm |
| Supply track | 0.3 mm |
| High-current track | 1.0 mm, plus a pour on its own layer up to 3.2 mm wide |
| Vias | 0.6 / 0.3 mm |
| Copper to edge | 0.4 mm |
| Hole clearance | 0.25 mm |

The 0.15 mm clearance on the Power class relies on conformal coating. IPC-2221 B4 (coated) allows
0.13 mm up to 100 V. Uncoated (B2) would need 0.6 mm on the clamp nets (injector, VVT and pump drains,
which reach up to 70 V). The DRV8873 HTSSOP pins are only 0.25 mm apart, so the coating is part of the
design.

**Floor plan:**

- Both headers face out of the bottom edge: J1 on the right, J2 on the left.
- The right half holds the power stage and the actuator drivers:
  - TVS, ideal diode, bulk capacitors, buck and LDO in a column on the right edge;
  - the injector switches with their drain tabs towards J1;
  - gate drivers, VVT and pump switches, the ETB bridge, the heater FET and CAN above them.
- The left half holds the analogue front ends: LM1815, trackers, TPIC8101, CJ125 and level shifters.
- The MCU sits top-centre, with USB-C on the top edge and the Tag-Connect top-left.

**Flow (`ecad/layout/run.sh`)**, adapted from drone-design's `ecad/layout`. The scripts are copied or
adapted; nothing in drone-design was modified.

| Step | Script | What it does |
|---|---|---|
| place | `autoplace.py` + `place.py` | The parts in a table are placed as given. Every other part is placed greedily at the nearest free spot to the pads it connects to (decoupling to the IC its name points at). `place.py` applies the placement and draws the outline, the In1 plane and the hole keep-outs. |
| pre | `preroute.py` | A ground via at every SMD ground pad. Exposed pads get a thermal-via grid. The low-side switch grounds get two vias. |
| auto | `autoroute.py` | drone-design's A* / PathFinder maze router on a 0.1 mm grid, generalised to a W x H board. It routes on L1, L3 and L4, with a search window per net. |
| finish | `route.py` | Ground pours on L1, L3 and L4, pours along the high-current tracks, dangling vias removed, silkscreen clean-up (`silk.py`), zone refill. |
| drc, fab | `kicad-cli` | DRC report; then gerbers, drill, position file and STEP into `ecad/fab/`. |

10 passed_DETAIL

## 6. Open items (VERIFY before fab)

Datasheets that could not be fetched from the build container:

- **ST DS12110:** the alternate-function table.
- **SMCJ24CA:** VC 38.9 V, VBR(min) 26.7 V, VRWM 24 V (Littelfuse site blocked).
- **Bosch CJ125:**
  - the datasheet is not public, so the whole pinout and the external network come from rusEFI;
  - the UB absolute maximum on the 12 V battery net.

Circuit values to check:

- LM1815 peak-detector capacitor (100 nF) against the SNOSBU8F application figures.
- VR amplitude at redline. 66 V is assumed for the 22k input resistor.
- Injector and VVT inductance (4 mH, 6 mH assumed) against the parts actually used.
- SWPA4030S3R3MT saturation current against the LMR36015 high-side current limit.
- Crystal load capacitors (12 pF, as drone-design).

Assumptions:

- Coil-on-plug coils with the igniter in the coil.
- Conformal coating (section 5).
- Sensor supply load of 40 mA per tracker.

Design caveats:

- LDL1117 is not AEC-Q100 qualified. An automotive LDO (for example TPS7B81-Q1) should replace it
  for a vehicle build.
- Single ground. Power returns (injectors, ETB) and sensor returns share the L2 plane. The sensor
  grounds on J2 are separate pins but are not split on the board.
- USB backfeed: with only USB connected, 5V is fed through D2. The buck's output then sees 5V with
  no input. The LMR36015 tolerates this, but **VERIFY** reverse current into VBAT through the buck's
  body diode.
- Current capacity on ETB_M and VBAT depends on the pours (section 5). Review the pour widths against
  the 7.5 A ITRIP before fab.

## 7. Blocked fetches and incidents during the build

- **Datasheet sites unreachable** from the container (proxy refusal): onsemi, Nexperia, Littelfuse,
  ST, Analog Devices, and the Bosch CJ125 document. The values above that depend on them are marked
  VERIFY.
- **Docker Hub rate limit:** the first `docker pull kicad/kicad:10.0.8` returned HTTP 429. The retry
  succeeded.
- **Prebuilt toolchain download:** on the first `pcb build`, the `pcb` launcher fetched a prebuilt
  toolchain (0.4.103) by itself. It was deleted (`~/.local/share/pcb/toolchains/0.4.103` and the
  downloads), and every run since uses the source-built toolchain via `pcb +local`.
