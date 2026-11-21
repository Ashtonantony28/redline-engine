# ADR 0006: Forged 2618-T61 piston with boss pillars, full-height side panels and under-crown oil-jet cooling

Status: Accepted
Date: 2026-11-21

## Context
The first M6 FEA pass failed the piston: Goodman SF 0.44 (target 1.5), static SF 0.71. The pin bosses hung
from a thin rib under the crown, and the crown cantilevered out from that rib. The material was a generic
cast alloy with an ASSUMED temperature derating. Budget-free probes reached only 1.31 at 0.442 kg, so that
M1 architecture could not meet the target within the 0.330 kg piston-assembly budget.

## Decision
- **Material.** The piston is forged AA 2618-T61 (`spec/engine.yaml` `al_2618_t61`). S_ut(T) and S_e(T)
  come from cited tables:
  - Han et al., Adv. Eng. Mater. 2026, doi:10.1002/adem.202503187. Overaged 1000 h at 230 C; Basquin fits
    at RT and 230 C, evaluated at 1e8 cycles.
  - ASM Handbook Vol. 2, 2618-T61 tensile strength after up to 10,000 h at temperature.
  The assumed derating (`spec/fea.yaml` `al_derating`) is retired. Two parts of the tables are
  assumptions, and the engine.yaml comments mark them: the linear interpolation between RT and 230 C, and
  the scaling above 230 C by the ASM tensile ratio.
- **Architecture** (`cad/piston.py`, `spec/cad.yaml` `piston`):
  - Each boss is carried straight up into the crown as a pillar whose width equals the boss diameter, with
    a tangent haunch into the crown (`boss_pillars`, shared by the CAD and the FEA).
  - The side panels run from the skirt bottom to the crown underside. Together with the crown, ring belt,
    skirt and pillars they form a box structure.
  - Compression height 28.5 -> 27.5 mm (top land 5.5 -> 5.0 mm). The oil groove still clears the pin bore.
  - The skirt flat moves from 32.5 to 30 mm, so the pin shortens from 61 to 56 mm.
  - The belt tie and crown bridge options were evaluated and are left off: in this model they made the
    result worse (`fea/design_history.json`).
- **FEA model.** The skirt and side panels are now in the piston FEA, with blends. The first pass removed
  them, which understated the crown support.
- **Budget.** `masses.piston_assembly` 0.330 -> 0.377 kg.
- **Under-crown oil-jet cooling** (piston cooling jets, standard on modern high-output engines). Chosen over
  more mass (capped at about 1.38 with the budget ignored) and over waiting for >230 C fatigue data, because
  the fatigue critical point was the hot crown top (524 K, S_e 39.7 MPa) where the cited 2618 data run out.
  - Hardware (`spec/thermo.yaml` `piston_cooling_jets`, MODELLED design values): one twin-orifice nozzle per
    cylinder, 2 x 2.0 mm orifices, 3.0 bar regulated supply, Cd 0.70 -> 18.6 m/s, 3.5 L/min per orifice
    (28.1 L/min engine), aimed at the crown underside on the thrust axis 22 mm off the piston axis, one jet
    each side of the conrod small end (6.3 mm clearance), check valve open above 2000 rpm. Oil 110 C
    (`spec/fea.yaml` T_oil), viscosity from `spec/dynamics.yaml`, rho / cp / k from Varghese & Agarwal.
  - Heat transfer: J. Stevens and B. W. Webb, J. Heat Transfer 113 (1991) 71-78 (free single-phase round
    liquid jet): Nu_0 = 2.67 Re^0.567 Pr^0.4 (z/d)^-0.0336 (u/d)^-0.237, Nu/Nu_0 = [1 + f^-9]^-1/9 with
    f = 1.13 exp(-0.23 r/d) (2.2 mm row), as used for piston oil-jet cooling by Varghese & Agarwal
    (SAE 2004-28-0061) and Varghese, Goyal & Agarwal (SAE 2005-01-1382). Re 4510 (range 4,000-52,000),
    Pr 113, stagnation h 14.4 kW/(m^2 K). Applied as max(splash, jet) on underside faces facing the jet;
    the splash film elsewhere is frozen at the M6 first-pass calibrated 1550 W/(m^2 K). ASSUMED parts: the
    decaying sign of the radial function (printed without it in the restatement), the cycle-mean relative
    velocity equal to the jet velocity (piston velocity averages out; Varghese & Agarwal), a flat-plate
    correlation on the crown underside, oblique incidence neglected (Sparrow & Lovell: mean h barely changes).
  - No cooling gallery: the jets alone bring the crown top below the 503 K end of the cited fatigue data.
  - Thermal consistency: the crown temperature is now an FEA output. `spec/thermo.yaml` T_piston (Woschni wall
    temperature of the cycle model) is set to the FEA area-mean crown temperature at the max-flux WOT speed:
    530 (assumed) -> 503 -> 505 K, converged (FEA 505.2 K vs 505 K, tolerance 3 K, checked by
    `tests/test_m6_fea.py`).
  - Cost: the jet flow is extra delivery of a variable-displacement pump, P = Q dp / eta = 281 W at
    eta 0.50 (ASSUMED), added to FMEP (`gasdyn/friction.py fmep_total`, `spec/gasdyn.yaml`
    `friction.piston_cooling_jets`): 0.084 bar at 2000 rpm, 0.037 bar at 4500 rpm, 0.024 bar at 7000 rpm,
    zero below 2000 rpm.
- **Structure re-balanced for the cool crown** (medium-mesh probes in `fea/design_history.json`): with the
  crown no longer limiting, the skirt-bottom panel junction (cold, purely mechanical, SF 1.32) and the
  pillar-to-crown haunch governed. Side panels 3.5 -> 3.0 mm, skirt wall 1.8 -> 2.02 mm (mass rebalanced,
  budget unchanged), pillar haunch radius 8 -> 10 mm, and a 4 mm panel-to-skirt-bore blend
  (`spec/fea.yaml` `panel_skirt_fillet`; a drawing requirement, the CAD omits blends as before, about 1 g).

## Consequences
- Modelled piston SF (finest mesh): fatigue 1.57, static 2.41 (first M6 pass 0.44 / 0.71; pillar
  architecture without jets 1.25 / 1.83). Peak stress changes 4.2 % between the two finest meshes. The
  limiting point is now the skirt bottom edge at the side-panel/skirt-bore blend (404 K, purely mechanical:
  sigma_a 41.5 MPa vs S_e 80 MPa). The hot crown top (now 510 K) is at 1.76; crown max 517.5 K, crown mean
  505.2 K. Caveat: the peak there still rises with refinement (95.8 / 99.2 / 103.5 MPa; SF 1.64 / 1.62 / 1.57),
  so the 5 % criterion is met but Richardson extrapolation is not applicable; a fourth mesh level or an edge
  break on the skirt bottom should confirm the margin before release.
- Reciprocating mass rises from 0.474 to 0.521 kg (budget). This raises conrod overrun tension, crank
  loads and the secondary shaking force. The conrod needs 1 mm more flange width to stay at 1.5.
- The piston now depends on the cooling jets: a blocked or failed jet (or running with them shut, e.g. a
  switchable-jet strategy at high load) puts the crown back above the cited data. A jet failure is a
  durability risk to be covered by an oil-pressure diagnosis (ECU, later milestone).
- The cooler piston (T_piston 530 -> 505 K) changes the cycle model slightly (wall heat loss +0.6 %, slightly
  denser charge: M2 3000 rpm WOT net IMEP 11.701 -> 11.716 bar). With the jet FMEP, the M4 WOT peaks move
  from 200.2 Nm / 108.8 kW to 199.7 Nm / 108.7 kW (same speeds, 4500 / 6500 rpm); the largest change is at
  2000 rpm (-0.7 % torque, jet FMEP 0.084 bar). Minimum BSFC of the M7 maps 243.6 -> 245.2 g/kWh. All within
  the M2 / M4 tolerances. Regenerated (their tests failed on the changed model): fixtures/wot/wot_coarse.csv,
  fixtures/maps/maps_m4.csv, mvem/data/thermo_maps.json and, derived from it, ecu/include/ecu/ecu_calib_gen.h.
- Oil pump load rises by about 28 L/min and 0.28 kW above 2000 rpm; the FMEP penalty is in every brake
  number (M4 WOT, M7 maps, turbo variant). The MVEM friction fit (`spec/mvem.yaml`) is not refitted.
- The quarter-model symmetry dictates a twin-orifice nozzle; a single off-axis jet would need a half model.
