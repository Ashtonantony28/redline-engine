# Redline CAD audit (MODELLED geometry)

All numbers describe the parametric CAD model; there is no physical engine. Regenerate with `python -m cad.audit --tag final` (writes out/cad/audit_final.{json,md}) and `python -m cad.audit --docs`.

## Method

* **Static clash** (crank 0 deg, VVT parked): every pair of the 349 placed part instances (60726 pairs) is screened with exact (optimal) bounding boxes; the 1478 pairs whose boxes come within the contact tolerance get an exact BRep boolean common volume and an exact minimum distance (BRepExtrema). Overlap > 0.01 mm^3 = overlap, distance <= 0.05 mm = contact. Every contact/overlap must be listed in spec/cad_fits.yaml with its build reason; anything else is a defect.
* **Moving sweep**: every pair with at least one moving part whose swept boxes come within 5 mm is swept. Crank-driven pairs: 2 deg grid over the full period (720 deg when a cam or valve is involved, at VVT parked and at full intake advance when only one side is phased), local refinement of the minimum down to 0.01 deg, then a branch and bound with a rigorous Lipschitz bound (largest surface speed of the bodies) certifies a lower bound. Valve parts against static parts are swept over the valve lift (0.5 mm grid, refined) instead of the crank angle. Running-contact pairs (bearings, sliding fits, cam on bucket) are checked for penetration (common volume) instead of clearance.
* **Support graph**: union-find over the contact/overlap pairs, rooted at the block; a part not connected is floating (and connections through unlisted contacts are reported separately).
* **Tessellation**: every source solid is meshed with the export settings and the chordal deviation of every triangle from the exact surface (exact point projection of the worst samples) and the largest angle between exact vertex normals inside a triangle are measured. OCCT's mesher alone does not honour the angular deflection across tori and swept B-spline faces, so cad/tess.py refines every curved face afterwards on the exact surface (red-green splits at surface midpoints) until the chord and facet-angle targets hold.

## Baseline (before the polish) vs final

| Check | Baseline | Final |
|---|---|---|
| instances | 35 | 349 |
| bbox candidates | 89 | 1478 |
| overlaps | 0 | 0 |
| unexplained overlaps | 0 | 0 |
| contacts | 45 | 580 |
| unexplained contacts | 1 | 0 |
| floating | 1 | 0 |
| moving pairs | 129 | 1140 |
| moving clashes | 1 | 0 |
| tessellation fail | 1 | 0 |
| worst chord mm | 0.0253 | 0.0193 |
| worst facet angle rad | 0.0914 | 0.0800 |

Baseline defects: block / crankshaft contact at 0.000 mm. Baseline floating: head; held only through defect contacts: crankshaft.

## Final summary

| Check | Count |
|---|---|
| instances | 349 |
| pairs | 60726 |
| bbox candidates | 1478 |
| overlaps | 0 |
| overlaps allowed | 0 |
| contacts | 580 |
| contacts allowed | 580 |
| unexplained overlaps | 0 |
| unexplained contacts | 0 |
| missing required contacts | 0 |
| floating | 0 |
| floating without unlisted contacts | 0 |
| moving pairs | 1140 |
| moving clashes | 0 |
| tessellation fail | 0 |

## Static part pairs (every pair with touching bounding boxes, aggregated by part kind)

All other 59248 pairs have exact bounding boxes more than 0.05 mm apart (clear).

| Part kinds | State | Pairs | Min distance mm | Max overlap mm^3 | Allowed | Reason |
|---|---|---|---|---|---|---|
| big_end_shell_lower / big_end_shell_upper | contact | 4 | 0.000 | 0.000 | yes | shell joint faces at the split |
| big_end_shell_lower / conrod | contact | 4 | 0.000 | 0.000 | yes | lower shell joint face meets the rod at the split plane |
| big_end_shell_lower / conrod_cap | contact | 4 | 0.000 | 0.000 | yes | big-end shell in the cap bore |
| big_end_shell_lower / crankshaft | contact | 4 | 0.022 | 0.000 | yes | big-end running clearance 22 um radial (spec/dynamics.yaml bearings.big_end) = oil film |
| big_end_shell_upper / conrod | contact | 4 | 0.000 | 0.000 | yes | big-end shell in the rod bore |
| big_end_shell_upper / conrod_cap | contact | 4 | 0.000 | 0.000 | yes | upper shell joint face meets the cap at the split plane |
| big_end_shell_upper / crankshaft | contact | 4 | 0.022 | 0.000 | yes | big-end running clearance 22 um radial (spec/dynamics.yaml bearings.big_end) = oil film |
| block / bracket_bolt_block | contact | 2 | 0.000 | 0.000 | yes | stay bolt thread engagement in the block rear wall |
| block / head_bolt | contact | 10 | 0.000 | 0.000 | yes | head bolt thread engagement in the block |
| block / head_gasket | contact | 1 | 0.000 | 0.000 | yes | gasket on the deck |
| block / idler | contact | 1 | 0.000 | 0.000 | yes | idler boss bolted to the block front face |
| block / intake_bracket | contact | 1 | 0.000 | 0.000 | yes | manifold stay plate on the block rear face |
| block / liner | contact | 4 | 0.000 | 0.000 | yes | cast-in liner (line-to-line in the barrel) |
| block / main_cap | contact | 5 | 0.000 | 0.000 | yes | cap on the bulkhead joint face |
| block / main_cap_bolt | contact | 10 | 0.000 | 0.000 | yes | cap bolt thread engagement in the bulkhead |
| block / main_shell_lower | contact | 5 | 0.000 | 0.000 | yes | lower shell joint face meets the saddle face at the split |
| block / main_shell_upper | contact | 5 | 0.000 | 0.000 | yes | upper main shell in the block saddle |
| block / oil_pan | contact | 1 | 0.000 | 0.000 | yes | pan flange on the block rail |
| block / pan_bolt | contact | 16 | 0.000 | 0.000 | yes | pan bolt thread engagement in the rail |
| block / tensioner | contact | 1 | 0.000 | 0.000 | yes | tensioner boss bolted to the block front face |
| bracket_bolt_block / intake_bracket | contact | 2 | 0.000 | 0.000 | yes | stay bolt head on the plate |
| bracket_bolt_pad / intake_bracket | contact | 2 | 0.000 | 0.000 | yes | stay bolt head under the tab |
| bracket_bolt_pad / intake_manifold | contact | 2 | 0.000 | 0.000 | yes | stay bolt thread engagement in the plenum pad |
| bucket_exhaust / camshaft_exhaust | contact | 4 | 0.000 | 0.000 | yes | lobe on the bucket crown while the valve lifts (lash on the base circle) |
| bucket_exhaust / head | contact | 8 | 0.030 | 0.000 | yes | bucket in its bore |
| bucket_exhaust / shim_exhaust | contact | 8 | 0.000 | 0.000 | yes | shim in the bucket crown pocket |
| bucket_intake / camshaft_intake | contact | 4 | 0.000 | 0.000 | yes | lobe on the bucket crown while the valve lifts (lash on the base circle) |
| bucket_intake / head | contact | 8 | 0.030 | 0.000 | yes | bucket in its bore |
| bucket_intake / shim_intake | contact | 8 | 0.000 | 0.000 | yes | shim in the bucket crown pocket |
| cam_cap / cam_cap_bolt | contact | 20 | 0.000 | 0.000 | yes | cap bolt head on the cap |
| cam_cap / camshaft_exhaust | contact | 5 | 0.040 | 0.000 | yes | cam journals under the caps |
| cam_cap / camshaft_intake | contact | 5 | 0.040 | 0.000 | yes | cam journals under the caps |
| cam_cap / head | contact | 10 | 0.000 | 0.000 | yes | cam cap on the journal pedestal joint face |
| cam_cap_bolt / head | contact | 20 | 0.000 | 0.000 | yes | cap bolt thread engagement in the pedestal |
| cam_pulley_exhaust / camshaft_exhaust | contact | 1 | 0.000 | 0.000 | yes | cam pulley on the cam nose |
| cam_pulley_exhaust / timing_belt | contact | 1 | 0.000 | 0.000 | yes | belt on the cam pulley rims |
| cam_pulley_intake / camshaft_intake | contact | 1 | 0.000 | 0.000 | yes | cam pulley on the cam nose |
| cam_pulley_intake / timing_belt | contact | 1 | 0.000 | 0.000 | yes | belt on the cam pulley rims |
| camshaft_exhaust / head | contact | 1 | 0.040 | 0.000 | yes | cam journals in the head saddles |
| camshaft_intake / head | contact | 1 | 0.040 | 0.000 | yes | cam journals in the head saddles |
| clip / piston | contact | 8 | 0.000 | 0.000 | yes | circlip seated in its groove in the pin bore |
| collet_exhaust / retainer_exhaust | contact | 16 | 0.000 | 0.000 | yes | collet cone in the retainer |
| collet_exhaust / valve_exhaust | contact | 16 | 0.000 | 0.000 | yes | collet bead in the stem groove |
| collet_intake / retainer_intake | contact | 16 | 0.000 | 0.000 | yes | collet cone in the retainer |
| collet_intake / valve_intake | contact | 16 | 0.000 | 0.000 | yes | collet bead in the stem groove |
| conrod / conrod_bolt_a | contact | 4 | 0.000 | 0.000 | yes | rod bolt thread engagement in the rod body |
| conrod / conrod_bolt_b | contact | 4 | 0.000 | 0.000 | yes | rod bolt thread engagement in the rod body |
| conrod / conrod_cap | contact | 4 | 0.000 | 0.000 | yes | rod/cap split face at the big-end axis |
| conrod / wrist_pin | contact | 4 | 0.000 | 0.000 | yes | pin in the conrod small end |
| conrod_bolt_a / conrod_cap | contact | 4 | 0.000 | 0.000 | yes | rod bolt head on the cap spot face |
| conrod_bolt_b / conrod_cap | contact | 4 | 0.000 | 0.000 | yes | rod bolt head on the cap spot face |
| crank_sprocket / crankshaft | contact | 1 | 0.000 | 0.000 | yes | timing pulley on the crank snout |
| crank_sprocket / damper_hub | contact | 1 | 0.000 | 0.000 | yes | damper hub clamps the pulley axially |
| crank_sprocket / timing_belt | contact | 1 | 0.000 | 0.000 | yes | belt on the crank pulley rim |
| crankshaft / damper_hub | contact | 1 | 0.000 | 0.000 | yes | damper hub on the crank snout |
| crankshaft / flywheel | contact | 1 | 0.000 | 0.000 | yes | flywheel on the crank flange (pilot face) |
| crankshaft / flywheel_bolt | contact | 6 | 0.000 | 0.000 | yes | flywheel bolt thread engagement in the flange |
| crankshaft / main_shell_lower | contact | 5 | 0.025 | 0.000 | yes | main running clearance 25 um radial (spec/dynamics.yaml bearings.main) = oil film |
| crankshaft / main_shell_upper | contact | 5 | 0.025 | 0.000 | yes | main running clearance 25 um radial (spec/dynamics.yaml bearings.main) = oil film |
| damper_hub / damper_rubber | contact | 1 | 0.000 | 0.000 | yes | rubber bonded to the hub rim |
| damper_ring / damper_rubber | contact | 1 | 0.000 | 0.000 | yes | inertia ring bonded to the rubber |
| exhaust_header / head | contact | 1 | 0.000 | 0.000 | yes | header flanges on the head exhaust face |
| exhaust_header / manifold_bolt_exhaust | contact | 8 | 0.000 | 0.000 | yes | bolt head on the header flange |
| flywheel / flywheel_bolt | contact | 6 | 0.000 | 0.000 | yes | flywheel bolt head on the flywheel hub |
| flywheel / ring_gear | contact | 1 | 0.000 | 0.000 | yes | ring gear shrunk on the flywheel seat |
| guide_exhaust / head | contact | 8 | 0.000 | 0.000 | yes | valve guide pressed in its bore |
| guide_exhaust / valve_exhaust | contact | 8 | 0.030 | 0.000 | yes | stem in the guide |
| guide_intake / head | contact | 8 | 0.000 | 0.000 | yes | valve guide pressed in its bore |
| guide_intake / valve_intake | contact | 8 | 0.030 | 0.000 | yes | stem in the guide |
| head / head_bolt | contact | 10 | 0.000 | 0.000 | yes | head bolt head on its seat in the head |
| head / head_gasket | contact | 1 | 0.000 | 0.000 | yes | head face on the gasket |
| head / intake_bracket | contact | 1 | 0.000 | 0.000 | yes | stay plate also bears on the flush head rear face |
| head / intake_manifold | contact | 1 | 0.000 | 0.000 | yes | intake flange on the head intake face |
| head / manifold_bolt_exhaust | contact | 8 | 0.000 | 0.000 | yes | manifold bolt thread engagement in the head |
| head / manifold_bolt_intake | contact | 8 | 0.000 | 0.000 | yes | manifold bolt thread engagement in the head |
| head / seat_insert_exhaust | contact | 8 | 0.000 | 0.000 | yes | seat insert pressed in its counterbore |
| head / seat_insert_intake | contact | 8 | 0.000 | 0.000 | yes | seat insert pressed in its counterbore |
| head / spark_plug | contact | 4 | 0.000 | 0.000 | yes | plug thread engagement and seat on the bore step |
| head / spring_seat_exhaust | contact | 8 | 0.000 | 0.000 | yes | spring seat on the pocket floor |
| head / spring_seat_intake | contact | 8 | 0.000 | 0.000 | yes | spring seat on the pocket floor |
| head_gasket / intake_bracket | contact | 1 | 0.000 | 0.000 | yes | stay plate across the flush gasket edge |
| head_gasket / liner | contact | 4 | 0.000 | 0.000 | yes | liner top flush with the deck under the gasket |
| idler / timing_belt | contact | 1 | 0.000 | 0.000 | yes | belt on the idler |
| intake_bracket / intake_manifold | contact | 1 | 0.000 | 0.000 | yes | stay tab under the plenum floor pad |
| intake_manifold / manifold_bolt_intake | contact | 8 | 0.000 | 0.000 | yes | bolt head on the intake flange |
| intake_manifold / throttle_body | contact | 1 | 0.000 | 0.000 | yes | throttle-body flange on the plenum neck flange |
| liner / piston | contact | 4 | 0.020 | 0.000 | yes | piston skirt in the liner bore |
| liner / ring | contact | 12 | 0.000 | 0.000 | yes | ring running face on the liner bore |
| main_cap / main_cap_bolt | contact | 10 | 0.000 | 0.000 | yes | cap bolt head on the cap |
| main_cap / main_shell_lower | contact | 5 | 0.000 | 0.000 | yes | lower main shell in the cap |
| main_cap / main_shell_upper | contact | 5 | 0.000 | 0.000 | yes | upper shell joint face meets the cap face |
| main_shell_lower / main_shell_upper | contact | 5 | 0.000 | 0.000 | yes | shell joint faces at the split |
| oil_pan / pan_bolt | contact | 16 | 0.000 | 0.000 | yes | pan bolt head under the pan flange |
| piston / ring | contact | 12 | 0.025 | 0.000 | yes | ring in its groove (25 um side clearance each side) |
| piston / wrist_pin | contact | 4 | 0.000 | 0.000 | yes | pin in the pin-boss bores |
| retainer_exhaust / spring_exhaust | contact | 8 | 0.000 | 0.000 | yes | spring top under the retainer (installed length) |
| retainer_intake / spring_intake | contact | 8 | 0.000 | 0.000 | yes | spring top under the retainer (installed length) |
| seat_insert_exhaust / valve_exhaust | contact | 4 | 0.000 | 0.000 | yes | valve face on its seat when closed |
| seat_insert_intake / valve_intake | contact | 4 | 0.000 | 0.000 | yes | valve face on its seat when closed |
| shim_exhaust / valve_exhaust | contact | 8 | 0.000 | 0.000 | yes | shim on the stem tip |
| shim_intake / valve_intake | contact | 8 | 0.000 | 0.000 | yes | shim on the stem tip |
| spring_exhaust / spring_seat_exhaust | contact | 8 | 0.000 | 0.000 | yes | spring bottom on its seat |
| spring_intake / spring_seat_intake | contact | 8 | 0.000 | 0.000 | yes | spring bottom on its seat |
| tensioner / timing_belt | contact | 1 | 0.000 | 0.000 | yes | belt on the tensioner |
| big_end_shell_lower / block | clear | 4 | 26.500 | 0.000 | - |  |
| big_end_shell_lower / oil_pan | clear | 2 | 35.000 | 0.000 | - |  |
| big_end_shell_upper / block | clear | 4 | 26.500 | 0.000 | - |  |
| block / clip | clear | 8 | 9.881 | 0.000 | - |  |
| block / conrod | clear | 4 | 24.500 | 0.000 | - |  |
| block / conrod_bolt_a | clear | 4 | 29.500 | 0.000 | - |  |
| block / conrod_bolt_b | clear | 4 | 29.500 | 0.000 | - |  |
| block / conrod_cap | clear | 4 | 24.500 | 0.000 | - |  |
| block / crankshaft | clear | 1 | 1.000 | 0.000 | - |  |
| block / exhaust_header | clear | 1 | 16.485 | 0.000 | - |  |
| block / intake_manifold | clear | 1 | 12.485 | 0.000 | - |  |
| block / piston | clear | 4 | 1.020 | 0.000 | - |  |
| block / ring | clear | 12 | 1.908 | 0.000 | - |  |
| block / valve_exhaust | clear | 2 | 5.739 | 0.000 | - |  |
| block / valve_intake | clear | 2 | 5.063 | 0.000 | - |  |
| block / wrist_pin | clear | 4 | 11.800 | 0.000 | - |  |
| bucket_exhaust / camshaft_exhaust | clear | 4 | 0.250 | 0.000 | - | lobe on the bucket crown while the valve lifts (lash on the base circle) |
| bucket_exhaust / collet_exhaust | clear | 16 | 1.031 | 0.000 | - |  |
| bucket_exhaust / guide_exhaust | clear | 2 | 12.352 | 0.000 | - |  |
| bucket_exhaust / intake_manifold | clear | 8 | 62.828 | 0.000 | - |  |
| bucket_exhaust / retainer_exhaust | clear | 8 | 1.000 | 0.000 | - |  |
| bucket_exhaust / spring_exhaust | clear | 8 | 3.061 | 0.000 | - |  |
| bucket_exhaust / valve_exhaust | clear | 8 | 2.250 | 0.000 | - |  |
| bucket_intake / camshaft_intake | clear | 4 | 0.200 | 0.000 | - | lobe on the bucket crown while the valve lifts (lash on the base circle) |
| bucket_intake / collet_intake | clear | 16 | 1.031 | 0.000 | - |  |
| bucket_intake / guide_intake | clear | 2 | 12.024 | 0.000 | - |  |
| bucket_intake / intake_manifold | clear | 8 | 34.204 | 0.000 | - |  |
| bucket_intake / retainer_intake | clear | 8 | 1.000 | 0.000 | - |  |
| bucket_intake / spring_intake | clear | 8 | 3.061 | 0.000 | - |  |
| bucket_intake / valve_intake | clear | 8 | 2.250 | 0.000 | - |  |
| cam_cap / intake_manifold | clear | 10 | 33.621 | 0.000 | - |  |
| cam_cap_bolt / camshaft_exhaust | clear | 10 | 1.500 | 0.000 | - |  |
| cam_cap_bolt / camshaft_intake | clear | 10 | 1.500 | 0.000 | - |  |
| cam_cap_bolt / intake_manifold | clear | 15 | 31.471 | 0.000 | - |  |
| camshaft_exhaust / collet_exhaust | clear | 12 | 8.253 | 0.000 | - |  |
| camshaft_exhaust / head_bolt | clear | 5 | 2.040 | 0.000 | - |  |
| camshaft_exhaust / intake_manifold | clear | 1 | 31.298 | 0.000 | - |  |
| camshaft_exhaust / retainer_exhaust | clear | 8 | 8.000 | 0.000 | - |  |
| camshaft_exhaust / shim_exhaust | clear | 8 | 3.250 | 0.000 | - |  |
| camshaft_exhaust / spring_exhaust | clear | 4 | 14.911 | 0.000 | - |  |
| camshaft_exhaust / timing_belt | clear | 1 | 42.790 | 0.000 | - |  |
| camshaft_exhaust / valve_exhaust | clear | 8 | 5.750 | 0.000 | - |  |
| camshaft_intake / collet_intake | clear | 12 | 8.203 | 0.000 | - |  |
| camshaft_intake / head_bolt | clear | 5 | 4.097 | 0.000 | - |  |
| camshaft_intake / intake_manifold | clear | 1 | 32.497 | 0.000 | - |  |
| camshaft_intake / retainer_intake | clear | 8 | 8.000 | 0.000 | - |  |
| camshaft_intake / shim_intake | clear | 8 | 3.200 | 0.000 | - |  |
| camshaft_intake / spring_intake | clear | 4 | 14.862 | 0.000 | - |  |
| camshaft_intake / timing_belt | clear | 1 | 42.790 | 0.000 | - |  |
| camshaft_intake / valve_intake | clear | 8 | 5.700 | 0.000 | - |  |
| clip / intake_bracket | clear | 1 | 44.800 | 0.000 | - |  |
| clip / liner | clear | 8 | 11.430 | 0.000 | - |  |
| clip / ring | clear | 8 | 9.508 | 0.000 | - |  |
| clip / wrist_pin | clear | 8 | 0.300 | 0.000 | - | circlip retains the pin end face |
| collet_exhaust / collet_exhaust | clear | 8 | 0.800 | 0.000 | - |  |
| collet_exhaust / head | clear | 16 | 11.780 | 0.000 | - |  |
| collet_exhaust / intake_manifold | clear | 16 | 76.412 | 0.000 | - |  |
| collet_exhaust / shim_exhaust | clear | 8 | 2.500 | 0.000 | - |  |
| collet_exhaust / spring_exhaust | clear | 16 | 4.775 | 0.000 | - |  |
| collet_intake / collet_intake | clear | 8 | 0.800 | 0.000 | - |  |
| collet_intake / head | clear | 16 | 11.780 | 0.000 | - |  |
| collet_intake / intake_manifold | clear | 16 | 47.874 | 0.000 | - |  |
| collet_intake / shim_intake | clear | 8 | 2.500 | 0.000 | - |  |
| collet_intake / spring_intake | clear | 16 | 4.775 | 0.000 | - |  |
| conrod / crankshaft | clear | 4 | 1.000 | 0.000 | - |  |
| conrod / intake_bracket | clear | 1 | 63.500 | 0.000 | - |  |
| conrod / liner | clear | 4 | 24.935 | 0.000 | - |  |
| conrod / piston | clear | 4 | 0.500 | 0.000 | - |  |
| conrod / ring | clear | 4 | 23.013 | 0.000 | - |  |
| conrod_bolt_a / crankshaft | clear | 4 | 3.750 | 0.000 | - |  |
| conrod_bolt_a / oil_pan | clear | 2 | 41.250 | 0.000 | - |  |
| conrod_bolt_b / crankshaft | clear | 4 | 3.750 | 0.000 | - |  |
| conrod_bolt_b / oil_pan | clear | 2 | 41.250 | 0.000 | - |  |
| conrod_cap / crankshaft | clear | 4 | 1.000 | 0.000 | - |  |
| conrod_cap / oil_pan | clear | 2 | 28.200 | 0.000 | - |  |
| crankshaft / damper_ring | clear | 1 | 39.000 | 0.000 | - |  |
| crankshaft / damper_rubber | clear | 1 | 35.000 | 0.000 | - |  |
| crankshaft / main_cap | clear | 5 | 1.000 | 0.000 | - |  |
| crankshaft / main_cap_bolt | clear | 10 | 4.343 | 0.000 | - |  |
| crankshaft / oil_pan | clear | 1 | 20.100 | 0.000 | - |  |
| crankshaft / ring_gear | clear | 1 | 84.890 | 0.000 | - |  |
| crankshaft / timing_belt | clear | 1 | 10.052 | 0.000 | - |  |
| damper_hub / damper_ring | clear | 1 | 4.000 | 0.000 | - |  |
| exhaust_header / oil_pan | clear | 1 | 45.564 | 0.000 | - |  |
| exhaust_header / pan_bolt | clear | 6 | 44.223 | 0.000 | - |  |
| flywheel_bolt / ring_gear | clear | 6 | 97.532 | 0.000 | - |  |
| guide_exhaust / intake_bracket | clear | 1 | 50.250 | 0.000 | - |  |
| guide_exhaust / intake_manifold | clear | 8 | 100.300 | 0.000 | - |  |
| guide_exhaust / spring_exhaust | clear | 8 | 3.263 | 0.000 | - |  |
| guide_exhaust / spring_seat_exhaust | clear | 8 | 0.300 | 0.000 | - |  |
| guide_intake / intake_manifold | clear | 8 | 45.898 | 0.000 | - |  |
| guide_intake / spring_intake | clear | 8 | 3.263 | 0.000 | - |  |
| guide_intake / spring_seat_intake | clear | 8 | 0.300 | 0.000 | - |  |
| head / retainer_exhaust | clear | 8 | 3.030 | 0.000 | - |  |
| head / retainer_intake | clear | 8 | 3.030 | 0.000 | - |  |
| head / shim_exhaust | clear | 8 | 11.530 | 0.000 | - |  |
| head / shim_intake | clear | 8 | 11.530 | 0.000 | - |  |
| head / spring_exhaust | clear | 8 | 0.561 | 0.000 | - |  |
| head / spring_intake | clear | 8 | 0.561 | 0.000 | - |  |
| head / valve_exhaust | clear | 8 | 1.475 | 0.000 | - |  |
| head / valve_intake | clear | 8 | 1.078 | 0.000 | - |  |
| head_bolt / head_gasket | clear | 10 | 0.500 | 0.000 | - |  |
| head_bolt / intake_manifold | clear | 10 | 24.149 | 0.000 | - |  |
| head_gasket / intake_manifold | clear | 1 | 12.485 | 0.000 | - |  |
| head_gasket / valve_exhaust | clear | 2 | 3.623 | 0.000 | - |  |
| head_gasket / valve_intake | clear | 2 | 2.956 | 0.000 | - |  |
| intake_bracket / liner | clear | 1 | 28.500 | 0.000 | - |  |
| intake_bracket / piston | clear | 1 | 31.020 | 0.000 | - |  |
| intake_bracket / ring | clear | 3 | 31.000 | 0.000 | - |  |
| intake_bracket / seat_insert_exhaust | clear | 1 | 39.750 | 0.000 | - |  |
| intake_bracket / seat_insert_intake | clear | 2 | 36.948 | 0.000 | - |  |
| intake_bracket / spark_plug | clear | 1 | 64.762 | 0.000 | - |  |
| intake_bracket / spring_seat_exhaust | clear | 1 | 44.336 | 0.000 | - |  |
| intake_bracket / valve_exhaust | clear | 1 | 41.250 | 0.000 | - |  |
| intake_bracket / valve_intake | clear | 2 | 38.351 | 0.000 | - |  |
| intake_bracket / wrist_pin | clear | 1 | 46.000 | 0.000 | - |  |
| intake_manifold / liner | clear | 4 | 44.500 | 0.000 | - |  |
| intake_manifold / retainer_exhaust | clear | 8 | 71.071 | 0.000 | - |  |
| intake_manifold / retainer_intake | clear | 8 | 38.662 | 0.000 | - |  |
| intake_manifold / seat_insert_exhaust | clear | 8 | 95.824 | 0.000 | - |  |
| intake_manifold / seat_insert_intake | clear | 8 | 53.700 | 0.000 | - |  |
| intake_manifold / shim_exhaust | clear | 8 | 72.499 | 0.000 | - |  |
| intake_manifold / shim_intake | clear | 8 | 49.697 | 0.000 | - |  |
| intake_manifold / spark_plug | clear | 4 | 83.764 | 0.000 | - |  |
| intake_manifold / spring_exhaust | clear | 8 | 78.347 | 0.000 | - |  |
| intake_manifold / spring_intake | clear | 8 | 38.485 | 0.000 | - |  |
| intake_manifold / spring_seat_exhaust | clear | 8 | 107.139 | 0.000 | - |  |
| intake_manifold / spring_seat_intake | clear | 8 | 42.458 | 0.000 | - |  |
| intake_manifold / valve_exhaust | clear | 8 | 75.862 | 0.000 | - |  |
| intake_manifold / valve_intake | clear | 8 | 47.954 | 0.000 | - |  |
| liner / valve_exhaust | clear | 2 | 3.346 | 0.000 | - |  |
| liner / valve_intake | clear | 2 | 2.632 | 0.000 | - |  |
| liner / wrist_pin | clear | 4 | 13.096 | 0.000 | - |  |
| retainer_exhaust / shim_exhaust | clear | 8 | 2.500 | 0.000 | - |  |
| retainer_exhaust / valve_exhaust | clear | 8 | 1.141 | 0.000 | - |  |
| retainer_intake / shim_intake | clear | 8 | 2.500 | 0.000 | - |  |
| retainer_intake / valve_intake | clear | 8 | 1.141 | 0.000 | - |  |
| seat_insert_exhaust / spark_plug | clear | 8 | 4.883 | 0.000 | - |  |
| seat_insert_exhaust / valve_exhaust | clear | 4 | 0.770 | 0.000 | - | valve face on its seat when closed |
| seat_insert_intake / spark_plug | clear | 8 | 3.190 | 0.000 | - |  |
| seat_insert_intake / valve_intake | clear | 4 | 0.836 | 0.000 | - | valve face on its seat when closed |
| spark_plug / valve_exhaust | clear | 8 | 6.383 | 0.000 | - |  |
| spark_plug / valve_intake | clear | 8 | 4.690 | 0.000 | - |  |
| spring_exhaust / valve_exhaust | clear | 8 | 6.013 | 0.000 | - |  |
| spring_intake / valve_intake | clear | 8 | 6.013 | 0.000 | - |  |
| spring_seat_exhaust / valve_exhaust | clear | 8 | 3.050 | 0.000 | - |  |
| spring_seat_intake / valve_intake | clear | 8 | 3.050 | 0.000 | - |  |

## Key clearances (selected pairs, exact minimum distance in the static pose)

Pairs whose boxes are far apart never reach the exact check above; these are measured anyway because they matter for packaging (spec/cad_fits.yaml report_clearances). Rotating partners (flywheel, ring gear) are axisymmetric, so the static distance holds over the cycle.

| Pair (kinds) | Min distance mm | Closest instances | Why it matters |
|---|---|---|---|
| exhaust_header / oil_pan | 45.6 | exhaust_header / oil_pan | hot 4-2-1 header and collector past the sump |
| exhaust_header / flywheel | 77.1 | exhaust_header / flywheel | collector and outlet near the flywheel |
| exhaust_header / ring_gear | 73.1 | exhaust_header / ring_gear | collector and outlet near the starter ring gear |
| exhaust_header / block | 16.5 | exhaust_header / block | primaries along the block exhaust side |
| exhaust_header / intake_* | 72.0 | exhaust_header / intake_bracket | hot header to the intake system and stay |
| intake_manifold / camshaft_* | 31.3 | intake_manifold / camshaft_exhaust | runners and plenum over the cams |
| intake_manifold / cam_cap* | 31.5 | intake_manifold / cam_cap_bolt_intake_2b | plenum floor over the cam caps and their bolts |
| intake_manifold / cam_pulley_* | 36.7 | intake_manifold / cam_pulley_exhaust | plenum front end to the cam pulleys |
| intake_manifold / timing_belt | 35.9 | intake_manifold / timing_belt | plenum front end to the timing belt |
| throttle_body / flywheel | 267.6 | throttle_body / flywheel | throttle body at the rear end |

## Moving pairs (minimum over the cycle or the lift range, aggregated by part family)

| Family A | Family B | Param | VVT deg | Mode | Pairs | Min mm | Worst pair at (deg or lift mm) | Certified >= mm | Result |
|---|---|---|---|---|---|---|---|---|---|
| crank_sprocket | timing_belt | crank | 0 | contact | 1 | 0.000 | crank_sprocket / timing_belt at 18.00 | - | running contact, no penetration |
| liner | ring | crank | 0 | contact | 12 | 0.000 | liner_1 / ring_1_1 at 0.00 | - | running contact, no penetration |
| conrod | wrist_pin | crank | 0 | contact | 4 | 0.000 | conrod_1 / wrist_pin_1 at 0.00 | - | running contact, no penetration |
| valve | seat_insert | lift | 0 | contact | 16 | 0.000 | valve_1_intake_1 / seat_insert_1_intake_1 at 0.00 | - | running contact, no penetration |
| bucket | camshaft_intake | crank | 0 | contact | 8 | 0.000 | bucket_2_intake_1 / camshaft_intake at 434.00 | - | running contact, no penetration |
| cam_pulley_intake | timing_belt | crank | 0 | contact | 1 | 0.000 | cam_pulley_intake / timing_belt at 0.00 | - | running contact, no penetration |
| cam_pulley_intake | timing_belt | crank | 40 | contact | 1 | 0.000 | cam_pulley_intake / timing_belt at 0.00 | - | running contact, no penetration |
| cam_pulley_exhaust | timing_belt | crank | 0 | contact | 1 | 0.000 | cam_pulley_exhaust / timing_belt at 0.00 | - | running contact, no penetration |
| bucket | camshaft_exhaust | crank | 0 | contact | 8 | 0.000 | bucket_3_exhaust_1 / camshaft_exhaust at 452.00 | - | running contact, no penetration |
| liner | piston | crank | 0 | contact | 4 | 0.020 | liner_1 / piston_1 at 2.00 | - | running contact, no penetration |
| crankshaft | big_end_shell_upper | crank | 0 | contact | 4 | 0.022 | crankshaft / big_end_shell_upper_2 at 352.00 | - | running contact, no penetration |
| crankshaft | big_end_shell_lower | crank | 0 | contact | 4 | 0.022 | crankshaft / big_end_shell_lower_2 at 352.00 | - | running contact, no penetration |
| crankshaft | main_shell_upper | crank | 0 | contact | 5 | 0.025 | crankshaft / main_shell_upper_1 at 188.00 | - | running contact, no penetration |
| crankshaft | main_shell_lower | crank | 0 | contact | 5 | 0.025 | crankshaft / main_shell_lower_5 at 144.00 | - | running contact, no penetration |
| head | bucket | lift | 0 | contact | 16 | 0.030 | head / bucket_1_exhaust_1 at 6.50 | - | running contact, no penetration |
| valve | guide | lift | 0 | contact | 16 | 0.030 | valve_1_exhaust_1 / guide_1_exhaust_1 at 3.00 | - | running contact, no penetration |
| head | camshaft_exhaust | crank | 0 | contact | 1 | 0.040 | head / camshaft_exhaust at 0.00 | - | running contact, no penetration |
| head | camshaft_intake | crank | 0 | contact | 1 | 0.040 | head / camshaft_intake at 0.00 | - | running contact, no penetration |
| head | camshaft_intake | crank | 40 | contact | 1 | 0.040 | head / camshaft_intake at 0.00 | - | running contact, no penetration |
| camshaft_intake | cam_cap_intake | crank | 0 | contact | 5 | 0.040 | camshaft_intake / cam_cap_intake_3 at 236.00 | - | running contact, no penetration |
| camshaft_intake | cam_cap_intake | crank | 40 | contact | 5 | 0.040 | camshaft_intake / cam_cap_intake_3 at 196.00 | - | running contact, no penetration |
| camshaft_exhaust | cam_cap_exhaust | crank | 0 | contact | 5 | 0.040 | camshaft_exhaust / cam_cap_exhaust_1 at 362.00 | - | running contact, no penetration |
| conrod | piston | crank | 0 | clearance | 4 | 0.500 | conrod_1 / piston_1 at 2.00 | 0.000 | clear |
| head_gasket | piston | crank | 0 | clearance | 4 | 0.721 | head_gasket / piston_2 at 180.00 | 0.388 | clear |
| crankshaft | conrod | crank | 0 | clearance | 4 | 1.000 | crankshaft / conrod_1 at 124.00 | 0.000 | clear |
| piston | valve | crank | 40 | clearance | 8 | 1.000 | piston_4 / valve_4_intake_2 at 0.00 | 0.524 | clear |
| block | crankshaft | crank | 0 | clearance | 1 | 1.000 | block / crankshaft at 4.00 | 0.216 | clear |
| crankshaft | main_cap | crank | 0 | clearance | 5 | 1.000 | crankshaft / main_cap_2 at 0.00 | 0.216 | clear |
| crankshaft | conrod_cap | crank | 0 | clearance | 4 | 1.000 | crankshaft / conrod_cap_2 at 10.00 | 0.000 | clear |
| block | piston | crank | 0 | clearance | 4 | 1.020 | block / piston_1 at 171.38 | 0.568 | clear |
| head | valve | lift | 0 | clearance | 16 | 1.078 | head / valve_4_intake_2 at 0.00 | 0.815 | clear |
| head | piston | crank | 0 | clearance | 4 | 1.200 | head / piston_1 at 0.00 | 0.793 | clear |
| camshaft_intake | cam_cap_bolt_intake | crank | 0 | clearance | 10 | 1.500 | camshaft_intake / cam_cap_bolt_intake_1b at 76.00 | 1.124 | clear |
| camshaft_intake | cam_cap_bolt_intake | crank | 40 | clearance | 10 | 1.500 | camshaft_intake / cam_cap_bolt_intake_1b at 36.00 | 1.124 | clear |
| camshaft_exhaust | cam_cap_bolt_exhaust | crank | 0 | clearance | 10 | 1.500 | camshaft_exhaust / cam_cap_bolt_exhaust_1a at 586.00 | 1.133 | clear |
| block | ring | crank | 0 | clearance | 12 | 1.908 | block / ring_3_4 at 180.00 | 1.094 | clear |
| block | flywheel_bolt | crank | 0 | clearance | 6 | 2.000 | block / flywheel_bolt_1 at 0.00 | 1.230 | clear |
| head_bolt | camshaft_exhaust | crank | 0 | clearance | 5 | 2.040 | head_bolt_2 / camshaft_exhaust at 0.00 | 1.673 | clear |
| valve | seat_insert | lift | 0 | clearance | 20 | 2.500 | valve_1_intake_1 / seat_insert_1_intake_2 at 0.00 | 1.836 | clear |
| seat_insert | valve | lift | 0 | clearance | 20 | 2.500 | seat_insert_1_intake_1 / valve_1_intake_2 at 0.00 | 1.836 | clear |
| liner | valve | lift | 0 | clearance | 16 | 2.629 | liner_1 / valve_1_intake_2 at 6.94 | 2.367 | clear |
| piston | valve | crank | 0 | clearance | 16 | 2.738 | piston_1 / valve_1_intake_1 at 366.64 | 1.838 | clear |
| head_gasket | valve | lift | 0 | clearance | 16 | 2.956 | head_gasket / valve_1_intake_1 at 6.00 | 2.694 | clear |
| shim | camshaft_intake | crank | 0 | clearance | 8 | 3.000 | shim_4_intake_1 / camshaft_intake at 253.17 | 2.491 | clear |
| shim | camshaft_exhaust | crank | 0 | clearance | 8 | 3.000 | shim_4_exhaust_1 / camshaft_exhaust at 466.98 | 2.515 | clear |
| head | retainer | lift | 0 | clearance | 16 | 3.030 | head / retainer_2_exhaust_1 at 2.00 | 2.767 | clear |
| valve | spring_seat | lift | 0 | clearance | 16 | 3.050 | valve_2_intake_2 / spring_seat_2_intake_2 at 10.25 | 2.787 | clear |
| bucket | bucket | crank | 0 | clearance | 8 | 3.500 | bucket_1_exhaust_1 / bucket_1_exhaust_2 at 0.00 | 3.262 | clear |
| crankshaft | conrod_bolt_b | crank | 0 | clearance | 4 | 3.750 | crankshaft / conrod_bolt_b_2 at 352.00 | 2.484 | clear |
| crankshaft | conrod_bolt_a | crank | 0 | clearance | 4 | 3.750 | crankshaft / conrod_bolt_a_1 at 252.00 | 2.484 | clear |
| damper_hub | timing_belt | crank | 0 | clearance | 1 | 3.856 | damper_hub / timing_belt at 226.00 | 2.534 | clear |
| crankshaft | main_cap_bolt | crank | 0 | clearance | 10 | 4.000 | crankshaft / main_cap_bolt_3b at 56.00 | 2.432 | clear |
| block | crank_sprocket | crank | 0 | clearance | 1 | 4.000 | block / crank_sprocket at 0.00 | 0.235 | clear |
| head | cam_pulley_intake | crank | 0 | clearance | 1 | 4.000 | head / cam_pulley_intake at 0.00 | 0.385 | clear |
| head | cam_pulley_intake | crank | 40 | clearance | 1 | 4.000 | head / cam_pulley_intake at 0.00 | 0.385 | clear |
| head | cam_pulley_exhaust | crank | 0 | clearance | 1 | 4.000 | head / cam_pulley_exhaust at 0.00 | 0.385 | clear |
| valve | valve | crank | 0 | clearance | 20 | 4.000 | valve_1_intake_1 / valve_1_intake_2 at 0.00 | 3.735 | clear |
| ring | valve | crank | 40 | clearance | 24 | 4.035 | ring_1_1 / valve_1_intake_2 at 364.77 | 3.131 | clear |
| head_bolt | camshaft_intake | crank | 0 | clearance | 5 | 4.097 | head_bolt_1 / camshaft_intake at 0.00 | 3.721 | clear |
| head_bolt | camshaft_intake | crank | 40 | clearance | 5 | 4.097 | head_bolt_1 / camshaft_intake at 0.00 | 3.721 | clear |
| head_bolt | bucket | lift | 0 | clearance | 16 | 4.366 | head_bolt_3 / bucket_1_intake_2 at 0.00 | 2.836 | clear |
| piston | seat_insert | crank | 0 | clearance | 16 | 4.477 | piston_1 / seat_insert_1_exhaust_1 at 0.00 | 0.366 | clear |
| spark_plug | valve | lift | 0 | clearance | 16 | 4.690 | spark_plug_1 / valve_1_intake_2 at 0.01 | 4.427 | clear |
| valve | valve | crank | 40 | clearance | 16 | 4.870 | valve_4_intake_2 / valve_4_exhaust_2 at 702.52 | 4.654 | clear |
| block | valve | lift | 0 | clearance | 16 | 5.063 | block / valve_4_intake_2 at 6.44 | 4.802 | clear |
| valve | camshaft_intake | crank | 0 | clearance | 8 | 5.500 | valve_4_intake_2 / camshaft_intake at 253.17 | 4.991 | clear |
| valve | camshaft_exhaust | crank | 0 | clearance | 8 | 5.500 | valve_1_exhaust_2 / camshaft_exhaust at 106.98 | 5.015 | clear |
| damper_ring | timing_belt | crank | 0 | clearance | 1 | 5.500 | damper_ring / timing_belt at 112.00 | 3.603 | clear |
| damper_rubber | timing_belt | crank | 0 | clearance | 1 | 5.500 | damper_rubber / timing_belt at 0.00 | 4.075 | clear |
| head_gasket | ring | crank | 0 | clearance | 4 | 5.548 | head_gasket / ring_1_1 at 0.00 | 1.414 | clear |
| head | ring | crank | 0 | clearance | 4 | 6.225 | head / ring_1_1 at 0.00 | 2.114 | clear |
| liner | ring | crank | 0 | clearance | 9 | 7.500 | liner_1 / ring_1_2 at 0.00 | 3.389 | clear |
| ring | liner | crank | 0 | clearance | 9 | 7.500 | ring_1_1 / liner_2 at 0.00 | 3.389 | clear |
| liner | piston | crank | 0 | clearance | 3 | 7.520 | liner_1 / piston_2 at 0.00 | 3.409 | clear |
| piston | liner | crank | 0 | clearance | 3 | 7.520 | piston_1 / liner_2 at 0.00 | 3.409 | clear |
| ring | valve | crank | 0 | clearance | 40 | 7.961 | ring_1_4 / valve_4_intake_2 at 6.45 | 7.066 | clear |
| retainer | camshaft_intake | crank | 0 | clearance | 8 | 8.000 | retainer_1_intake_1 / camshaft_intake at 613.17 | 7.491 | clear |
| collet_b | camshaft_intake | crank | 0 | clearance | 8 | 8.000 | collet_b_1_intake_2 / camshaft_intake at 613.17 | 7.491 | clear |
| collet | camshaft_intake | crank | 0 | clearance | 8 | 8.000 | collet_2_intake_1 / camshaft_intake at 142.83 | 7.491 | clear |
| retainer | camshaft_exhaust | crank | 0 | clearance | 8 | 8.000 | retainer_4_exhaust_1 / camshaft_exhaust at 466.98 | 7.515 | clear |
| collet_b | camshaft_exhaust | crank | 0 | clearance | 8 | 8.000 | collet_b_4_exhaust_1 / camshaft_exhaust at 37.02 | 7.515 | clear |
| collet | camshaft_exhaust | crank | 0 | clearance | 8 | 8.000 | collet_3_exhaust_1 / camshaft_exhaust at 286.98 | 7.515 | clear |
| ring_gear | intake_bracket | crank | 0 | clearance | 1 | 8.944 | ring_gear / intake_bracket at 69.98 | 5.161 | clear |
| block | clip | crank | 0 | clearance | 8 | 9.881 | block / clip_2_2 at 0.00 | 4.240 | clear |
| ring | ring | crank | 0 | clearance | 27 | 10.000 | ring_1_1 / ring_1_2 at 90.00 | 8.356 | clear |
| ring | piston | crank | 0 | clearance | 9 | 10.020 | ring_1_1 / piston_2 at 98.00 | 8.376 | clear |
| piston | ring | crank | 0 | clearance | 9 | 10.020 | piston_1 / ring_1_2 at 82.00 | 8.376 | clear |
| piston | piston | crank | 0 | clearance | 3 | 10.040 | piston_1 / piston_2 at 90.00 | 8.396 | clear |
| crankshaft | timing_belt | crank | 0 | clearance | 1 | 10.052 | crankshaft / timing_belt at 132.00 | 8.484 | clear |
| bucket | guide | lift | 0 | clearance | 16 | 10.249 | bucket_1_intake_2 / guide_1_intake_2 at 10.50 | 10.233 | clear |
| block | conrod_bolt_b | crank | 0 | clearance | 4 | 10.646 | block / conrod_bolt_b_1 at 144.20 | 10.006 | clear |
| block | conrod_bolt_a | crank | 0 | clearance | 4 | 10.646 | block / conrod_bolt_a_1 at 215.80 | 10.006 | clear |
| liner | clip | crank | 0 | clearance | 8 | 11.430 | liner_1 / clip_2_1 at 0.00 | 5.998 | clear |
| head | collet_b | lift | 0 | clearance | 16 | 11.780 | head / collet_b_4_exhaust_2 at 0.00 | 3.862 | clear |
| head | collet | lift | 0 | clearance | 16 | 11.780 | head / collet_4_exhaust_2 at 0.00 | 3.862 | clear |
| head | shim | lift | 0 | clearance | 16 | 11.798 | head / shim_4_exhaust_2 at 0.00 | 4.858 | clear |
| block | wrist_pin | crank | 0 | clearance | 4 | 11.800 | block / wrist_pin_2 at 0.00 | 6.115 | clear |
| block | conrod_cap | crank | 0 | clearance | 4 | 12.361 | block / conrod_cap_1 at 90.00 | 11.357 | clear |
| liner | wrist_pin | crank | 0 | clearance | 4 | 13.096 | liner_1 / wrist_pin_1 at 0.00 | 7.722 | clear |
| block | ring_gear | crank | 0 | clearance | 1 | 14.000 | block / ring_gear at 0.00 | 10.217 | clear |
| block | flywheel | crank | 0 | clearance | 1 | 14.000 | block / flywheel at 0.00 | 10.372 | clear |
| oil_pan | flywheel | crank | 0 | clearance | 1 | 14.000 | oil_pan / flywheel at 70.00 | 10.372 | clear |
| liner | conrod | crank | 0 | clearance | 4 | 14.838 | liner_1 / conrod_1 at 305.88 | 13.933 | clear |
| flywheel | intake_bracket | crank | 0 | clearance | 1 | 16.125 | flywheel / intake_bracket at 2.00 | 12.496 | clear |
| valve | bucket | crank | 0 | clearance | 8 | 17.250 | valve_1_exhaust_1 / bucket_1_exhaust_2 at 0.00 | 17.012 | clear |
| bucket | valve | crank | 0 | clearance | 8 | 17.250 | bucket_1_exhaust_1 / valve_1_exhaust_2 at 0.00 | 17.012 | clear |
| block | conrod | crank | 0 | clearance | 4 | 17.251 | block / conrod_2 at 125.68 | 16.345 | clear |
| cam_pulley_intake | cam_pulley_exhaust | crank | 0 | clearance | 1 | 17.802 | cam_pulley_intake / cam_pulley_exhaust at 319.75 | 16.356 | clear |
| cam_pulley_intake | cam_pulley_exhaust | crank | 40 | clearance | 1 | 17.802 | cam_pulley_intake / cam_pulley_exhaust at 365.99 | 16.356 | clear |
| oil_pan | ring_gear | crank | 0 | clearance | 1 | 18.000 | oil_pan / ring_gear at 82.00 | 14.217 | clear |
| ring_gear | bracket_bolt_block_a | crank | 0 | clearance | 1 | 18.045 | ring_gear / bracket_bolt_block_a at 26.02 | 14.266 | clear |
| head_bolt | valve | lift | 0 | clearance | 8 | 18.586 | head_bolt_1 / valve_1_intake_1 at 0.00 | 2.336 | clear |
| oil_pan | crankshaft | crank | 0 | clearance | 1 | 20.100 | oil_pan / crankshaft at 254.00 | 18.532 | clear |
| conrod | valve | crank | 40 | clearance | 8 | 20.237 | conrod_1 / valve_1_intake_1 at 365.99 | 19.276 | clear |
| block | big_end_shell_lower | crank | 0 | clearance | 4 | 21.914 | block / big_end_shell_lower_1 at 228.25 | 20.960 | clear |
| retainer | guide | lift | 0 | clearance | 16 | 22.000 | retainer_1_exhaust_1 / guide_1_exhaust_1 at 0.00 | 3.266 | clear |
| conrod | ring | crank | 0 | clearance | 12 | 23.013 | conrod_4 / ring_3_4 at 128.00 | 21.277 | clear |
| bucket | spring_seat | lift | 0 | clearance | 16 | 24.659 | bucket_1_exhaust_2 / spring_seat_1_exhaust_2 at 0.00 | 2.110 | clear |
| conrod | valve | crank | 0 | clearance | 16 | 24.868 | conrod_2 / valve_2_intake_1 at 187.88 | 23.890 | clear |
| flywheel | bracket_bolt_block_a | crank | 0 | clearance | 1 | 27.955 | flywheel / bracket_bolt_block_a at 20.00 | 24.327 | clear |
| camshaft_exhaust | intake_bracket | crank | 0 | clearance | 1 | 28.102 | camshaft_exhaust / intake_bracket at 120.71 | 27.742 | clear |
| ring | head_bolt | crank | 0 | clearance | 24 | 29.504 | ring_1_1 / head_bolt_1 at 0.00 | 8.769 | clear |
| piston | head_bolt | crank | 0 | clearance | 8 | 29.524 | piston_4 / head_bolt_9 at 0.00 | 8.790 | clear |
| crankshaft | piston | crank | 0 | clearance | 4 | 30.000 | crankshaft / piston_2 at 0.00 | 27.619 | clear |
| camshaft_exhaust | intake_manifold | crank | 0 | clearance | 1 | 31.298 | camshaft_exhaust / intake_manifold at 0.00 | 5.978 | clear |
| ring | intake_bracket | crank | 0 | clearance | 3 | 31.673 | ring_1_4 / intake_bracket at 0.00 | 26.889 | clear |
| piston | intake_bracket | crank | 0 | clearance | 1 | 31.693 | piston_4 / intake_bracket at 0.00 | 26.909 | clear |
| camshaft_intake | intake_manifold | crank | 0 | clearance | 1 | 32.497 | camshaft_intake / intake_manifold at 0.00 | 8.444 | clear |
| camshaft_intake | intake_manifold | crank | 40 | clearance | 1 | 32.875 | camshaft_intake / intake_manifold at 0.00 | 8.444 | clear |
| oil_pan | conrod_cap | crank | 0 | clearance | 4 | 33.770 | oil_pan / conrod_cap_1 at 221.61 | 32.777 | clear |
| crankshaft | clip | crank | 0 | clearance | 8 | 34.000 | crankshaft / clip_2_2 at 0.00 | 31.619 | clear |
| bucket | intake_manifold | lift | 0 | clearance | 16 | 34.652 | bucket_1_intake_1 / intake_manifold at 0.00 | 21.146 | clear |
| camshaft_exhaust | bracket_bolt_pad_b | crank | 0 | clearance | 1 | 35.376 | camshaft_exhaust / bracket_bolt_pad_b at 0.00 | 0.174 | clear |
| crankshaft | wrist_pin | crank | 0 | clearance | 4 | 35.500 | crankshaft / wrist_pin_1 at 180.00 | 33.119 | clear |
| valve | intake_bracket | lift | 0 | clearance | 4 | 38.446 | valve_4_intake_2 / intake_bracket at 0.00 | 37.336 | clear |
| retainer | intake_manifold | lift | 0 | clearance | 16 | 40.925 | retainer_1_intake_1 / intake_manifold at 0.00 | 30.207 | clear |
| camshaft_intake | timing_belt | crank | 0 | clearance | 1 | 42.790 | camshaft_intake / timing_belt at 616.00 | 42.414 | clear |
| camshaft_intake | timing_belt | crank | 40 | clearance | 1 | 42.790 | camshaft_intake / timing_belt at 576.00 | 42.414 | clear |
| camshaft_exhaust | timing_belt | crank | 0 | clearance | 1 | 42.790 | camshaft_exhaust / timing_belt at 78.00 | 42.423 | clear |
| ring_gear | bracket_bolt_block_b | crank | 0 | clearance | 1 | 42.861 | ring_gear / bracket_bolt_block_b at 38.01 | 39.080 | clear |
| bucket | intake_bracket | lift | 0 | clearance | 2 | 43.177 | bucket_4_exhaust_2 / intake_bracket at 0.00 | 38.679 | clear |
| clip | intake_bracket | crank | 0 | clearance | 1 | 44.800 | clip_2_4 / intake_bracket at 0.00 | 40.689 | clear |
| wrist_pin | intake_bracket | crank | 0 | clearance | 1 | 46.000 | wrist_pin_4 / intake_bracket at 0.00 | 41.889 | clear |
| piston | intake_manifold | crank | 0 | clearance | 4 | 47.021 | piston_1 / intake_manifold at 0.00 | 42.909 | clear |
| ring | intake_manifold | crank | 0 | clearance | 4 | 47.293 | ring_1_1 / intake_manifold at 0.00 | 42.904 | clear |
| valve | intake_manifold | lift | 0 | clearance | 16 | 47.954 | valve_1_intake_1 / intake_manifold at 0.00 | 21.467 | clear |
| collet_b | intake_manifold | lift | 0 | clearance | 16 | 49.653 | collet_b_1_intake_1 / intake_manifold at 0.00 | 40.034 | clear |
| flywheel | bracket_bolt_block_b | crank | 0 | clearance | 1 | 51.011 | flywheel / bracket_bolt_block_b at 36.00 | 47.383 | clear |
| shim | intake_manifold | lift | 0 | clearance | 16 | 52.607 | shim_1_intake_1 / intake_manifold at 0.00 | 43.861 | clear |
| collet | intake_manifold | lift | 0 | clearance | 16 | 53.718 | collet_1_intake_1 / intake_manifold at 0.00 | 45.094 | clear |
| oil_pan | conrod_bolt_a | crank | 0 | clearance | 4 | 56.250 | oil_pan / conrod_bolt_a_2 at 0.00 | 4.776 | clear |
| oil_pan | conrod_bolt_b | crank | 0 | clearance | 4 | 56.250 | oil_pan / conrod_bolt_b_2 at 0.00 | 4.776 | clear |
| conrod | intake_bracket | crank | 0 | clearance | 1 | 63.500 | conrod_4 / intake_bracket at 0.00 | 58.932 | clear |
| camshaft_intake | cam_pulley_exhaust | crank | 0 | clearance | 1 | 63.592 | camshaft_intake / cam_pulley_exhaust at 504.00 | 62.493 | clear |
| camshaft_intake | cam_pulley_exhaust | crank | 40 | clearance | 1 | 63.592 | camshaft_intake / cam_pulley_exhaust at 464.00 | 62.493 | clear |
| block | big_end_shell_upper | crank | 0 | clearance | 4 | 67.798 | block / big_end_shell_upper_1 at 0.00 | 16.191 | clear |
| oil_pan | big_end_shell_lower | crank | 0 | clearance | 4 | 69.557 | oil_pan / big_end_shell_lower_2 at 0.00 | 14.191 | clear |
| head_gasket | ring_gear | crank | 0 | clearance | 1 | 71.386 | head_gasket / ring_gear at 69.98 | 67.604 | clear |
| head | ring_gear | crank | 0 | clearance | 1 | 72.073 | head / ring_gear at 69.98 | 68.290 | clear |
| head_gasket | flywheel | crank | 0 | clearance | 1 | 80.324 | head_gasket / flywheel at 4.00 | 76.696 | clear |
| head | flywheel | crank | 0 | clearance | 1 | 80.987 | head / flywheel at 0.00 | 0.568 | clear |
| crankshaft | exhaust_header | crank | 0 | clearance | 1 | 84.304 | crankshaft / exhaust_header at 0.00 | 38.592 | clear |
| conrod_cap | exhaust_header | crank | 0 | clearance | 4 | 104.832 | conrod_cap_1 / exhaust_header at 0.00 | 39.922 | clear |

## Valve springs and piston-to-valve

* intake spring at full lift: 25.50 mm vs solid 22.00 mm, coil-bind margin 3.50 mm (limit 1.50 mm); nearest parts at full lift: retainer_1_intake_1 0.00 mm, bucket_1_intake_1 3.06 mm, guide_1_intake_1 3.26 mm, collet_1_intake_1 4.78 mm.
* exhaust spring at full lift: 26.50 mm vs solid 22.00 mm, coil-bind margin 4.50 mm (limit 1.50 mm); nearest parts at full lift: retainer_1_exhaust_1 0.00 mm, bucket_1_exhaust_1 3.06 mm, guide_1_exhaust_1 3.26 mm, collet_1_exhaust_1 4.78 mm.

The CAD minimum is the 3D distance over the 720-deg sweep; M3 is the gap along the valve axis to the pocket floor. The pocket wall is valve_relief_margin/2 from the head all round, so the CAD minimum must lie in [min(M3, wall), M3]: below M3 when the head passes the pocket rim or wall before the floor.

| Valve / VVT advance | CAD min mm | M3 model min mm | Pocket wall mm | Expected range mm |
|---|---|---|---|---|
| intake_adv0 | 2.738 | 5.864 | 1.000 | 1.000 to 5.864 |
| intake_adv40 | 1.000 | 1.233 | 1.000 | 1.000 to 1.233 |
| exhaust_adv0 | 2.852 | 3.538 | 1.000 | 1.000 to 3.538 |

## Tessellation

Export settings: absolute chordal deflection 0.02 mm, angular 0.08 rad. Worst measured chordal deviation 0.0193 mm (limit 0.02), worst facet angle 0.0800 rad (limit 0.1), 3189972 triangles over all source solids. glTF exports carry exact per-vertex surface normals (smooth shading).

Baseline (build123d exporters, relative deflection 0.001, 0.1 rad): worst chord 0.0253 mm, facet angle 0.0914 rad, no exported normals.

## Simplifications and known gaps (honest list)

* Threads are plain cylinders at the minor/major diameter with a line-to-line fit (listed as contacts), no helix.
* No head water jacket, oil drain-backs, cam cover or timing cover; the block jacket is an open-deck annulus.
* Valve springs are static solids at the installed length in the assembly; the audit rebuilds each spring at the lift of its valve for the static pose and checks coil bind and neighbour distances at full lift instead of sweeping the deforming coil.
* Head bolts sit under the cam axes in closed pockets below the journals (bolts cannot be removed with the cams fitted, as on some production heads); the pocket roof is head_bolt_pocket_roof below the journal bore.
* Intake: swan-neck runners curl over the head into a plenum above the cam area; the 550 mm gasdyn runner length (port + flange + runner to the bellmouth lip) is met exactly, and there is no cam cover modelled under the plenum. The throttle body is held by flange contact (no flange bolts modelled) and its blade is fixed at the rendered angle.
* spec/dynamics.yaml bearing diameters (56/48 mm) differ from the CAD journals (52/45 mm, spec/cad.yaml); the CAD shells use the dynamics radial clearance and width but the CAD diameters.
* 12 clearance pairs keep a constant gap while they move: crank vs conrod and conrod cap (1.0 mm side gaps) and conrod vs piston (0.5 mm small-end side gap). The faces involved are planes normal to the common rotation axis, so the gap is exactly constant and the grid minimum is the true minimum; the Lipschitz certificate cannot see that (it bounds the surface speed, not the gap's derivative) and reports a lower bound of 0 for these pairs (Certified column).
* Running-contact pairs (bearing shells, sliding fits, cam lobe on bucket) are checked for penetration by exact common volume on a subsampled grid (spec/cad_fits.yaml sweep.contact_volume_step_deg/_mm), not certified between samples; clearance pairs carry the Lipschitz certificate.
* Belts, the oil pan gasket and sealants are not modelled as compliant parts; the belt is a rigid band at its static path.
