/* The sync path of the portable core as ONE translation unit for the RV64 target (M9.1, ADR 0005):
 * crank decoder, angle scheduler, knock control and the capture entry points, sources unchanged.
 * Compiled together, GCC inlines ecu_angle_wrap (12 calls per tooth), the decoder and the scheduler
 * into ecu_crank_edge. On Vortex-V every return is a JALR that stops dispatch until it resolves,
 * so the calls cost more than the work. Integer only: this object is the one the `nofloat` check
 * inspects. */
#include "../../src/decoder.c"
#include "../../src/sched.c"
#include "../../src/knock.c"
