/* Redline ECU: hardware abstraction boundary.
 *
 * The portable core (ecu/src) calls ONLY the functions declared here. A target provides them:
 *   ecu/host/hal_host.c   host build (ctypes bridge for the Python plant, C unit tests)
 *   (M9) an RV64 HAL for Vortex-V co-sim, same core sources, no changes.
 * The platform in turn calls the core entry points in ecu.h from its capture interrupts
 * (crank / cam edge), its knock-window-complete interrupt and its slow-loop timer.
 */
#ifndef ECU_HAL_H
#define ECU_HAL_H

#include "ecu/ecu_types.h"

/* Arm an output compare: the event must take effect at absolute timer tick ev->tick
 * (wrap-safe 32-bit). Called from the angle scheduler, i.e. from the crank-edge context. */
void hal_emit_event(const ecu_event_t *ev);

#endif
