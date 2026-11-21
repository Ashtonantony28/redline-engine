/* Per-cylinder knock control and the edge entry points. Fixed point only.
 * Knocking window (counts > threshold): retard that cylinder by one step (capped).
 * Clean window: recover by a small step toward zero retard. */
#include "ecu/ecu.h"

void ecu_knock_window(ecu_t *e, uint8_t cyl, uint32_t counts)
{
    if (cyl >= ECU_NCYL) {
        return;
    }
    int32_t r = e->knock_retard_u[cyl];
    if (counts > ECU_KNOCK_THRESHOLD) {
        r += ECU_KNOCK_STEP_U;
        e->n_knock++;
    } else {
        r -= ECU_KNOCK_RECOVER_U;
    }
    if (r > ECU_KNOCK_MAX_U) {
        r = ECU_KNOCK_MAX_U;
    }
    if (r < 0) {
        r = 0;
    }
    e->knock_retard_u[cyl] = r;
}

void ecu_crank_edge(ecu_t *e, uint32_t tick)
{
    ecu_tooth_t t;
    ecu_dec_crank(&e->dec, tick, &t);
    ecu_sched_tooth(&e->sch, &e->dec, &t);
}

void ecu_cam_edge(ecu_t *e, uint32_t tick)
{
    ecu_dec_cam(&e->dec, tick);
}
