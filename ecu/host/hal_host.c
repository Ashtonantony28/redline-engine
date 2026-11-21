/* Host HAL + ctypes bridge (libecu_host.so). One static ECU instance; output compare events are
 * queued in a fixed ring buffer that the Python plant drains after each edge call. No heap. */
#include "ecu/ecu.h"
#include "ecu/ecu_hal.h"

#define HIL_RING 256u

static ecu_t g_ecu;
static ecu_event_t g_ring[HIL_RING];
static uint32_t g_head, g_tail, g_overflow;

void hal_emit_event(const ecu_event_t *ev)
{
    uint32_t next = (g_head + 1u) % HIL_RING;
    if (next == g_tail) {
        g_overflow++;
        return;
    }
    g_ring[g_head] = *ev;
    g_head = next;
}

void hil_reset(void)
{
    ecu_init(&g_ecu);
    g_head = g_tail = g_overflow = 0u;
}

uint32_t hil_sizeof_inputs(void) { return (uint32_t)sizeof(ecu_inputs_t); }
uint32_t hil_sizeof_outputs(void) { return (uint32_t)sizeof(ecu_outputs_t); }
void hil_crank(uint32_t tick) { ecu_crank_edge(&g_ecu, tick); }
void hil_cam(uint32_t tick) { ecu_cam_edge(&g_ecu, tick); }
void hil_knock(uint32_t cyl, uint32_t counts) { ecu_knock_window(&g_ecu, (uint8_t)cyl, counts); }
void hil_slow(const ecu_inputs_t *in, ecu_outputs_t *out) { ecu_slow_step(&g_ecu, in, out); }

/* Drain up to `max` queued events into parallel arrays; returns the count. */
int hil_drain(uint32_t *tick, uint32_t *arg, uint8_t *type, uint8_t *cyl, int max)
{
    int n = 0;
    while (g_tail != g_head && n < max) {
        tick[n] = g_ring[g_tail].tick;
        arg[n] = g_ring[g_tail].arg;
        type[n] = g_ring[g_tail].type;
        cyl[n] = g_ring[g_tail].cyl;
        g_tail = (g_tail + 1u) % HIL_RING;
        ++n;
    }
    return n;
}

/* Test hook: drive the scheduler with fixed advance / pulse (bypasses the slow loop). */
void hil_sched_set(uint32_t cyl, int32_t adv_u, uint32_t pw_ticks)
{
    if (cyl < ECU_NCYL) {
        g_ecu.sch.adv_u[cyl] = adv_u;
        g_ecu.sch.pw_ticks[cyl] = pw_ticks;
    }
}

/* Decoder / scheduler status: [sync, tooth, half, period, rev_ticks, n_sync, n_loss, n_noise,
 * n_dropout, n_cam_err, cam_angle_u, rev_cut, n_cut, n_knock, overflow, limp] */
void hil_status(uint32_t *o)
{
    const ecu_decoder_t *d = &g_ecu.dec;
    o[0] = d->sync;
    o[1] = d->tooth;
    o[2] = d->half;
    o[3] = d->period;
    o[4] = d->rev_ticks;
    o[5] = d->n_sync;
    o[6] = d->n_loss;
    o[7] = d->n_noise;
    o[8] = d->n_dropout;
    o[9] = d->n_cam_err;
    o[10] = (uint32_t)d->cam_angle_u;
    o[11] = g_ecu.sch.rev_cut;
    o[12] = g_ecu.sch.n_cut;
    o[13] = g_ecu.n_knock;
    o[14] = g_overflow;
    o[15] = g_ecu.limp;
}

/* Generated constants the plant needs to agree on (so Python never retypes them). */
void hil_consts(uint32_t *o)
{
    o[0] = ECU_TIMER_HZ;
    o[1] = ECU_ANGLE_U_PER_DEG;
    o[2] = ECU_NCYL;
    o[3] = ECU_TEETH_TOTAL;
    o[4] = ECU_TEETH_MISSING;
    o[5] = (uint32_t)ECU_TOOTH0_U;
    o[6] = (uint32_t)ECU_CAM_EDGE_U;
    o[7] = ECU_KNOCK_THRESHOLD;
    o[8] = ECU_EV_INJ;
    o[9] = ECU_EV_DWELL;
    o[10] = ECU_EV_SPARK;
}
