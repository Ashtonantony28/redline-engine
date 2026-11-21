/* RV64 HAL + firmware main loop for the Vortex-V SoC (M9). The portable core in ecu/src is built
 * unchanged; this file provides hal_emit_event() and the platform side of the calling contract in
 * ecu.h, with the capture events arriving through the mailbox (hil_mbox.h) instead of interrupts.
 *
 * SoC facilities used (register map from $VORTEX_V_ROOT/sw/common/vortex.h, not retyped here):
 *   CYCLE   cycle counter, read around every handler to report its cost
 *   LEDS    heartbeat: the low bits of the command count
 * No heap, no libc: memset / memcpy below serve the compiler's own struct copies. */
#include <stddef.h>
#include <stdint.h>

#include "vortex.h"

#include "ecu/ecu.h"
#include "ecu/ecu_hal.h"
#include "hil_mbox.h"
#include "knock_kernel.h"

#define LED_MASK 0xFu

volatile hil_mbox_t g_mbox;
static ecu_t g_ecu;
static uint32_t g_ncmd;

void *memset(void *dst, int c, size_t n)
{
    uint8_t *d = (uint8_t *)dst;
    while (n--) {
        *d++ = (uint8_t)c;
    }
    return dst;
}

void *memcpy(void *dst, const void *src, size_t n)
{
    uint8_t *d = (uint8_t *)dst;
    const uint8_t *s = (const uint8_t *)src;
    while (n--) {
        *d++ = *s++;
    }
    return dst;
}

/* Output compare: queue the event for the plant (drained after the acknowledge). */
void hal_emit_event(const ecu_event_t *ev)
{
    uint32_t n = g_mbox.n_events;
    if (n >= HIL_MBOX_EVENTS) {
        g_mbox.overflow++;
        return;
    }
    volatile ecu_event_t *slot = &g_mbox.events[n];
    slot->tick = ev->tick;
    slot->arg = ev->arg;
    slot->type = ev->type;
    slot->cyl = ev->cyl;
    g_mbox.n_events = n + 1u;
}

/* Same order as hil_status() / hil_consts() in ecu/host/hal_host.c (hil/bridge.py names them). */
static void fill_status(void)
{
    const ecu_decoder_t *d = &g_ecu.dec;
    volatile uint32_t *o = g_mbox.status;
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
    o[14] = g_mbox.overflow;
    o[15] = g_ecu.limp;
}

static void fill_consts(void)
{
    volatile uint32_t *o = g_mbox.consts;
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

static void fill_layout(void)
{
    volatile uint32_t *l = g_mbox.layout;
    l[HIL_L_CMD] = (uint32_t)offsetof(hil_mbox_t, cmd);
    l[HIL_L_ARG] = (uint32_t)offsetof(hil_mbox_t, arg);
    l[HIL_L_N_EVENTS] = (uint32_t)offsetof(hil_mbox_t, n_events);
    l[HIL_L_OVERFLOW] = (uint32_t)offsetof(hil_mbox_t, overflow);
    l[HIL_L_CYCLES] = (uint32_t)offsetof(hil_mbox_t, cycles);
    l[HIL_L_EVENTS] = (uint32_t)offsetof(hil_mbox_t, events);
    l[HIL_L_EVENT_SIZE] = (uint32_t)sizeof(ecu_event_t);
    l[HIL_L_EVENT_TICK] = (uint32_t)offsetof(ecu_event_t, tick);
    l[HIL_L_EVENT_ARG] = (uint32_t)offsetof(ecu_event_t, arg);
    l[HIL_L_EVENT_TYPE] = (uint32_t)offsetof(ecu_event_t, type);
    l[HIL_L_EVENT_CYL] = (uint32_t)offsetof(ecu_event_t, cyl);
    l[HIL_L_IN] = (uint32_t)offsetof(hil_mbox_t, in);
    l[HIL_L_IN_SIZE] = (uint32_t)sizeof(ecu_inputs_t);
    l[HIL_L_OUT] = (uint32_t)offsetof(hil_mbox_t, out);
    l[HIL_L_OUT_SIZE] = (uint32_t)sizeof(ecu_outputs_t);
    l[HIL_L_STATUS] = (uint32_t)offsetof(hil_mbox_t, status);
    l[HIL_L_CONSTS] = (uint32_t)offsetof(hil_mbox_t, consts);
    l[HIL_L_KNOCK_X] = (uint32_t)offsetof(hil_mbox_t, knock_x);
    l[HIL_L_KNOCK_N_IN] = KNOCK_N_IN;
    l[HIL_L_KNOCK_LOGIT] = (uint32_t)offsetof(hil_mbox_t, knock_logit);
}

static void slow_step(void)
{
    ecu_inputs_t in;
    ecu_outputs_t out;
    in = *(const ecu_inputs_t *)&g_mbox.in;
    ecu_slow_step(&g_ecu, &in, &out);
    *(ecu_outputs_t *)&g_mbox.out = out;
}

int main(void)
{
    _Static_assert(KNOCK_N_IN <= HIL_KNOCK_MAX_IN, "knock window larger than the mailbox buffer");
    ecu_init(&g_ecu);
    g_mbox.cmd = HIL_CMD_NONE;
    g_mbox.n_events = 0u;
    g_mbox.overflow = 0u;
    fill_layout();
    g_mbox.magic = HIL_MBOX_MAGIC;
    for (;;) {
        uint32_t cmd = g_mbox.cmd;
        if (cmd == HIL_CMD_NONE) {
            continue;
        }
        g_mbox.n_events = 0u;
        uint32_t a0 = g_mbox.arg[0], a1 = g_mbox.arg[1], a2 = g_mbox.arg[2];
        uint64_t c0 = CYCLE;
        switch (cmd) {
        case HIL_CMD_RESET:
            ecu_init(&g_ecu);
            g_mbox.overflow = 0u;
            break;
        case HIL_CMD_CRANK:
            ecu_crank_edge(&g_ecu, a0);
            break;
        case HIL_CMD_CAM:
            ecu_cam_edge(&g_ecu, a0);
            break;
        case HIL_CMD_KNOCK:
            ecu_knock_window(&g_ecu, (uint8_t)a0, a1);
            break;
        case HIL_CMD_SLOW:
            slow_step();
            break;
        case HIL_CMD_SCHED_SET:
            if (a0 < ECU_NCYL) {
                g_ecu.sch.adv_u[a0] = (int32_t)a1;
                g_ecu.sch.pw_ticks[a0] = a2;
            }
            break;
        case HIL_CMD_STATUS:
            fill_status();
            break;
        case HIL_CMD_CONSTS:
            fill_consts();
            break;
        case HIL_CMD_KNOCK_INFER:
            g_mbox.knock_logit = knock_infer((const int8_t *)g_mbox.knock_x);
            break;
        case HIL_CMD_EXIT:
            return 0;                       /* crt0 stores PASS to tohost */
        default:
            break;
        }
        g_mbox.cycles = CYCLE - c0;
        LEDS = ++g_ncmd & LED_MASK;
        g_mbox.cmd = HIL_CMD_NONE;          /* acknowledge: last store of the command */
    }
}
