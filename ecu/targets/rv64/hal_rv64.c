/* RV64 HAL + firmware main loop for the Vortex-V SoC (M9, M9.1). The portable core in ecu/src is
 * built unchanged; this file provides hal_emit_event() and the platform side of the calling
 * contract in ecu.h. Two modes share it (hil_mbox.h):
 *
 *   frozen-time command mode (M9)   every entry point is a polled mailbox command; no interrupts.
 *   real-time mode (M9.1, ADR 0005) the sync path preempts the slow loop:
 *     - irq_handler() (machine timer interrupt, Vortex-V ADR 0015) drains the capture FIFO and
 *       runs ecu_crank_edge / ecu_cam_edge / ecu_knock_window. It then arms MTIMECMP for the next
 *       predicted crank edge minus a lead, waits a bounded time inside the handler if it woke
 *       ahead of the edge, and falls back to a fixed poll period when no edge can be predicted.
 *     - the background loop runs ecu_slow_step once per ECU_DT_SLOW on a snapshot: the
 *       handler-owned state (decoder, scheduler, knock retard) is copied into a second ecu_t
 *       with a retry-if-preempted read, the step runs on the copy with interrupts enabled, and
 *       its outputs (advance, pulse width, limp flag, stall verdict) are committed to the live
 *       instance afterwards. The handler therefore never sees a half-updated slow loop and the
 *       slow loop never sees a decoder that changed under it (a tooth newer than `now_tick`
 *       would read as a stall).
 *
 * SoC facilities used (register map from $VORTEX_V_ROOT/sw/common/vortex.h, not retyped here):
 *   CYCLE     cycle counter (mtime): handler cost, slow-loop wall time, idle calibration
 *   MTIMECMP  machine timer compare, mie / mstatus (irq_enable, irq_disable)
 *   LEDS      heartbeat: the low bits of the command count
 * No heap, no libc: memset / memcpy below serve the compiler's own struct copies. */
#include <stddef.h>
#include <stdint.h>

#include "vortex.h"

#include "ecu/ecu.h"
#include "ecu/ecu_hal.h"
#include "hil_mbox.h"
#include "knock_kernel.h"

#define LED_MASK 0xFu

#define Q16 16                          /* fixed point of the cycles-per-tick factor */
#define SLOW_TICKS ((uint32_t)(ECU_DT_SLOW * ECU_TIMER_HZ_F + 0.5f))    /* slow-loop period [ticks] */
#define SLOW_FAR 0x40000000u            /* "not due" for the idle calibration */
#define IDLE_CAL_ITERS 2048u
#define NOW_TICK_GUARD 1u               /* a captured tick is rounded, the timer truncated */
#define Q8_SHIFT 8                      /* log2(ECU_Q8) */
#define Q8_HALF (ECU_Q8 / 2u)

typedef uint32_t __attribute__((may_alias)) word_alias_t;

volatile hil_mbox_t g_mbox;
static ecu_t g_ecu;                     /* live instance: the capture context owns dec, sch, knock */
static ecu_t g_slow;                    /* real-time mode: the slow loop's snapshot and its own state */
static uint32_t g_ncmd;

/* real-time mode state */
static struct {
    uint8_t active;                     /* hal_emit_event arms rt output compares */
    uint8_t irq;                        /* captures serviced from the timer interrupt */
    uint8_t predict;
    uint8_t expecting;                  /* the timer was armed for a predicted edge */
    uint64_t cyc_per_tick_q16;
    uint64_t poll, lead, spin;          /* [cycles], from rt.cfg */
    uint64_t n_irq, isr_sum, isr_max, n_spin_timeout;
} g_rt;
static volatile uint32_t g_serviced;    /* captures serviced: the snapshot's sequence number */

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

/* Output compare. Real-time mode: arm it (store the event, then move the write index; the model
 * stamps the cycle). Command mode: queue the event for the plant (drained after the acknowledge). */
void hal_emit_event(const ecu_event_t *ev)
{
    if (g_rt.active) {
        const uint32_t w = g_mbox.rt.evt_wr;
        volatile ecu_event_t *slot = &g_mbox.events[w % HIL_MBOX_EVENTS];
        slot->tick = ev->tick;
        slot->arg = ev->arg;
        slot->type = ev->type;
        slot->cyl = ev->cyl;
        g_mbox.rt.evt_wr = w + 1u;
        return;
    }
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
    l[HIL_L_RT_CFG] = (uint32_t)offsetof(hil_mbox_t, rt.cfg);
    l[HIL_L_RT_READY] = (uint32_t)offsetof(hil_mbox_t, rt.ready);
    l[HIL_L_RT_GO] = (uint32_t)offsetof(hil_mbox_t, rt.go);
    l[HIL_L_RT_STOP] = (uint32_t)offsetof(hil_mbox_t, rt.stop);
    l[HIL_L_RT_TIMER] = (uint32_t)offsetof(hil_mbox_t, rt.timer);
    l[HIL_L_RT_CAP_WR] = (uint32_t)offsetof(hil_mbox_t, rt.cap_wr);
    l[HIL_L_RT_CAP_RD] = (uint32_t)offsetof(hil_mbox_t, rt.cap_rd);
    l[HIL_L_RT_EVT_WR] = (uint32_t)offsetof(hil_mbox_t, rt.evt_wr);
    l[HIL_L_RT_IN_SEL] = (uint32_t)offsetof(hil_mbox_t, rt.in_sel);
    l[HIL_L_RT_OUT_SEL] = (uint32_t)offsetof(hil_mbox_t, rt.out_sel);
    l[HIL_L_RT_SLOW_SEQ] = (uint32_t)offsetof(hil_mbox_t, rt.slow_seq);
    l[HIL_L_RT_CAP] = (uint32_t)offsetof(hil_mbox_t, rt.cap);
    l[HIL_L_RT_CAP_SIZE] = (uint32_t)sizeof(hil_cap_t);
    l[HIL_L_RT_IN] = (uint32_t)offsetof(hil_mbox_t, rt.in);
    l[HIL_L_RT_OUT] = (uint32_t)offsetof(hil_mbox_t, rt.out);
    l[HIL_L_RT_STATS] = (uint32_t)offsetof(hil_mbox_t, rt.stats);
}

static void slow_step(void)
{
    ecu_inputs_t in;
    ecu_outputs_t out;
    in = *(const ecu_inputs_t *)&g_mbox.in;
    ecu_slow_step(&g_ecu, &in, &out);
    *(ecu_outputs_t *)&g_mbox.out = out;
}

/* ---------------------------------------------------------------------------------------------
 * real-time mode (M9.1, ADR 0005)
 * ------------------------------------------------------------------------------------------- */

/* Word copy through volatile lvalues: the source or the destination changes under an interrupt
 * or under the plant. Both are 4-byte aligned and a multiple of 4 bytes long. */
static void copy_words(volatile void *dst, const volatile void *src, size_t bytes)
{
    volatile word_alias_t *d = (volatile word_alias_t *)dst;
    const volatile word_alias_t *s = (const volatile word_alias_t *)src;
    for (size_t i = 0; i < bytes / sizeof(word_alias_t); ++i) {
        d[i] = s[i];
    }
}

/* Drain the capture FIFO in arrival order. Returns the number of entries serviced. */
static uint32_t rt_service(void)
{
    uint32_t rd = g_mbox.rt.cap_rd;
    uint32_t n = 0u;
    while (rd != g_mbox.rt.cap_wr) {
        const volatile hil_cap_t *c = &g_mbox.rt.cap[rd % HIL_CAP_SLOTS];
        const uint32_t kind = c->kind, a0 = c->a0;
        if (kind == HIL_CAP_CRANK) {
            ecu_crank_edge(&g_ecu, a0);
        } else if (kind == HIL_CAP_CAM) {
            ecu_cam_edge(&g_ecu, a0);
        } else if (kind == HIL_CAP_KNOCK) {
            ecu_knock_window(&g_ecu, (uint8_t)a0, c->a1);
        }
        g_serviced = g_serviced + 1u;
        g_mbox.rt.cap_rd = ++rd;            /* serviced: the model stamps this store */
        ++n;
    }
    return n;
}

/* Cycles from now until the timer should fire next. With a period estimate the next crank edge is
 * expected one pitch (or the gap) after the last accepted one; wake `lead` cycles before that.
 * An edge overdue by more than the spin window is not waited for: poll. */
static uint64_t rt_next_wake(void)
{
    const ecu_decoder_t *d = &g_ecu.dec;
    g_rt.expecting = 0u;
    if (!g_rt.predict || !d->have_last || d->period_q8 == 0u) {
        return g_rt.poll;
    }
    const uint32_t steps = (d->sync != ECU_SYNC_NONE && d->tooth == ECU_TEETH_REAL - 1) ? ECU_TEETH_MISSING + 1u : 1u;
    const uint32_t next_tick = d->last_tick + ((d->period_q8 * steps + Q8_HALF) >> Q8_SHIFT);
    const int64_t ahead = (int64_t)(int32_t)(next_tick - g_mbox.rt.timer) * (int64_t)g_rt.cyc_per_tick_q16 / (1 << Q16);
    if (ahead < -(int64_t)g_rt.spin) {
        return g_rt.poll;
    }
    /* no poll in between: cam and knock captures wait in the FIFO, in order, for the tooth's wake */
    const int64_t wake = ahead - (int64_t)g_rt.lead;
    g_rt.expecting = 1u;
    return wake > 0 ? (uint64_t)wake : 0u;
}

/* Machine timer interrupt (called from Vortex-V's crt0 with the caller-saved registers stacked,
 * mstatus.MIE clear). Must rewrite MTIMECMP before returning. */
void irq_handler(void)
{
    const uint64_t t0 = CYCLE;
    uint32_t n = rt_service();
    if (n == 0u && g_rt.expecting) {        /* woke ahead of the predicted edge: wait for it */
        const uint64_t limit = t0 + g_rt.spin;
        while (g_mbox.rt.cap_wr == g_mbox.rt.cap_rd && CYCLE < limit) {
        }
        n = rt_service();
        if (n == 0u) {
            g_rt.n_spin_timeout++;
        }
    }
    const uint64_t wake = rt_next_wake();
    const uint64_t t1 = CYCLE;
    MTIMECMP = t1 + wake;
    g_rt.n_irq++;
    g_rt.isr_sum += t1 - t0;
    if (t1 - t0 > g_rt.isr_max) {
        g_rt.isr_max = t1 - t0;
    }
}

/* One slow-loop step on a snapshot (see the file header). */
static void rt_slow(void)
{
    volatile uint64_t *st = g_mbox.rt.stats;
    const uint64_t t0 = CYCLE;
    const size_t handler_owned = offsetof(ecu_t, stft_i);     /* dec, sch, knock_retard_u, n_knock */
    uint32_t seq;
    for (;;) {
        seq = g_serviced;
        copy_words(&g_slow, &g_ecu, handler_owned);
        if (seq == g_serviced) {
            break;
        }
        st[HIL_RT_S_SNAP_RETRIES]++;
    }
    ecu_inputs_t in;
    ecu_outputs_t out;
    copy_words(&in, &g_mbox.rt.in[g_mbox.rt.in_sel & 1u], sizeof(in));
    in.now_tick = g_mbox.rt.timer + NOW_TICK_GUARD;           /* not older than any snapshot tooth */
    const uint8_t had_last = g_slow.dec.have_last;
    ecu_slow_step(&g_slow, &in, &out);

    /* commit: single-word stores, each atomic against the handler */
    for (uint32_t c = 0u; c < ECU_NCYL; ++c) {
        g_ecu.sch.adv_u[c] = g_slow.sch.adv_u[c];
        g_ecu.sch.pw_ticks[c] = g_slow.sch.pw_ticks[c];
    }
    g_ecu.sch.limp = g_slow.sch.limp;
    if (had_last && !g_slow.dec.have_last) {                  /* the step declared a stall */
        if (g_rt.irq) {
            irq_disable();
        }
        if (g_ecu.dec.n_edges == g_slow.dec.n_edges) {        /* still true: no edge since the snapshot */
            ecu_dec_stall(&g_ecu.dec);
            st[HIL_RT_S_STALLS]++;
        }
        if (g_rt.irq) {
            irq_enable();
        }
    }
    const uint32_t sel = (g_mbox.rt.out_sel ^ 1u) & 1u;
    copy_words(&g_mbox.rt.out[sel], &out, sizeof(out));
    g_mbox.rt.out_sel = sel;
    g_mbox.rt.slow_seq = g_mbox.rt.slow_seq + 1u;
    const uint64_t wall = CYCLE - t0;
    st[HIL_RT_S_SLOW_STEPS]++;
    st[HIL_RT_S_SLOW_WALL_SUM] += wall;
    if (wall > st[HIL_RT_S_SLOW_WALL_MAX]) {
        st[HIL_RT_S_SLOW_WALL_MAX] = wall;
    }
}

/* Background loop: the slow loop when its period is due, captures too in polled mode, else idle.
 * Returns the idle passes made; ends on rt.stop or after `max_idle` idle passes (calibration). */
static uint64_t rt_background(uint32_t next_slow, uint64_t max_idle)
{
    uint64_t idle = 0u;
    while (!g_mbox.rt.stop && idle < max_idle) {
        if ((int32_t)(g_mbox.rt.timer - next_slow) >= 0) {
            rt_slow();
            next_slow += SLOW_TICKS;
            while ((int32_t)(g_mbox.rt.timer - next_slow) >= 0) {     /* overran: skip, stay on the grid */
                next_slow += SLOW_TICKS;
                g_mbox.rt.stats[HIL_RT_S_SLOW_OVERRUNS]++;
            }
            continue;
        }
        if (!g_rt.irq && g_mbox.rt.cap_rd != g_mbox.rt.cap_wr) {
            rt_service();
            continue;
        }
        ++idle;
    }
    return idle;
}

static void rt_run(void)
{
    volatile hil_rt_t *rt = &g_mbox.rt;
    const uint32_t flags = rt->cfg[HIL_RT_FLAGS];
    const size_t handler_owned = offsetof(ecu_t, stft_i);
    g_rt.irq = (uint8_t)((flags & HIL_RT_F_IRQ) != 0u);
    g_rt.predict = (uint8_t)((flags & HIL_RT_F_PREDICT) != 0u);
    g_rt.expecting = 0u;
    g_rt.cyc_per_tick_q16 = ((uint64_t)rt->cfg[HIL_RT_CLK_HZ] << Q16) / ECU_TIMER_HZ;
    g_rt.poll = rt->cfg[HIL_RT_POLL_CYCLES];
    g_rt.lead = rt->cfg[HIL_RT_LEAD_CYCLES];
    g_rt.spin = rt->cfg[HIL_RT_SPIN_CYCLES];
    g_rt.n_irq = g_rt.isr_sum = g_rt.isr_max = g_rt.n_spin_timeout = 0u;
    for (uint32_t i = 0u; i < HIL_RT_NSTATS; ++i) {
        rt->stats[i] = 0u;
    }
    rt->timer = 0u;
    rt->cap_rd = rt->cap_wr;
    rt->stop = 0u;
    g_slow = g_ecu;                         /* the slow loop's own state moves to its instance */

    /* idle calibration: nothing due, nothing pending, interrupts off */
    const uint64_t c0 = CYCLE;
    (void)rt_background(SLOW_FAR, IDLE_CAL_ITERS);
    rt->stats[HIL_RT_S_IDLE_CAL_CYCLES] = CYCLE - c0;
    rt->stats[HIL_RT_S_IDLE_CAL_ITERS] = IDLE_CAL_ITERS;

    rt->go = 0u;
    rt->ready = 1u;
    while (!rt->go) {
    }
    const uint64_t t0 = CYCLE;
    g_rt.active = 1u;
    if (g_rt.irq) {
        MTIMECMP = t0 + g_rt.poll;
        irq_enable();
    }
    const uint64_t idle = rt_background(0u, UINT64_MAX);
    if (g_rt.irq) {
        irq_disable();
        timer_disarm();
    }
    g_rt.active = 0u;
    rt->stats[HIL_RT_S_RUN_CYCLES] = CYCLE - t0;
    rt->stats[HIL_RT_S_IDLE_ITERS] = idle;
    rt->stats[HIL_RT_S_IRQS] = g_rt.n_irq;
    rt->stats[HIL_RT_S_ISR_CYCLES] = g_rt.isr_sum;
    rt->stats[HIL_RT_S_ISR_MAX] = g_rt.isr_max;
    rt->stats[HIL_RT_S_SPIN_TIMEOUTS] = g_rt.n_spin_timeout;
    rt->stats[HIL_RT_S_SERVICED] = g_serviced;
    /* back to one instance: the slow loop's state returns to the live one */
    copy_words((uint8_t *)&g_ecu + handler_owned, (const uint8_t *)&g_slow + handler_owned,
               sizeof(ecu_t) - handler_owned);
    rt->ready = 0u;
}

int main(void)
{
    _Static_assert(offsetof(ecu_t, stft_i) % sizeof(word_alias_t) == 0 && sizeof(ecu_t) % sizeof(word_alias_t) == 0 &&
                   sizeof(ecu_inputs_t) % sizeof(word_alias_t) == 0 && sizeof(ecu_outputs_t) % sizeof(word_alias_t) == 0,
                   "copy_words needs whole words");
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
        case HIL_CMD_RT:
            rt_run();
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
