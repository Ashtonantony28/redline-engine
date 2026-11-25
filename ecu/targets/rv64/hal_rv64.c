/* RV64 HAL + firmware main loop for the Vortex-V SoC (M9, M9.1, M9.2). The portable core in ecu/src
 * is built unchanged; this file provides hal_emit_event() and the platform side of the calling
 * contract in ecu.h. Three modes share it (hil_mbox.h):
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
 *   capture/compare mode (M9.2, ADR 0007; real-time mode + HIL_RT_F_CAPCMP) uses the SoC's
 *   capture/compare timer (Vortex-V ADR 0017) instead of the played mailbox peripheral:
 *     - crank and cam edges are timestamped by the unit's input capture; irq_handler() is the
 *       machine external interrupt (level, while the FIFO is not empty). It drains the FIFO in
 *       order, merged by tick with the knock windows that stay in the mailbox (the unit has no
 *       knock input). No MTIMECMP, no prediction, no spin.
 *     - hal_emit_event() pushes one one-shot compare per event: INJ = pulse of `arg` ticks on
 *       output cyl, DWELL = set and SPARK = clear on output ECU_NCYL + cyl. The slot count is
 *       checked before every push: `free` counts slots known to be free and is decremented per
 *       push; only when it reaches 0 is CMP_STATUS read again (slots only free on their own, so
 *       the count is a lower bound and the read stays off the usual path). An event with no
 *       free slot is held in a backlog, pushed on the next capture and counted: the unit never
 *       drops a push unseen.
 *     - the slow loop commits its outputs with interrupts masked, so a capture sees either the
 *       old or the new advance / pulse width of every cylinder, never a mix.
 *
 * SoC facilities used (register map from $VORTEX_V_ROOT/sw/common/vortex.h, not retyped here):
 *   CYCLE     cycle counter (mtime): handler cost, slow-loop wall time, idle calibration
 *   MTIMECMP  machine timer compare, mie / mstatus (irq_enable, irq_disable)
 *   CC_*      capture/compare timer (M9.2): TIME, PRESC, CTRL, CAP_*, CMP_*, mie.MEIE
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

/* capture/compare timer (Vortex-V ADR 0017, SoC parameters N_OUT = 8, N_CMP = 8) */
#define CC_N_CMP 8u                     /* one-shot compare slots in vortex_capcmp */
#define CC_PENDING_MASK 0xFFu           /* CMP_STATUS bits 7:0: slots busy */
#define CC_DROP_SHIFT 16                /* CMP_STATUS bits 47:16: pushes dropped */
#define CC_DROP_MASK 0xFFFFFFFFu
#define CC_OUT_INJ0 0u                  /* injector of cylinder c on cmp_out[c] */
#define CC_OUT_COIL0 ECU_NCYL           /* coil of cylinder c on cmp_out[ECU_NCYL + c] */
#define CC_KIND_CRANK 1u                /* capture kind = channel + 1: cap_in[0] crank */
#define CC_KIND_CAM 2u                  /*                              cap_in[1] cam */
#define CC_CTRL_EDGES ((CC_EDGE_RISE << 0) | (CC_EDGE_RISE << 2))   /* rising edge on both channels */
#define CC_BACKLOG 8u                   /* events held back for want of a slot (power of 2) */

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
    uint8_t capcmp;                     /* M9.2: the SoC capture/compare timer */
} g_rt;

/* capture/compare mode state (M9.2) */
static struct {
    uint32_t free;                      /* slots known free (a lower bound: slots free on their own) */
    uint64_t cap_rd;                    /* shadow of CC_CAP_RD */
    uint32_t bl_rd, bl_wr;              /* backlog ring */
    ecu_event_t backlog[CC_BACKLOG];
    uint64_t n_push, n_deferred, n_spurious;
} g_cc;
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

/* One compare per event (ADR 0007; mapping as in the Vortex-V register map for Redline). */
static void cc_push_event(const ecu_event_t *ev)
{
    if (ev->type == ECU_EV_INJ) {
        cc_pulse(ev->tick, CC_OUT_INJ0 + ev->cyl, ev->arg);
    } else {
        cc_push(ev->tick, CC_OUT_COIL0 + ev->cyl, ev->type == ECU_EV_DWELL ? CC_ACT_SET : CC_ACT_CLEAR);
    }
    g_cc.n_push++;
}

/* Arm one compare, or hold the event back when every slot is busy. */
static void cc_emit(const ecu_event_t *ev)
{
    if (g_cc.free == 0u) {
        g_cc.free = CC_N_CMP - (uint32_t)(CC_CMP_STATUS & CC_PENDING_MASK);
    }
    if (g_cc.free == 0u || g_cc.bl_rd != g_cc.bl_wr) {    /* keep arming order behind a backlog */
        if (g_cc.bl_wr - g_cc.bl_rd < CC_BACKLOG) {
            g_cc.backlog[g_cc.bl_wr++ % CC_BACKLOG] = *ev;
        }
        g_cc.n_deferred++;                  /* counted either way; a full backlog loses the event */
        return;
    }
    g_cc.free--;
    cc_push_event(ev);
}

static void cc_backlog_retry(void)
{
    while (g_cc.bl_rd != g_cc.bl_wr) {
        g_cc.free = CC_N_CMP - (uint32_t)(CC_CMP_STATUS & CC_PENDING_MASK);
        if (g_cc.free == 0u) {
            return;
        }
        cc_push_event(&g_cc.backlog[g_cc.bl_rd++ % CC_BACKLOG]);
    }
}

/* Output compare. Capture/compare mode: push a compare slot. Real-time mode: arm it (store the
 * event, then move the write index; the model stamps the cycle). Command mode: queue the event for
 * the plant (drained after the acknowledge). */
void hal_emit_event(const ecu_event_t *ev)
{
    if (g_rt.capcmp) {
        cc_emit(ev);
        return;
    }
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
    l[HIL_L_RT_T0] = (uint32_t)offsetof(hil_mbox_t, rt.t0);
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

/* The ECU timer: the SoC capture/compare timer in M9.2, the model-played mailbox word otherwise. */
static inline uint32_t rt_now(void)
{
    return g_rt.capcmp ? (uint32_t)CC_TIME : g_mbox.rt.timer;
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

/* Capture/compare mode (M9.2): the machine external interrupt. Drains the hardware capture FIFO
 * and the mailbox knock FIFO in tick order (a knock window strictly older than the oldest
 * hardware capture goes first; the host build sees them in that order too). Every crank / cam
 * entry is popped after it has been serviced, so the pop stamps the end of its service. */
static void cc_isr(void)
{
    const uint64_t t0 = CYCLE;
    uint32_t n = 0u;
    cc_backlog_retry();
    for (;;) {
        const uint64_t e = CC_CAP_DATA;
        const uint32_t rd = g_mbox.rt.cap_rd;
        if (rd != g_mbox.rt.cap_wr) {
            const volatile hil_cap_t *c = &g_mbox.rt.cap[rd % HIL_CAP_SLOTS];
            if (!CC_CAP_VALID(e) || (int32_t)(c->reserved - CC_CAP_TICK(e)) < 0) {
                if (c->kind == HIL_CAP_KNOCK) {
                    ecu_knock_window(&g_ecu, (uint8_t)c->a0, c->a1);
                }
                g_serviced = g_serviced + 1u;
                g_mbox.rt.cap_rd = rd + 1u;
                ++n;
                continue;
            }
        }
        if (!CC_CAP_VALID(e)) {
            break;
        }
        if (CC_CAP_KIND(e) == CC_KIND_CRANK) {
            ecu_crank_edge(&g_ecu, CC_CAP_TICK(e));
        } else if (CC_CAP_KIND(e) == CC_KIND_CAM) {
            ecu_cam_edge(&g_ecu, CC_CAP_TICK(e));
        }
        g_serviced = g_serviced + 1u;
        CC_CAP_RD = ++g_cc.cap_rd;          /* pop: the interrupt drops once the FIFO is empty */
        ++n;
    }
    if (n == 0u) {
        g_cc.n_spurious++;                  /* re-taken before the last pop reached the unit */
    }
    const uint64_t t1 = CYCLE;
    g_rt.n_irq++;
    g_rt.isr_sum += t1 - t0;
    if (t1 - t0 > g_rt.isr_max) {
        g_rt.isr_max = t1 - t0;
    }
}

/* Interrupt entry (Vortex-V crt0: caller-saved registers stacked, mstatus.MIE clear).
 * M9.1: the machine timer interrupt; must rewrite MTIMECMP before returning.
 * M9.2: the machine external interrupt of the capture/compare timer. */
void irq_handler(uint64_t mcause)
{
    if (g_rt.capcmp) {
        if (mcause == MCAUSE_EXT) {
            cc_isr();
        } else {
            timer_disarm();                 /* not used in this mode */
        }
        return;
    }
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
    in.now_tick = rt_now() + NOW_TICK_GUARD;                  /* not older than any snapshot tooth */
    const uint8_t had_last = g_slow.dec.have_last;
    ecu_slow_step(&g_slow, &in, &out);

    /* commit: single-word stores, each atomic against the handler; M9.2 masks the interrupt over
     * the whole commit, so no capture sees one cylinder's new advance with its old pulse width */
    if (g_rt.capcmp) {
        irq_disable();
    }
    for (uint32_t c = 0u; c < ECU_NCYL; ++c) {
        g_ecu.sch.adv_u[c] = g_slow.sch.adv_u[c];
        g_ecu.sch.pw_ticks[c] = g_slow.sch.pw_ticks[c];
    }
    g_ecu.sch.limp = g_slow.sch.limp;
    if (had_last && !g_slow.dec.have_last) {                  /* the step declared a stall */
        if (g_rt.irq && !g_rt.capcmp) {
            irq_disable();
        }
        if (g_ecu.dec.n_edges == g_slow.dec.n_edges) {        /* still true: no edge since the snapshot */
            ecu_dec_stall(&g_ecu.dec);
            st[HIL_RT_S_STALLS]++;
        }
        if (g_rt.irq && !g_rt.capcmp) {
            irq_enable();
        }
    }
    if (g_rt.capcmp && g_rt.irq) {
        CSR_SET(0x300, MSTATUS_MIE);        /* mie.MEIE stays set from rt_run */
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
        if ((int32_t)(rt_now() - next_slow) >= 0) {
            rt_slow();
            next_slow += SLOW_TICKS;
            while ((int32_t)(rt_now() - next_slow) >= 0) {            /* overran: skip, stay on the grid */
                next_slow += SLOW_TICKS;
                g_mbox.rt.stats[HIL_RT_S_SLOW_OVERRUNS]++;
            }
            continue;
        }
        if (!g_rt.irq && !g_rt.capcmp && g_mbox.rt.cap_rd != g_mbox.rt.cap_wr) {
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
    g_rt.capcmp = (uint8_t)((flags & HIL_RT_F_CAPCMP) != 0u);
    if (g_rt.capcmp) {                      /* unit to a known state; the timer runs at ECU_TIMER_HZ */
        CC_CTRL = 0u;
        CC_CMP_FLUSH = 0u;
        CC_OUT = 0u;
        CC_PRESC = rt->cfg[HIL_RT_CLK_HZ] / ECU_TIMER_HZ - 1u;
        g_cc.cap_rd = CC_CAP_WR;
        CC_CAP_RD = g_cc.cap_rd;            /* pop everything */
        CC_CAP_OVF = 0u;
        CC_CMP_STATUS = 0u;                 /* clear the drop count */
        CC_CTRL = CC_CTRL_EDGES;            /* capturing, interrupt still off */
        g_cc.free = 0u;
        g_cc.bl_rd = g_cc.bl_wr = 0u;
        g_cc.n_push = g_cc.n_deferred = g_cc.n_spurious = 0u;
    }
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
    uint32_t first_slow = 0u;
    if (g_rt.capcmp) {
        first_slow = rt->t0;                /* the plant's time 0 on the hardware timer */
        if (g_rt.irq) {
            timer_disarm();
            CC_CTRL = CC_CTRL_EDGES | CC_CTRL_CAP_IE;
            CSR_SET(0x304, MIE_MEIE);
            CSR_CLEAR(0x304, MIE_MTIE);
            CSR_SET(0x300, MSTATUS_MIE);
        }
    } else if (g_rt.irq) {
        MTIMECMP = t0 + g_rt.poll;
        irq_enable();
    }
    const uint64_t idle = rt_background(first_slow, UINT64_MAX);
    if (g_rt.irq) {
        irq_disable();
        timer_disarm();
    }
    if (g_rt.capcmp) {
        CSR_CLEAR(0x304, MIE_MEIE);
        CC_CTRL = CC_CTRL_EDGES;
        rt->stats[HIL_RT_S_CC_PUSHES] = g_cc.n_push;
        rt->stats[HIL_RT_S_CC_DEFERRED] = g_cc.n_deferred;
        rt->stats[HIL_RT_S_CC_DROPS] = (CC_CMP_STATUS >> CC_DROP_SHIFT) & CC_DROP_MASK;
        rt->stats[HIL_RT_S_CC_CAP_OVF] = CC_CAP_OVF;
        rt->stats[HIL_RT_S_CC_SPURIOUS] = g_cc.n_spurious;
        g_rt.capcmp = 0u;                   /* command mode: events go to the mailbox again */
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
