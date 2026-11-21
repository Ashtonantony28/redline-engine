/* Angle-domain event scheduler. Fixed point only (built with -mgeneral-regs-only).
 *
 * On every tooth with full sync, each cylinder's three angle events (injection start, dwell
 * start, spark) are tested against the window [tooth angle, next real tooth). An event inside the
 * window is armed as an absolute timer compare: tick + delta_u * period / pitch_u, i.e. linear
 * interpolation on the per-pitch period averaged over the last few teeth (Q8 sub-tick). Time-defined quantities (injector pulse, coil dwell)
 * are converted to angle with the same period. A repeat guard stops an event firing twice per cycle
 * when its target moves; an event whose target slipped up to one pitch behind the current tooth fires
 * immediately (late catch-up) instead of being lost. The rev limiter (fuel + spark cut with hysteresis) compares the
 * revolution period against generated tick thresholds.
 */
#include "ecu/ecu.h"
#include "ecu/ecu_hal.h"

enum { EV_INJ = 0, EV_DWELL = 1, EV_SPARK = 2, EV_N = 3 };

void ecu_sched_init(ecu_sched_t *s)
{
    uint8_t *p = (uint8_t *)s;
    for (uint32_t i = 0; i < sizeof(*s); ++i) {
        p[i] = 0u;
    }
}

static uint32_t ticks_to_u(uint32_t ticks, uint32_t period_q8, uint32_t cap)
{
    uint64_t u = (uint64_t)ticks * ECU_PITCH_U * ECU_Q8 / period_q8;
    return u > cap ? cap : (uint32_t)u;
}

void ecu_sched_tooth(ecu_sched_t *s, const ecu_decoder_t *d, const ecu_tooth_t *t)
{
    if (!t->valid || !t->full || t->period_q8 == 0u) {
        return;
    }
    uint32_t cut = s->limp ? ECU_LIMP_CUT_TICKS : ECU_REV_CUT_TICKS;
    uint32_t resume = s->limp ? ECU_LIMP_RESUME_TICKS : ECU_REV_RESUME_TICKS;
    if (d->rev_ticks != 0u) {
        if (d->rev_ticks < cut) {
            s->rev_cut = 1u;
        } else if (d->rev_ticks > resume) {
            s->rev_cut = 0u;
        }
    }
    uint32_t dwell_u = ticks_to_u(ECU_DWELL_TICKS, t->period_q8, ECU_REV_U / 2u);
    for (uint8_t c = 0u; c < ECU_NCYL; ++c) {
        for (uint8_t k = 0u; k < EV_N; ++k) {
            int32_t target;
            uint32_t pw = s->pw_ticks[c];
            if (k == EV_INJ) {
                if (pw == 0u) {
                    continue;
                }
                target = ECU_TDC_U[c] + ECU_EOI_U - (int32_t)ticks_to_u(pw, t->period_q8, ECU_MAX_PW_U);
            } else if (k == EV_DWELL) {
                target = ECU_TDC_U[c] - s->adv_u[c] - (int32_t)dwell_u;
            } else {
                target = ECU_TDC_U[c] - s->adv_u[c];
            }
            uint32_t delta = (uint32_t)ecu_angle_wrap(target - t->angle_u);
            if (delta >= t->span_u) {
                /* late by less than one pitch (target jittered across the window edge, or moved
                 * back): fire now rather than lose the event; the repeat guard blocks doubles */
                if ((uint32_t)ECU_CYCLE_U - delta > (uint32_t)ECU_PITCH_U) {
                    continue;
                }
                delta = 0u;
            }
            if (s->fired_once[c][k] && (d->abs_u - s->last_fire_abs[c][k]) < ECU_REPEAT_GUARD_U) {
                continue;
            }
            s->fired_once[c][k] = 1u;
            s->last_fire_abs[c][k] = d->abs_u;
            if (s->rev_cut && (k == EV_INJ || ECU_REV_SPARK_CUT)) {
                s->n_cut++;
                continue;
            }
            ecu_event_t ev;
            ev.tick = t->tick + (uint32_t)(((uint64_t)delta * t->period_q8 + (uint64_t)ECU_PITCH_U * ECU_Q8 / 2u) /
                                           ((uint64_t)ECU_PITCH_U * ECU_Q8));
            ev.cyl = c;
            if (k == EV_INJ) {
                ev.type = ECU_EV_INJ;
                ev.arg = pw;
                s->inj_pw[c] = pw;
                s->inj_count[c]++;
            } else {
                ev.type = (k == EV_DWELL) ? ECU_EV_DWELL : ECU_EV_SPARK;
                ev.arg = 0u;
            }
            hal_emit_event(&ev);
        }
    }
}
