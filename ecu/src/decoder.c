/* 60-2 crank decoder + single-edge cam phase. Fixed point only (built with -mgeneral-regs-only).
 *
 * Every accepted edge is classified against the per-pitch period estimate:
 *   dt <  noise_ratio   x period : noise, rejected (edge ignored, period kept)
 *   dt >= gap_ratio     x period : the missing-tooth gap (missing + 1 pitches)
 *   dt >= dropout_ratio x period : one tooth lost (2 pitches), counted, index advanced by 2
 *   otherwise                    : normal tooth
 * Two gaps in a row are impossible on the wheel: the second is re-classified as a normal tooth so a
 * period estimate corrupted by start-up noise re-baselines instead of locking in at 1/3.
 * While searching, a gap is only believed after ECU_SYNC_MIN_GOOD normal teeth.
 * Sync is declared on the first such gap (crank angle known within one revolution of the gap); a gap
 * at the wrong tooth count or no gap where expected is a loss of sync. A gap at the wrong count
 * resynchronises on that gap immediately; a missed gap drops to search. Cam phase (720 deg) comes
 * from the cam edge: of the two candidate revolution halves, the one that puts the edge nearest
 * ECU_CAM_EDGE_U wins.
 */
#include "ecu/ecu.h"

enum { CLS_NORMAL, CLS_GAP, CLS_DROPOUT, CLS_NOISE };

#define GOOD_SAT 255u

int32_t ecu_angle_wrap(int32_t a)
{
    a %= (int32_t)ECU_CYCLE_U;
    return a < 0 ? a + (int32_t)ECU_CYCLE_U : a;
}

void ecu_dec_init(ecu_decoder_t *d)
{
    uint8_t *p = (uint8_t *)d;
    for (uint32_t i = 0; i < sizeof(*d); ++i) {
        p[i] = 0u;
    }
}

static void lose_sync(ecu_decoder_t *d)
{
    d->sync = ECU_SYNC_NONE;
    d->rev_seen = 0u;
    d->rev_ticks = 0u;
}

void ecu_dec_stall(ecu_decoder_t *d)
{
    if (d->sync != ECU_SYNC_NONE) {
        d->n_loss++;
    }
    lose_sync(d);
    d->have_last = 0u;
    d->have_period = 0u;
    d->period = 0u;
    d->period_q8 = 0u;
    d->hist_n = 0u;
}

static int classify(uint32_t dt, uint32_t period)
{
    uint64_t a = (uint64_t)dt * ECU_RATIO_DEN;
    uint64_t p = (uint64_t)period;
    if (a >= p * ECU_GAP_RATIO_NUM) {
        return CLS_GAP;
    }
    if (a >= p * ECU_DROPOUT_RATIO_NUM) {
        return CLS_DROPOUT;
    }
    if (a < p * ECU_NOISE_RATIO_NUM) {
        return CLS_NOISE;
    }
    return CLS_NORMAL;
}

static int32_t tooth_angle(uint8_t tooth, uint8_t half)
{
    return ecu_angle_wrap((int32_t)ECU_TOOTH0_U + (int32_t)tooth * ECU_PITCH_U + (int32_t)half * (int32_t)ECU_REV_U);
}

void ecu_dec_crank(ecu_decoder_t *d, uint32_t tick, ecu_tooth_t *out)
{
    out->valid = 0u;
    out->full = 0u;
    d->n_edges++;
    if (!d->have_last) {
        d->have_last = 1u;
        d->last_tick = tick;
        return;
    }
    uint32_t dt = tick - d->last_tick;       /* wrap-safe */
    if (!d->have_period) {
        d->have_period = 1u;
        d->period = dt;
        d->last_tick = tick;
        return;
    }
    int cls = classify(dt, d->period);
    if (cls == CLS_NOISE) {
        d->n_noise++;
        return;
    }
    if (cls == CLS_GAP && d->prev_gap) {
        cls = CLS_NORMAL;
    }
    d->prev_gap = (uint8_t)(cls == CLS_GAP);
    uint8_t good_before = d->good;
    d->good = (cls == CLS_NORMAL) ? (uint8_t)(d->good < GOOD_SAT ? d->good + 1u : GOOD_SAT) : 0u;
    d->last_prev = d->last_tick;
    d->last_tick = tick;
    uint32_t steps = 1u;
    if (cls == CLS_GAP) {
        steps = ECU_TEETH_MISSING + 1u;
    } else if (cls == CLS_DROPOUT) {
        steps = 2u;
    }
    d->period = dt / steps;
    /* averaged sub-tick period for the scheduler: span of the last ECU_PERIOD_HIST intervals */
    d->hist_steps[d->hist_i] = (uint8_t)steps;
    d->hist_tick[d->hist_i] = d->last_prev;     /* start tick of this interval */
    d->hist_i = (uint8_t)((d->hist_i + 1u) % ECU_PERIOD_HIST);
    if (d->hist_n < ECU_PERIOD_HIST) {
        d->hist_n++;
        d->period_q8 = d->period * ECU_Q8;
    } else {
        uint32_t nsteps = 0u;
        for (uint32_t i = 0u; i < ECU_PERIOD_HIST; ++i) {
            nsteps += d->hist_steps[i];
        }
        uint32_t oldest_start = d->hist_tick[d->hist_i];   /* slot after the newest = oldest kept */
        d->period_q8 = (uint32_t)((uint64_t)(tick - oldest_start) * ECU_Q8 / nsteps);
    }

    if (d->sync == ECU_SYNC_NONE) {
        if (cls != CLS_GAP || good_before < ECU_SYNC_MIN_GOOD) {
            return;
        }
        d->sync = ECU_SYNC_CRANK;
        d->tooth = 0u;
        d->rev_seen = 0u;
        d->n_sync++;
    } else if (cls == CLS_GAP) {
        if (d->tooth != ECU_TEETH_REAL - 1) {
            d->n_loss++;              /* gap at the wrong count: resync on this gap */
            d->sync = ECU_SYNC_CRANK;
            d->rev_seen = 0u;
            d->rev_ticks = 0u;
            d->n_sync++;
        } else {
            d->half ^= 1u;
        }
        d->tooth = 0u;
    } else {
        if (cls == CLS_DROPOUT) {
            d->n_dropout++;
        }
        uint32_t t = (uint32_t)d->tooth + steps;
        if (t >= ECU_TEETH_REAL) {    /* expected gap never came */
            d->n_loss++;
            lose_sync(d);
            return;
        }
        d->tooth = (uint8_t)t;
    }
    d->abs_u += steps * (uint32_t)ECU_PITCH_U;

    uint64_t bit = (uint64_t)1u << d->tooth;
    if (d->rev_seen & bit) {
        d->rev_ticks = tick - d->rev_tick[d->tooth];
    }
    d->rev_tick[d->tooth] = tick;
    d->rev_seen |= bit;

    out->valid = 1u;
    out->full = (uint8_t)(d->sync == ECU_SYNC_FULL);
    out->tick = tick;
    out->period_q8 = d->period_q8;
    out->angle_u = tooth_angle(d->tooth, d->half);
    out->span_u = (d->tooth == ECU_TEETH_REAL - 1) ? (ECU_TEETH_MISSING + 1u) * ECU_PITCH_U : (uint32_t)ECU_PITCH_U;
}

void ecu_dec_cam(ecu_decoder_t *d, uint32_t tick)
{
    if (d->sync == ECU_SYNC_NONE || d->period_q8 == 0u) {
        return;
    }
    uint32_t span = (d->tooth == ECU_TEETH_REAL - 1) ? (ECU_TEETH_MISSING + 1u) * ECU_PITCH_U : (uint32_t)ECU_PITCH_U;
    uint64_t frac = (uint64_t)(tick - d->last_tick) * ECU_PITCH_U * ECU_Q8 / d->period_q8;
    if (frac > span) {
        frac = span;
    }
    int32_t best_err = (int32_t)ECU_CYCLE_U;
    uint8_t best_h = 0u;
    int32_t best_a = 0;
    for (uint8_t h = 0u; h < 2u; ++h) {
        int32_t a = ecu_angle_wrap(tooth_angle(d->tooth, h) + (int32_t)frac);
        int32_t e = ecu_angle_wrap((int32_t)ECU_CAM_EDGE_U - a);
        if (e > (int32_t)ECU_CYCLE_U / 2) {
            e = (int32_t)ECU_CYCLE_U - e;
        }
        if (e < best_err) {
            best_err = e;
            best_h = h;
            best_a = a;
        }
    }
    if (d->sync == ECU_SYNC_FULL && best_h != d->half) {
        d->n_cam_err++;
    }
    d->half = best_h;
    d->sync = ECU_SYNC_FULL;
    d->cam_angle_u = best_a;
}
