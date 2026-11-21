/* Slow loop (ECU_DT_SLOW): faults and limp mode, speed-density fueling with inverse x-tau wall
 * wetting and injector dead time, closed-loop lambda PI with STFT/LTFT, spark map with knock
 * retard and idle spark, pedal map, idle speed control, electronic throttle PID with spring and
 * friction feed-forward. Float arithmetic is allowed here (not in the angle scheduler). */
#include <math.h>

#include "ecu/ecu.h"
#include "tables.h"

#define SEC_PER_MIN 60.0f
#define STROKES_PER_CYCLE_REV 2.0f     /* crank revolutions per four-stroke cycle */
#define RPM_FLOOR 100.0f               /* film time-constant evaluation floor */
#define DEG_PER_REV 360.0f

enum { F_MAP = 0, F_TPS, F_LAMBDA, F_IAT, F_CLT, F_PEDAL, F_N };

void ecu_init(ecu_t *e)
{
    uint8_t *p = (uint8_t *)e;
    for (uint32_t i = 0; i < sizeof(*e); ++i) {
        p[i] = 0u;
    }
    ecu_dec_init(&e->dec);
    ecu_sched_init(&e->sch);
}

static float clampf(float x, float lo, float hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

static float sgnf(float x)
{
    return x > 0.0f ? 1.0f : (x < 0.0f ? -1.0f : 0.0f);
}

static void fault_check(ecu_t *e, int idx, int bad)
{
    if (bad) {
        if (e->fault_cnt[idx] < ECU_FAULT_DEBOUNCE) {
            e->fault_cnt[idx]++;
        }
        if (e->fault_cnt[idx] >= ECU_FAULT_DEBOUNCE) {
            e->faults |= (uint16_t)(1u << idx);      /* latched until ecu_init */
        }
    } else if (e->fault_cnt[idx] > 0u && !(e->faults & (1u << idx))) {
        e->fault_cnt[idx]--;
    }
}

static float rpm_of(const ecu_decoder_t *d)
{
    if (d->sync == ECU_SYNC_NONE || !d->have_period) {
        return 0.0f;
    }
    if (d->rev_ticks) {
        return SEC_PER_MIN * ECU_TIMER_HZ_F / (float)d->rev_ticks;
    }
    return SEC_PER_MIN * ECU_TIMER_HZ_F / ((float)d->period * (float)ECU_TEETH_TOTAL);
}

static int32_t to_units(float deg)
{
    float u = deg * (float)ECU_ANGLE_U_PER_DEG;
    return (int32_t)(u >= 0.0f ? u + 0.5f : u - 0.5f);
}

void ecu_slow_step(ecu_t *e, const ecu_inputs_t *in, ecu_outputs_t *out)
{
    const float dt = ECU_DT_SLOW;
    ecu_decoder_t *d = &e->dec;
    ecu_sched_t *s = &e->sch;

    if (d->have_last && (uint32_t)(in->now_tick - d->last_tick) > ECU_STALL_TICKS) {
        ecu_dec_stall(d);
    }
    const float rpm = rpm_of(d);
    const int full = d->sync == ECU_SYNC_FULL;

    /* ---- sensor faults and limp mode ---- */
    fault_check(e, F_MAP, in->map_pa < ECU_MAP_MIN || in->map_pa > ECU_MAP_MAX);
    fault_check(e, F_TPS, fabsf(in->tps1_deg + in->tps2_deg - ECU_THR_MAX_DEG) > ECU_TPS_SUM_TOL ||
                              in->tps1_deg < -ECU_TPS_SUM_TOL || in->tps1_deg > ECU_THR_MAX_DEG + ECU_TPS_SUM_TOL);
    fault_check(e, F_LAMBDA, in->lambda < ECU_LAMBDA_MIN || in->lambda > ECU_LAMBDA_MAX);
    fault_check(e, F_IAT, in->iat_k < ECU_TEMP_MIN || in->iat_k > ECU_TEMP_MAX);
    fault_check(e, F_CLT, in->clt_k < ECU_TEMP_MIN || in->clt_k > ECU_TEMP_MAX);
    fault_check(e, F_PEDAL, in->pedal < -ECU_IDLE_PEDAL_MAX || in->pedal > 1.0f + ECU_IDLE_PEDAL_MAX);
    e->limp = (uint8_t)((e->faults & (ECU_FAULT_TPS | ECU_FAULT_PEDAL)) != 0u);
    s->limp = e->limp;
    const float thr_pos = (e->faults & ECU_FAULT_TPS) ? ECU_ETC_LIMP_DEG : in->tps1_deg;
    const float map = (e->faults & ECU_FAULT_MAP)
                          ? ecu_interp1(ECU_LIMP_OPEN_X, ECU_LIMP_MAP_PA, ECU_LIMP_OPEN_X_N, thr_pos)
                          : in->map_pa;
    const float iat = (e->faults & ECU_FAULT_IAT) ? ECU_IAT_DEFAULT : in->iat_k;
    const float pedal = (e->faults & ECU_FAULT_PEDAL) ? 0.0f : clampf(in->pedal, 0.0f, 1.0f);

    /* ---- MAP rate (filtered) for the transient fuel prediction and the lambda freeze ---- */
    if (!e->map_init) {
        e->map_prev = map;
        e->map_init = 1u;
    }
    e->dmap += ((map - e->map_prev) / dt - e->dmap) * dt / (ECU_MAP_DERIV_TAU + dt);
    e->map_prev = map;
    e->freeze_timer = (fabsf(e->dmap) > ECU_CL_FREEZE_DMAP) ? ECU_CL_FREEZE_HOLD
                                                            : (e->freeze_timer > dt ? e->freeze_timer - dt : 0.0f);

    /* ---- closed-loop lambda: PI -> STFT, slow transfer into LTFT ---- */
    int cl_ok = full && !(e->faults & ECU_FAULT_LAMBDA) && rpm > ECU_CL_MIN_RPM && !s->rev_cut;
    e->cl_timer = cl_ok ? e->cl_timer + dt : 0.0f;
    e->closed_loop = (uint8_t)(e->cl_timer >= ECU_CL_DELAY);
    if (e->closed_loop && e->freeze_timer <= 0.0f) {
        float err = in->lambda / ECU_LAMBDA_TARGET - 1.0f;        /* lean -> positive -> add fuel */
        e->stft_i = clampf(e->stft_i + ECU_LAM_KI * err * dt, -ECU_STFT_LIM, ECU_STFT_LIM);
        float xfer = ECU_LTFT_RATE * e->stft_i * dt;
        e->ltft = clampf(e->ltft + xfer, -ECU_LTFT_LIM, ECU_LTFT_LIM);
        e->stft_i -= xfer;
        e->stft = clampf(ECU_LAM_KP * err + e->stft_i, -ECU_STFT_LIM, ECU_STFT_LIM);
    } else if (e->faults & ECU_FAULT_LAMBDA) {
        e->stft = 0.0f;
        e->stft_i = 0.0f;
    }

    /* ---- speed-density fueling, inverse x-tau, dead time ----
     * Event form of the Aquino film inverse. Port fuel inducted at IVC k is
     *   M(IVC k-1) - M(IVC k) + m_k,   M(IVC k) = (M(IVC k-1) e1 + x m_k) e2,
     * with e1 the film decay from IVC k-1 to the injection, e2 from the injection (EOI) to IVC k and
     * e = e1 e2 the decay over one cycle, so
     *   m_k = (m_des - M(IVC k-1) (1 - e)) / (1 - x e2),   M(IVC k-1) = M_pre e2 / e
     * gives exactly m_des in the cylinder when the film model is right. M_pre is the film just
     * before this cylinder's injection: the pulse is recomputed every slow tick, and the value in
     * force when the scheduler commits the injection is the one used. */
    const float rpm_f = rpm > RPM_FLOOR ? rpm : RPM_FLOOR;
    const float tau = ECU_FILM_TAU_REF * powf(rpm_f / ECU_FILM_RPM_REF, ECU_FILM_TAU_EXP);
    const float t_cycle = SEC_PER_MIN * STROKES_PER_CYCLE_REV / rpm_f;
    const float t_eoi_ivc = ECU_EOI_TO_IVC_DEG / DEG_PER_REV * SEC_PER_MIN / rpm_f;
    const float e_cycle = expf(-t_cycle / tau);
    const float e2 = expf(-t_eoi_ivc / tau);
    const float dead = ecu_interp1(ECU_DEAD_VBATT, ECU_DEAD_S, ECU_DEAD_VBATT_N, in->vbatt);
    const float pw_prev = (float)s->pw_ticks[0] / ECU_TIMER_HZ_F;
    float map_ivc = map + ECU_MAP_PRED_GAIN * e->dmap * (pw_prev + t_eoi_ivc);   /* MAP expected at IVC */
    map_ivc = clampf(map_ivc, ECU_MAP_MIN, ECU_MAP_MAX);
    const float ve = ecu_interp2(ECU_TAB_RPM, ECU_TAB_RPM_N, ECU_TAB_MAP, ECU_TAB_MAP_N, &ECU_VE_TABLE[0][0], rpm, map_ivc);
    const float air = ve * map_ivc * ECU_VD_CYL / (ECU_R_AIR * iat);
    const float m_des = air / (ECU_AFR_STOICH * ECU_LAMBDA_TARGET) * (1.0f + e->stft + e->ltft);
    for (int c = 0; c < ECU_NCYL; ++c) {
        if (s->inj_count[c] != e->inj_seen[c]) {                 /* film gets x of what was injected */
            e->inj_seen[c] = s->inj_count[c];
            float m = ((float)s->inj_pw[c] / ECU_TIMER_HZ_F - dead) * ECU_INJ_FLOW;
            e->film[c] += ECU_FILM_X * (m > 0.0f ? m : 0.0f);
        }
        e->film[c] -= e->film[c] * dt / tau;
        float m_inj = (m_des - e->film[c] * e2 / e_cycle * (1.0f - e_cycle)) / (1.0f - ECU_FILM_X * e2);
        float pw = (m_inj > 0.0f && full && rpm > 0.0f) ? m_inj / ECU_INJ_FLOW + dead : 0.0f;
        s->pw_ticks[c] = (uint32_t)(pw * ECU_TIMER_HZ_F + 0.5f);
        out->pw_s[c] = pw;
        out->film_kg[c] = e->film[c];
    }

    /* ---- idle speed control (air PID + spark P) ---- */
    e->idle_active = (uint8_t)(full && !e->limp && pedal < ECU_IDLE_PEDAL_MAX && rpm < ECU_IDLE_RPM + ECU_IDLE_RPM_BAND);
    const float idle_err = ECU_IDLE_RPM - rpm;
    float idle_spark = 0.0f;
    float target;
    if (e->idle_active) {
        e->idle_i = clampf(e->idle_i + ECU_IDLE_KI * idle_err * dt, -ECU_IDLE_TRIM_LIM, ECU_IDLE_TRIM_LIM);
        float trim = clampf(ECU_IDLE_KP * idle_err + e->idle_i + ECU_IDLE_KD * (idle_err - e->idle_prev_e) / dt,
                            -ECU_IDLE_TRIM_LIM, ECU_IDLE_TRIM_LIM);
        target = ECU_IDLE_BASE_DEG + trim;
        idle_spark = clampf(ECU_IDLE_SPARK_KP * idle_err, -ECU_IDLE_SPARK_LIM, ECU_IDLE_SPARK_LIM) - ECU_IDLE_SPARK_RESERVE;
    } else {
        target = ecu_interp1(ECU_PEDAL_X, ECU_PEDAL_DEG, ECU_PEDAL_X_N, pedal);
    }
    e->idle_prev_e = idle_err;
    target = clampf(target, 0.0f, ECU_THR_MAX_DEG);

    /* ---- spark: table - knock retard (+ idle spark) ---- */
    const float base = ecu_interp2(ECU_TAB_RPM, ECU_TAB_RPM_N, ECU_TAB_MAP, ECU_TAB_MAP_N, &ECU_SPARK_TABLE[0][0], rpm, map);
    for (int c = 0; c < ECU_NCYL; ++c) {
        float kr = (float)e->knock_retard_u[c] / (float)ECU_ANGLE_U_PER_DEG;
        float adv = clampf(base + idle_spark - kr, ECU_ADV_MIN, ECU_ADV_MAX);
        s->adv_u[c] = to_units(adv);
        out->spark_deg[c] = adv;
        out->knock_retard_deg[c] = kr;
    }

    /* ---- electronic throttle: PID + spring/friction feed-forward; unpowered in limp ---- */
    float v = 0.0f;
    if (!e->limp) {
        float err = target - thr_pos;
        e->etc_i = clampf(e->etc_i + ECU_ETC_KI * err * dt, -ECU_ETC_ILIM, ECU_ETC_ILIM);
        float der = (err - e->etc_prev_e) / dt;
        e->etc_prev_e = err;
        float off = target - ECU_ETC_LIMP_DEG;
        float torque = ECU_ETC_SPRING_RATE * off / ECU_DEG_PER_RAD + ECU_ETC_PRELOAD * sgnf(off);
        if (fabsf(err) > ECU_ETC_FF_DEADBAND) {
            torque += ECU_ETC_COULOMB * sgnf(err);
        }
        float ff = ECU_ETC_FF_SCALE * ECU_ETC_R_OVER_KT * torque;
        float vmax = in->vbatt < ECU_ETC_VMAX ? in->vbatt : ECU_ETC_VMAX;
        v = clampf(ECU_ETC_KP * err + e->etc_i + ECU_ETC_KD * der + ff, -vmax, vmax);
    } else {
        e->etc_i = 0.0f;
        target = ECU_ETC_LIMP_DEG;
    }

    out->etc_voltage = v;
    out->rpm = rpm;
    out->map_used = map;
    out->throttle_target_deg = target;
    out->stft = e->stft;
    out->ltft = e->ltft;
    out->idle_air_deg = e->idle_i;
    out->closed_loop = (float)e->closed_loop;
    out->idle_active = (float)e->idle_active;
    out->limp = (float)e->limp;
    out->faults = (float)e->faults;
    out->sync = (float)d->sync;
    out->rev_cut = (float)s->rev_cut;
}
