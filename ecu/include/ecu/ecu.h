/* Redline ECU portable core: public API. C11, deterministic, no heap, no floating point in the
 * crank decoder / angle scheduler / knock controller (those translation units are compiled with
 * -mgeneral-regs-only in the host build to prove it). Floats are used only in the slow loop.
 *
 * Calling contract (any target, see ecu_hal.h):
 *   ecu_init()          once
 *   ecu_crank_edge()    every crank-tooth edge, with the captured timer tick
 *   ecu_cam_edge()      every cam edge
 *   ecu_knock_window()  when the knock window of a cylinder closes (integrated sensor counts)
 *   ecu_slow_step()     at ECU_DT_SLOW (1 kHz), sensor inputs in engineering units
 * Output compare events leave through hal_emit_event(). Angles: ADR 0001 (0 = cyl 1 firing TDC).
 */
#ifndef ECU_H
#define ECU_H

#include <stdint.h>

#include "ecu/ecu_calib_gen.h"
#include "ecu/ecu_types.h"

#define ECU_Q8 256u           /* fixed-point scale of period_q8 */

typedef struct {
    uint8_t sync;            /* ecu_sync_t */
    uint8_t have_last;
    uint8_t have_period;
    uint8_t tooth;           /* slot of the last accepted tooth, 0..ECU_TEETH_REAL-1 (0 = after gap) */
    uint8_t half;            /* revolution half of `tooth` once ECU_SYNC_FULL */
    uint8_t good;            /* consecutive normal teeth (gap plausibility while searching) */
    uint8_t prev_gap;        /* last accepted edge was classified as the gap */
    uint32_t last_tick;
    uint32_t last_prev;      /* tick of the accepted edge before last_tick */
    uint32_t period;         /* ticks per tooth pitch from the last interval (classification) */
    uint32_t period_q8;      /* ticks per pitch x 256, averaged over the last ECU_PERIOD_HIST intervals */
    uint32_t hist_tick[ECU_PERIOD_HIST];
    uint8_t hist_steps[ECU_PERIOD_HIST];
    uint8_t hist_n, hist_i;
    uint32_t rev_ticks;      /* ticks for the last full revolution (same tooth), 0 = unknown */
    uint32_t rev_tick[ECU_TEETH_REAL];
    uint64_t rev_seen;       /* bit per tooth slot: rev_tick[] valid */
    uint32_t abs_u;          /* monotonic crank angle [angle units], wraps */
    int32_t cam_angle_u;     /* measured cycle angle of the last cam edge */
    uint32_t n_edges, n_sync, n_loss, n_noise, n_dropout, n_cam_err;
} ecu_decoder_t;

typedef struct {
    int32_t adv_u[ECU_NCYL];         /* spark advance BTDC [angle units] (slow loop writes) */
    uint32_t pw_ticks[ECU_NCYL];     /* injector electrical pulse [ticks], 0 = no injection */
    uint32_t last_fire_abs[ECU_NCYL][3];
    uint8_t fired_once[ECU_NCYL][3];
    uint8_t rev_cut;                 /* rev limiter: fuel (and spark) cut active */
    uint8_t limp;                    /* use the limp-mode rev limit */
    uint32_t inj_count[ECU_NCYL];    /* injections committed (slow loop film model reads) */
    uint32_t inj_pw[ECU_NCYL];
    uint32_t n_cut;
} ecu_sched_t;

typedef struct {
    uint32_t now_tick;       /* timer now (stall detection) */
    float map_pa;
    float iat_k;
    float clt_k;
    float tps1_deg;          /* throttle opening from the closed stop, sensor 1 */
    float tps2_deg;          /* complementary sensor: max opening - opening */
    float pedal;             /* 0..1 */
    float lambda;            /* wideband reading */
    float vbatt;
} ecu_inputs_t;

typedef struct {
    float etc_voltage;       /* throttle motor command [V] */
    float rpm;
    float map_used;
    float throttle_target_deg;
    float spark_deg[ECU_NCYL];
    float knock_retard_deg[ECU_NCYL];
    float pw_s[ECU_NCYL];
    float film_kg[ECU_NCYL];
    float stft;
    float ltft;
    float idle_air_deg;
    float closed_loop;
    float idle_active;
    float limp;
    float faults;
    float sync;
    float rev_cut;
} ecu_outputs_t;

typedef struct {
    ecu_decoder_t dec;
    ecu_sched_t sch;
    int32_t knock_retard_u[ECU_NCYL];
    uint32_t n_knock;
    /* slow-loop state (float) */
    float stft_i, stft, ltft, cl_timer;
    float map_prev, dmap, freeze_timer;
    uint8_t map_init;
    float film[ECU_NCYL];
    uint32_t inj_seen[ECU_NCYL];
    float etc_i, etc_prev_e;
    float idle_i, idle_prev_e;
    uint16_t fault_cnt[6];
    uint16_t faults;
    uint8_t limp;
    uint8_t closed_loop;
    uint8_t idle_active;
} ecu_t;

void ecu_init(ecu_t *e);
void ecu_crank_edge(ecu_t *e, uint32_t tick);
void ecu_cam_edge(ecu_t *e, uint32_t tick);
void ecu_knock_window(ecu_t *e, uint8_t cyl, uint32_t counts);
void ecu_slow_step(ecu_t *e, const ecu_inputs_t *in, ecu_outputs_t *out);

/* fixed-point helpers (exposed for tests) */
typedef struct {
    uint8_t valid;           /* accepted tooth with known crank angle */
    uint8_t full;            /* cycle angle valid (cam phase known) */
    uint32_t tick;
    uint32_t period_q8;      /* ticks per pitch x ECU_Q8 */
    int32_t angle_u;         /* cycle angle of this tooth (full) [angle units] */
    uint32_t span_u;         /* angle to the next real tooth */
} ecu_tooth_t;

void ecu_dec_init(ecu_decoder_t *d);
void ecu_dec_crank(ecu_decoder_t *d, uint32_t tick, ecu_tooth_t *out);
void ecu_dec_cam(ecu_decoder_t *d, uint32_t tick);
void ecu_dec_stall(ecu_decoder_t *d);
void ecu_sched_init(ecu_sched_t *s);
void ecu_sched_tooth(ecu_sched_t *s, const ecu_decoder_t *d, const ecu_tooth_t *t);
int32_t ecu_angle_wrap(int32_t a);   /* into [0, ECU_CYCLE_U) */

#endif
