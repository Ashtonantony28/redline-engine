/* Redline ECU C unit tests: decoder, cam phase, scheduler accuracy, knock, rev limiter, slow loop.
 * Synthetic tooth streams come from an exact kinematic wheel model (double precision, test only). */
#include <stdint.h>
#include <string.h>

#include "ecu/ecu.h"
#include "tables.h"
#include "tharness.h"

int hil_drain(uint32_t *tick, uint32_t *arg, uint8_t *type, uint8_t *cyl, int max);
void hil_reset(void);

#define DEG_PER_REV 360.0
#define DEG_PER_CYCLE 720.0
#define US 1.0e6
#define MAXEV 4096

/* ---------------- synthetic wheel ---------------- */
typedef struct {
    double theta0;      /* cycle angle at t = 0 [deg], unwrapped */
    double omega;       /* deg/s at t = 0 */
    double alpha;       /* deg/s^2 */
    uint32_t tick0;
    int noise_every;    /* insert a noise edge at noise_frac pitch after every n-th tooth (0 = off) */
    double noise_frac;
    long drop_k;        /* drop the tooth with this absolute index (-1 = none) */
    long extra_k;       /* insert an extra (accepted) edge half a pitch after this tooth (-1 = none) */
    int cam;            /* deliver cam edges */
} wheel_t;

static double pitch_deg(void) { return DEG_PER_REV / ECU_TEETH_TOTAL; }
static double tooth0_deg(void) { return (double)ECU_TOOTH0_U / ECU_ANGLE_U_PER_DEG; }
static double cam_deg(void) { return (double)ECU_CAM_EDGE_U / ECU_ANGLE_U_PER_DEG; }

static double time_of(const wheel_t *w, double theta)
{
    double d = theta - w->theta0;
    if (w->alpha == 0.0) {
        return d / w->omega;
    }
    return (-w->omega + sqrt(w->omega * w->omega + 2.0 * w->alpha * d)) / w->alpha;
}

static double theta_of(const wheel_t *w, double t) { return w->theta0 + w->omega * t + 0.5 * w->alpha * t * t; }

static uint32_t tick_of(const wheel_t *w, double t) { return w->tick0 + (uint32_t)llround(t * US); }

static double t_of_tick(const wheel_t *w, uint32_t tick) { return (double)(uint32_t)(tick - w->tick0) / US; }

typedef enum { E_TOOTH, E_NOISE, E_CAM } ekind_t;
typedef struct { double theta; ekind_t kind; long k; } edge_t;

/* Build the sorted edge list in (theta_start, theta_end]. */
static int build_edges(const wheel_t *w, double th_s, double th_e, edge_t *out, int max)
{
    int n = 0;
    double p = pitch_deg();
    long k0 = (long)floor((th_s - tooth0_deg()) / p) + 1;
    for (long k = k0;; ++k) {
        double th = tooth0_deg() + (double)k * p;
        if (th > th_e) break;
        long slot = ((k % ECU_TEETH_TOTAL) + ECU_TEETH_TOTAL) % ECU_TEETH_TOTAL;
        if (slot >= ECU_TEETH_REAL || k == w->drop_k) continue;
        if (n < max) out[n++] = (edge_t){th, E_TOOTH, k};
        if (k == w->extra_k && n < max) out[n++] = (edge_t){th + 0.5 * p, E_NOISE, k};
        if (w->noise_every && (k % w->noise_every) == 0 && n < max) out[n++] = (edge_t){th + w->noise_frac * p, E_NOISE, k};
    }
    if (w->cam) {
        long j0 = (long)floor((th_s - cam_deg()) / DEG_PER_CYCLE) + 1;
        for (long j = j0;; ++j) {
            double th = cam_deg() + (double)j * DEG_PER_CYCLE;
            if (th > th_e) break;
            if (n < max) out[n++] = (edge_t){th, E_CAM, j};
        }
    }
    /* insertion sort by angle (small lists) */
    for (int i = 1; i < n; ++i) {
        edge_t x = out[i];
        int j = i - 1;
        while (j >= 0 && out[j].theta > x.theta) { out[j + 1] = out[j]; --j; }
        out[j + 1] = x;
    }
    return n;
}

static edge_t g_edges[20000];

typedef struct {
    int teeth_checked, bad_angle, synced_at_edge;
    double first_sync_theta;
} dec_stats_t;

/* Drive the decoder only; verify the reported angle of every accepted tooth once synced. */
static dec_stats_t drive_decoder(ecu_decoder_t *d, const wheel_t *w, double th_s, double th_e)
{
    dec_stats_t st = {0, 0, -1, -1.0};
    int n = build_edges(w, th_s, th_e, g_edges, (int)(sizeof g_edges / sizeof g_edges[0]));
    for (int i = 0; i < n; ++i) {
        uint32_t tk = tick_of(w, time_of(w, g_edges[i].theta));
        if (g_edges[i].kind == E_CAM) {
            ecu_dec_cam(d, tk);
            continue;
        }
        ecu_tooth_t t;
        ecu_dec_crank(d, tk, &t);
        if (d->sync != ECU_SYNC_NONE && st.synced_at_edge < 0) {
            st.synced_at_edge = i;
            st.first_sync_theta = g_edges[i].theta;
        }
        if (t.valid && g_edges[i].kind == E_TOOTH) {
            double truth = fmod(g_edges[i].theta, DEG_PER_CYCLE) * ECU_ANGLE_U_PER_DEG;
            int32_t rev_u = (int32_t)((ECU_TOOTH0_U + (int32_t)d->tooth * ECU_PITCH_U) % ECU_REV_U);
            int32_t truth_rev = (int32_t)llround(fmod(truth, (double)ECU_REV_U));
            st.teeth_checked++;
            if (rev_u != truth_rev) st.bad_angle++;
            if (t.full && t.angle_u != (int32_t)llround(truth)) st.bad_angle++;
        }
    }
    return st;
}

static wheel_t wheel(double rpm, double theta0)
{
    wheel_t w;
    memset(&w, 0, sizeof w);
    w.theta0 = theta0;
    w.omega = rpm * DEG_PER_REV / 60.0;
    w.drop_k = -1;
    w.extra_k = -1;
    w.cam = 1;
    return w;
}

/* ---------------- decoder tests ---------------- */
static void test_sync_within_one_rev_after_gap(void)
{
    for (double start = 0.0; start < DEG_PER_CYCLE; start += 37.0) {
        ecu_decoder_t d;
        ecu_dec_init(&d);
        wheel_t w = wheel(3000.0, start);
        dec_stats_t st = drive_decoder(&d, &w, start, start + 3.0 * DEG_PER_CYCLE);
        /* the gap ends at tooth 0 (tooth0_deg mod 360); sync must come at the first tooth-0 edge */
        double first_gap_end = tooth0_deg();
        while (first_gap_end <= start + 5.0 * pitch_deg()) first_gap_end += DEG_PER_REV;   /* >= 2 normal teeth first */
        CHECK(st.synced_at_edge >= 0);
        CHECK(st.first_sync_theta <= first_gap_end + 1e-9);
        CHECK(st.first_sync_theta - start <= DEG_PER_REV + 6.0 * pitch_deg() + 1e-9);
        CHECK(st.bad_angle == 0);
        CHECK(st.teeth_checked > 100);
        CHECK(d.sync == ECU_SYNC_FULL);
        CHECK(d.n_loss == 0u);
    }
}

static void test_cam_phase_and_rpm(void)
{
    ecu_decoder_t d;
    ecu_dec_init(&d);
    wheel_t w = wheel(2400.0, 10.0);
    dec_stats_t st = drive_decoder(&d, &w, 10.0, 10.0 + 4.0 * DEG_PER_CYCLE);
    CHECK(st.bad_angle == 0);
    CHECK(d.sync == ECU_SYNC_FULL);
    CHECK(d.n_cam_err == 0u);
    CHECK_NEAR(d.cam_angle_u, ECU_CAM_EDGE_U, ECU_ANGLE_U_PER_DEG / 4);
    CHECK_NEAR(60.0 * ECU_TIMER_HZ / (double)d.rev_ticks, 2400.0, 1.0);
}

static void test_noise_rejection(void)
{
    ecu_decoder_t d;
    ecu_dec_init(&d);
    wheel_t w = wheel(4000.0, 100.0);
    w.noise_every = 7;
    w.noise_frac = 0.2;
    dec_stats_t st = drive_decoder(&d, &w, 100.0, 100.0 + 4.0 * DEG_PER_CYCLE);
    CHECK(d.n_noise > 50u);
    CHECK(d.n_loss == 0u);
    CHECK(st.bad_angle == 0);
    CHECK(d.sync == ECU_SYNC_FULL);
}

static void test_single_tooth_dropout(void)
{
    ecu_decoder_t d;
    ecu_dec_init(&d);
    wheel_t w = wheel(3000.0, 0.0);
    w.drop_k = 150;     /* mid-wheel tooth in the 3rd revolution */
    dec_stats_t st = drive_decoder(&d, &w, 0.0, 4.0 * DEG_PER_CYCLE);
    CHECK(d.n_dropout == 1u);
    CHECK(d.n_loss == 0u);
    CHECK(st.bad_angle == 0);
}

static void test_loss_of_sync_and_recovery(void)
{
    ecu_decoder_t d;
    ecu_dec_init(&d);
    wheel_t w = wheel(3000.0, 0.0);
    w.extra_k = 200;    /* accepted spurious edge: count is off by one at the next gap */
    dec_stats_t st = drive_decoder(&d, &w, 0.0, 6.0 * DEG_PER_CYCLE);
    CHECK(d.n_loss == 1u);
    CHECK(d.n_sync == 2u);
    CHECK(d.sync == ECU_SYNC_FULL);
    /* teeth between the spurious edge and the resync report a shifted angle; after resync all good */
    CHECK(st.bad_angle > 0);
    CHECK(st.bad_angle < ECU_TEETH_REAL * 2);
    ecu_decoder_t d2 = d;
    wheel_t w2 = wheel(3000.0, 6.0 * DEG_PER_CYCLE);
    w2.tick0 = tick_of(&w, time_of(&w, 6.0 * DEG_PER_CYCLE));
    dec_stats_t st2 = drive_decoder(&d2, &w2, 6.0 * DEG_PER_CYCLE, 8.0 * DEG_PER_CYCLE);
    CHECK(st2.bad_angle == 0);
    CHECK(d2.n_loss == 1u);
}

static void test_stall_and_timer_wrap(void)
{
    ecu_decoder_t d;
    ecu_dec_init(&d);
    wheel_t w = wheel(1500.0, 0.0);
    w.tick0 = 0xFFFF0000u;   /* timer wraps during the stream */
    dec_stats_t st = drive_decoder(&d, &w, 0.0, 6.0 * DEG_PER_CYCLE);
    CHECK(st.bad_angle == 0);
    CHECK(d.n_loss == 0u);
    CHECK_NEAR(60.0 * ECU_TIMER_HZ / (double)d.rev_ticks, 1500.0, 1.0);
    ecu_dec_stall(&d);
    CHECK(d.sync == ECU_SYNC_NONE);
    CHECK(d.rev_ticks == 0u);
}

/* ---------------- scheduler ---------------- */
static uint32_t ev_tick[MAXEV], ev_arg[MAXEV];
static uint8_t ev_type[MAXEV], ev_cyl[MAXEV];

static int drive_ecu(ecu_t *e, const wheel_t *w, double th_s, double th_e)
{
    int n = build_edges(w, th_s, th_e, g_edges, (int)(sizeof g_edges / sizeof g_edges[0]));
    int nev = 0;
    for (int i = 0; i < n; ++i) {
        uint32_t tk = tick_of(w, time_of(w, g_edges[i].theta));
        if (g_edges[i].kind == E_CAM) ecu_cam_edge(e, tk);
        else ecu_crank_edge(e, tk);
        nev += hil_drain(ev_tick + nev, ev_arg + nev, ev_type + nev, ev_cyl + nev, MAXEV - nev);
    }
    return nev;
}

static double wrap180(double a)
{
    a = fmod(a, DEG_PER_CYCLE);
    if (a > DEG_PER_CYCLE / 2) a -= DEG_PER_CYCLE;
    if (a <= -DEG_PER_CYCLE / 2) a += DEG_PER_CYCLE;
    return a;
}

static void check_schedule(double rpm, double alpha_rpm_s, double cycles, double *max_err_out)
{
    hil_reset();
    ecu_t e;
    ecu_init(&e);
    wheel_t w = wheel(rpm, 0.0);
    w.alpha = alpha_rpm_s * DEG_PER_REV / 60.0;
    const double adv = 30.0;
    const uint32_t pw = 4000u;
    for (int c = 0; c < ECU_NCYL; ++c) {
        e.sch.adv_u[c] = (int32_t)(adv * ECU_ANGLE_U_PER_DEG);
        e.sch.pw_ticks[c] = pw;
    }
    int nev = drive_ecu(&e, &w, 0.0, cycles * DEG_PER_CYCLE);
    int nspark = 0, ninj = 0, ndwell = 0;
    double max_err = 0.0, max_dwell_err = 0.0;
    double last_spark_th[ECU_NCYL];
    int bad_spacing = 0;
    for (int c = 0; c < ECU_NCYL; ++c) last_spark_th[c] = -1.0;
    for (int i = 0; i < nev; ++i) {
        double t = t_of_tick(&w, ev_tick[i]);
        double th = theta_of(&w, t);
        double tdc = (double)ECU_TDC_U[ev_cyl[i]] / ECU_ANGLE_U_PER_DEG;
        double om = w.omega + w.alpha * t;
        if (ev_type[i] == ECU_EV_SPARK) {
            double err = fabs(wrap180(th - (tdc - adv)));
            if (err > max_err) max_err = err;
            ++nspark;
            /* exactly one spark per cylinder per cycle: no misses, no doubles */
            if (last_spark_th[ev_cyl[i]] >= 0.0 && fabs(th - last_spark_th[ev_cyl[i]] - DEG_PER_CYCLE) > 1.0) ++bad_spacing;
            last_spark_th[ev_cyl[i]] = th;
            /* matching dwell start: 3 ms earlier (interpolated) */
            for (int j = i - 1; j >= 0; --j) {
                if (ev_type[j] == ECU_EV_DWELL && ev_cyl[j] == ev_cyl[i]) {
                    double dw = (double)(uint32_t)(ev_tick[i] - ev_tick[j]) / US;
                    double derr = fabs(dw / ((double)ECU_DWELL_TICKS / US) - 1.0);   /* relative */
                    if (derr > max_dwell_err) max_dwell_err = derr;
                    break;
                }
            }
        } else if (ev_type[i] == ECU_EV_INJ) {
            /* injection start = EOI - pulse converted to angle at the current speed */
            double target = tdc + (double)ECU_EOI_U / ECU_ANGLE_U_PER_DEG - (double)ev_arg[i] / US * om;
            double err = fabs(wrap180(th - target));
            if (err > max_err) max_err = err;
            ++ninj;
        } else if (ev_type[i] == ECU_EV_DWELL) {
            ++ndwell;
        }
    }
    CHECK(nspark >= 3 * ECU_NCYL);
    CHECK(bad_spacing == 0);
    CHECK(ninj >= 3 * ECU_NCYL);
    CHECK(ndwell >= 3 * ECU_NCYL);
    printf("    dwell duration error %.3f %%\n", 100.0 * max_dwell_err);
    CHECK(max_dwell_err < 0.01);
    *max_err_out = max_err;
}

static void test_scheduler_accuracy_7000rpm(void)
{
    double err;
    check_schedule(7000.0, 0.0, 6.0, &err);
    printf("    7000 rpm steady: max event angle error %.4f deg\n", err);
    CHECK(err < 0.5);
    check_schedule(5800.0, 12000.0, 4.0, &err);   /* WOT-like free-rev acceleration toward redline */
    printf("    5800 rpm + 12000 rpm/s: max event angle error %.4f deg\n", err);
    CHECK(err < 0.5);
}

static void test_rev_limiter(void)
{
    hil_reset();
    ecu_t e;
    ecu_init(&e);
    for (int c = 0; c < ECU_NCYL; ++c) {
        e.sch.adv_u[c] = 20 * ECU_ANGLE_U_PER_DEG;
        e.sch.pw_ticks[c] = 3000u;
    }
    wheel_t w = wheel(7100.0, 0.0);
    drive_ecu(&e, &w, 0.0, 2.0 * DEG_PER_CYCLE);
    CHECK(e.sch.rev_cut == 1u);
    int nev = drive_ecu(&e, &w, 2.0 * DEG_PER_CYCLE, 4.0 * DEG_PER_CYCLE);
    CHECK(nev == 0);
    CHECK(e.sch.n_cut > 0u);
    wheel_t w2 = wheel(6700.0, 4.0 * DEG_PER_CYCLE);
    w2.tick0 = tick_of(&w, time_of(&w, 4.0 * DEG_PER_CYCLE));
    nev = drive_ecu(&e, &w2, 4.0 * DEG_PER_CYCLE, 7.0 * DEG_PER_CYCLE);
    CHECK(e.sch.rev_cut == 0u);
    CHECK(nev > 12);
}

/* ---------------- knock ---------------- */
static void test_knock_retard_and_recover(void)
{
    ecu_t e;
    ecu_init(&e);
    for (int i = 0; i < 3; ++i) ecu_knock_window(&e, 1, ECU_KNOCK_THRESHOLD + 1u);
    CHECK(e.knock_retard_u[1] == 3 * ECU_KNOCK_STEP_U);
    CHECK(e.knock_retard_u[0] == 0);
    for (int i = 0; i < 10; ++i) ecu_knock_window(&e, 1, ECU_KNOCK_THRESHOLD + 1u);
    CHECK(e.knock_retard_u[1] == ECU_KNOCK_MAX_U);
    int clean = 0;
    while (e.knock_retard_u[1] > 0 && clean < 10000) { ecu_knock_window(&e, 1, 0u); ++clean; }
    CHECK(clean == (ECU_KNOCK_MAX_U + ECU_KNOCK_RECOVER_U - 1) / ECU_KNOCK_RECOVER_U);
}

/* ---------------- slow loop ---------------- */
static ecu_inputs_t nominal_inputs(uint32_t now)
{
    ecu_inputs_t in;
    in.now_tick = now;
    in.map_pa = 0.5e5f;
    in.iat_k = 315.0f;
    in.clt_k = 360.0f;
    in.tps1_deg = 5.0f;
    in.tps2_deg = ECU_THR_MAX_DEG - 5.0f;
    in.pedal = 0.3f;
    in.lambda = 1.0f;
    in.vbatt = 13.8f;
    return in;
}

static void test_speed_density_pulse_width(void)
{
    hil_reset();
    ecu_t e;
    ecu_init(&e);
    wheel_t w = wheel(2000.0, 0.0);
    drive_ecu(&e, &w, 0.0, 2.0 * DEG_PER_CYCLE);
    CHECK(e.dec.sync == ECU_SYNC_FULL);
    uint32_t now = tick_of(&w, time_of(&w, 2.0 * DEG_PER_CYCLE));
    ecu_inputs_t in = nominal_inputs(now);
    ecu_outputs_t out;
    ecu_slow_step(&e, &in, &out);
    CHECK_NEAR(out.rpm, 2000.0, 1.0);
    float ve = ecu_interp2(ECU_TAB_RPM, ECU_TAB_RPM_N, ECU_TAB_MAP, ECU_TAB_MAP_N, &ECU_VE_TABLE[0][0], out.rpm, in.map_pa);
    double air = ve * in.map_pa * ECU_VD_CYL / (ECU_R_AIR * in.iat_k);
    double tau = ECU_FILM_TAU_REF * pow(out.rpm / ECU_FILM_RPM_REF, ECU_FILM_TAU_EXP);
    double e2 = exp(-ECU_EOI_TO_IVC_DEG / 360.0 * 60.0 / out.rpm / tau);
    double m_inj = air / ECU_AFR_STOICH / (1.0 - ECU_FILM_X * e2);   /* empty film */
    double dead = ecu_interp1(ECU_DEAD_VBATT, ECU_DEAD_S, ECU_DEAD_VBATT_N, in.vbatt);
    CHECK_NEAR(out.pw_s[0], m_inj / ECU_INJ_FLOW + dead, 1e-7);
    CHECK(e.sch.pw_ticks[0] > 0u);
    CHECK(out.closed_loop == 0.0f);
    /* throttle: target from the pedal map, PID pushes the plate open */
    CHECK(out.etc_voltage > 0.0f);
    CHECK_NEAR(out.throttle_target_deg, ecu_interp1(ECU_PEDAL_X, ECU_PEDAL_DEG, ECU_PEDAL_X_N, 0.3f), 1e-5);
}

static void test_fault_debounce_and_limp(void)
{
    hil_reset();
    ecu_t e;
    ecu_init(&e);
    ecu_inputs_t in = nominal_inputs(0u);
    ecu_outputs_t out;
    in.tps2_deg = 40.0f;             /* sensors disagree */
    for (int i = 0; i < ECU_FAULT_DEBOUNCE - 1; ++i) ecu_slow_step(&e, &in, &out);
    CHECK(out.limp == 0.0f);
    ecu_slow_step(&e, &in, &out);
    CHECK(out.limp == 1.0f);
    CHECK(((uint32_t)out.faults & ECU_FAULT_TPS) != 0u);
    CHECK(out.etc_voltage == 0.0f);
    CHECK(e.sch.limp == 1u);
    /* MAP out of range: alpha-N substitute */
    in.map_pa = 0.0f;
    for (int i = 0; i < ECU_FAULT_DEBOUNCE; ++i) ecu_slow_step(&e, &in, &out);
    CHECK(((uint32_t)out.faults & ECU_FAULT_MAP) != 0u);
    CHECK(out.map_used > ECU_MAP_MIN);
}

int main(void)
{
    RUN(test_sync_within_one_rev_after_gap);
    RUN(test_cam_phase_and_rpm);
    RUN(test_noise_rejection);
    RUN(test_single_tooth_dropout);
    RUN(test_loss_of_sync_and_recovery);
    RUN(test_stall_and_timer_wrap);
    RUN(test_scheduler_accuracy_7000rpm);
    RUN(test_rev_limiter);
    RUN(test_knock_retard_and_recover);
    RUN(test_speed_density_pulse_width);
    RUN(test_fault_debounce_and_limp);
    return TH_SUMMARY();
}
