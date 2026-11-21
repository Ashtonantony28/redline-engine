/* Integer division helpers for the sync-path units of the RV64 target (M9.1, ADR 0005).
 *
 * decoder.c / sched.c / knock.c are compiled with -mno-div: GCC then turns every division by a
 * constant (ecu_angle_wrap's `% ECU_CYCLE_U`, the event-tick interpolation) into a multiply, and
 * every division by a variable into a call to __udivdi3 / __divdi3, which are defined here instead
 * of in libgcc. On Vortex-V a `divu` costs 64 cycles and M9's crank-tooth handler spent about
 * 1300 of its 1430 cycles in 20 of them. Here a quotient is a reciprocal multiply plus an exact
 * correction (rv_recip.h), and the reciprocal of the last divisor is kept: the scheduler divides
 * by the same tooth period (period_q8) five times per tooth.
 *
 * Exact: the results equal RISC-V div / divu for every operand pair, including the ISA's
 * divide-by-zero and overflow results (checked on the host in tests/softfp_check.c).
 * The one-entry cache makes these functions non-reentrant. They are called only from the
 * crank/cam/knock capture context (interrupts masked in the handler) and from the frozen-time
 * command loop while interrupts are off; the slow loop's units are built with hardware division. */
#include <stdint.h>

#include "rv_recip.h"

#ifndef RV_DIV_NAME
#define RV_DIV_NAME(x) __##x
#endif

static uint64_t g_div_d;        /* cached divisor (0 = none) */
static uint64_t g_div_x;        /* rv_recip(g_div_d << g_div_sh) */
static unsigned g_div_sh;

uint64_t RV_DIV_NAME(udivdi3)(uint64_t n, uint64_t d)
{
    if (d <= 1u) {
        return d ? n : ~(uint64_t)0;            /* RISC-V divu by zero: all ones */
    }
    if (d != g_div_d) {
        g_div_sh = rv_clz64(d);
        g_div_x = rv_recip(d << g_div_sh);
        g_div_d = d;
    }
    /* n / d = n 2^sh / (d 2^sh) ~ n x / 2^(126 - sh); x never overestimates, so r stays >= 0 */
    const unsigned shift = RV_RECIP_FRAC_BITS - g_div_sh;      /* d >= 2: sh <= 62 */
    uint64_t q = 0u;
    uint64_t r = n;
    while (r >= d) {
        uint64_t t = rv_mulhu(r, g_div_x) >> shift;
        if (t == 0u) {
            t = 1u;
        }
        q += t;
        r -= t * d;
    }
    return q;
}

int64_t RV_DIV_NAME(divdi3)(int64_t n, int64_t d)
{
    const uint64_t un = n < 0 ? 0u - (uint64_t)n : (uint64_t)n;
    const uint64_t ud = d < 0 ? 0u - (uint64_t)d : (uint64_t)d;
    if (ud == 0u) {
        return -1;                               /* RISC-V div by zero */
    }
    const uint64_t q = RV_DIV_NAME(udivdi3)(un, ud);
    return (int64_t)(((n < 0) != (d < 0)) ? 0u - q : q);       /* INT64_MIN / -1 wraps to INT64_MIN */
}
