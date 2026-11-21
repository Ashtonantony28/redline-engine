/* Division by multiplication for the RV64 target (M9.1, ADR 0005). Vortex-V multiplies in one
 * cycle but divides one quotient bit per clock (rtl/core/vortex_muldiv.sv: 64 iterations), so the
 * target's integer-division helpers (rv_div.c) and its soft-float divisions (rv_softfp.c) use a
 * Newton-Raphson reciprocal and finish with an exact remainder correction: the quotient is the
 * mathematically exact one, bit for bit what a `divu` would return.
 *
 * rv_recip(d), d normalised (bit 63 set, i.e. d / 2^64 in [1/2, 1)), returns X with
 *     2^126 / d - RV_RECIP_SLACK  <=  X  <=  2^126 / d                      (never above)
 * Derivation: x' = x (2 - d x) converges quadratically from below; the seed 48/17 - 32/17 d has
 * relative error <= 1/17, four steps leave < 2^-66. Each step truncates twice: the product d x
 * (which can raise x' by at most 2 units of 2^-62) and x' itself (which lowers it by < 4), so the
 * last iterate lies in [2^126/d - 7, 2^126/d + 2]; subtracting 2 gives the bound above. */
#ifndef RV_RECIP_H
#define RV_RECIP_H

#include <stdint.h>

#define RV_RECIP_SEED_A 0xB4B4B4B4B4B4B4B4ull   /* 48/17 x 2^62 */
#define RV_RECIP_SEED_B 0x7878787878787878ull   /* 32/17 x 2^62 */
#define RV_RECIP_TWO 0x8000000000000000ull      /* 2.0 in Q2.62 */
#define RV_RECIP_STEPS 4
#define RV_RECIP_OVERSHOOT 2u
#define RV_RECIP_SLACK 9u
#define RV_RECIP_FRAC_BITS 62

static inline uint64_t rv_mulhu(uint64_t a, uint64_t b)
{
    return (uint64_t)(((__uint128_t)a * b) >> 64);
}

static inline uint64_t rv_recip(uint64_t d)
{
    uint64_t x = RV_RECIP_SEED_A - rv_mulhu(RV_RECIP_SEED_B, d);
    for (int i = 0; i < RV_RECIP_STEPS; ++i) {
        x = rv_mulhu(x, RV_RECIP_TWO - rv_mulhu(d, x)) << 2;
    }
    return x - RV_RECIP_OVERSHOOT;
}

/* Count leading zeros, x != 0 (no Zbb on Vortex-V; libgcc's table-driven __clzdi2 is a call). */
static inline unsigned rv_clz64(uint64_t x)
{
    unsigned n = 0;
    if (!(x >> 32)) { n += 32; x <<= 32; }
    if (!(x >> 48)) { n += 16; x <<= 16; }
    if (!(x >> 56)) { n += 8; x <<= 8; }
    if (!(x >> 60)) { n += 4; x <<= 4; }
    if (!(x >> 62)) { n += 2; x <<= 2; }
    if (!(x >> 63)) { n += 1; }
    return n;
}

static inline unsigned rv_clz32(uint32_t x)
{
    return rv_clz64((uint64_t)x) - 32u;
}

#endif
