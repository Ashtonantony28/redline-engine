/* IEEE-754 binary32 / binary64 soft-float for the RV64 target (M9.1, ADR 0005): the libgcc entry
 * points the ECU slow loop and ecu_libm.c reference, rewritten for a 64-bit integer core with a
 * one-cycle multiplier and a 64-cycle divider. Linked ahead of libgcc, it replaces libgcc's
 * generic soft-fp (which spent 15.7k of the slow loop's 44k cycles in __divdf3 alone).
 *
 * Results are the correctly rounded IEEE-754 ones (round to nearest, ties to even; no exception
 * flags, no other rounding mode: the ECU uses neither), so the slow loop stays BIT-EXACT against
 * the host build, whose SSE2 arithmetic is correctly rounded too. tests/softfp_check.c compiles
 * this file on the host and compares every operation with the FPU on special-value pairs and
 * seeded random operands; the frozen-time co-simulation compares whole slow-loop outputs.
 *
 * Conventions (the structure follows Berkeley SoftFloat 3; the divisions are rv_recip.h's):
 *   significands carry the leading one at bit 30 (binary32) / bit 62 (binary64) into round-pack,
 *   with 7 / 10 rounding bits below the result's last place and a sticky bit ORed into bit 0;
 *   pack(sign, exp, sig) adds sig to the exponent field, so a carry out of the fraction bumps it.
 * NaN: a NaN operand is returned quieted (first operand first); an invalid operation returns the
 *   RISC-V canonical NaN. NaN payloads are the only results not compared with the host (x86's
 *   default NaN has the sign bit set); the ECU never produces one.
 * Float -> integer conversion truncates and saturates out of range, as libgcc does. */
#include <stdint.h>

#include "rv_recip.h"

#ifndef SOFTFP_NAME
#define SOFTFP_NAME(x) __##x
#endif

/* Every helper is inlined into its entry point: on Vortex-V a return (JALR) stops dispatch until
 * it resolves, so a nested call costs a pipeline drain. */
#define SOFTFP_INLINE static inline __attribute__((always_inline))

typedef union { float f; uint32_t u; } f32_t;
typedef union { double f; uint64_t u; } f64_t;

/* ---- binary32 layout ---- */
#define F32_SIGN 0x80000000u
#define F32_EXP_MAX 0xFF
#define F32_FRAC_BITS 23
#define F32_FRAC_MASK 0x007FFFFFu
#define F32_IMPLICIT 0x00800000u
#define F32_INF 0x7F800000u
#define F32_QUIET 0x00400000u
#define F32_DEFAULT_NAN 0x7FC00000u
#define F32_EXP_BIAS 0x7F
#define F32_ROUND_BITS 7
#define F32_ROUND_MASK 0x7Fu
#define F32_ROUND_HALF 0x40u
#define F32_SIG_TOP 0x40000000u         /* leading one of a round-pack significand */
#define F32_EXP_LIMIT 0xFD              /* round-pack exponent above which the result overflows */
#define F32_INT_EXP 0x9E                /* biased exponent of 2^31 */

/* ---- binary64 layout ---- */
#define F64_SIGN 0x8000000000000000ull
#define F64_EXP_MAX 0x7FF
#define F64_FRAC_BITS 52
#define F64_FRAC_MASK 0x000FFFFFFFFFFFFFull
#define F64_IMPLICIT 0x0010000000000000ull
#define F64_INF 0x7FF0000000000000ull
#define F64_QUIET 0x0008000000000000ull
#define F64_DEFAULT_NAN 0x7FF8000000000000ull
#define F64_EXP_BIAS 0x3FF
#define F64_ROUND_BITS 10
#define F64_ROUND_MASK 0x3FFull
#define F64_ROUND_HALF 0x200ull
#define F64_SIG_TOP 0x4000000000000000ull
#define F64_EXP_LIMIT 0x7FD
#define F64_INT_EXP 0x433               /* biased exponent at which the significand is an integer */

#define I32_BITS 32

static inline uint32_t f32_bits(float f) { f32_t v; v.f = f; return v.u; }
static inline float f32_from(uint32_t u) { f32_t v; v.u = u; return v.f; }
static inline uint64_t f64_bits(double f) { f64_t v; v.f = f; return v.u; }
static inline double f64_from(uint64_t u) { f64_t v; v.u = u; return v.f; }

static inline int f32_exp(uint32_t a) { return (int)((a >> F32_FRAC_BITS) & F32_EXP_MAX); }
static inline uint32_t f32_frac(uint32_t a) { return a & F32_FRAC_MASK; }
static inline int f32_is_nan(uint32_t a) { return (a & ~F32_SIGN) > F32_INF; }
static inline uint32_t f32_pack(uint32_t sign, int exp, uint32_t sig)
{
    return sign + ((uint32_t)exp << F32_FRAC_BITS) + sig;
}

static inline int f64_exp(uint64_t a) { return (int)((a >> F64_FRAC_BITS) & F64_EXP_MAX); }
static inline uint64_t f64_frac(uint64_t a) { return a & F64_FRAC_MASK; }
static inline int f64_is_nan(uint64_t a) { return (a & ~F64_SIGN) > F64_INF; }
static inline uint64_t f64_pack(uint64_t sign, int exp, uint64_t sig)
{
    return sign + ((uint64_t)exp << F64_FRAC_BITS) + sig;
}

static inline uint32_t f32_nan2(uint32_t a, uint32_t b)
{
    return (f32_is_nan(a) ? a : b) | F32_QUIET;
}

static inline uint64_t f64_nan2(uint64_t a, uint64_t b)
{
    return (f64_is_nan(a) ? a : b) | F64_QUIET;
}

/* x >> n with every bit shifted out ORed into bit 0 (n >= 1) */
static inline uint32_t shr_jam32(uint32_t x, unsigned n)
{
    return n < I32_BITS - 1 ? (x >> n) | (uint32_t)((uint32_t)(x << (I32_BITS - n)) != 0u) : (uint32_t)(x != 0u);
}

static inline uint64_t shr_jam64(uint64_t x, unsigned n)
{
    return n < 63u ? (x >> n) | (uint64_t)((x << (64u - n)) != 0u) : (uint64_t)(x != 0u);
}

/* ---- round and pack ---- */
SOFTFP_INLINE uint32_t f32_round_pack(uint32_t sign, int exp, uint32_t sig)
{
    if ((unsigned)exp >= F32_EXP_LIMIT) {
        if (exp < 0) {
            sig = shr_jam32(sig, (unsigned)-exp);            /* subnormal result */
            exp = 0;
        } else if (exp > F32_EXP_LIMIT || sig + F32_ROUND_HALF >= F32_SIGN) {
            return sign | F32_INF;
        }
    }
    const uint32_t round = sig & F32_ROUND_MASK;
    sig = (sig + F32_ROUND_HALF) >> F32_ROUND_BITS;
    if (round == F32_ROUND_HALF) {
        sig &= ~(uint32_t)1;                                  /* tie: to even */
    }
    if (!sig) {
        exp = 0;
    }
    return f32_pack(sign, exp, sig);
}

SOFTFP_INLINE uint32_t f32_norm_round_pack(uint32_t sign, int exp, uint32_t sig)
{
    const int shift = (int)rv_clz32(sig) - 1;                 /* sig != 0 */
    exp -= shift;
    if (shift >= F32_ROUND_BITS && (unsigned)exp < F32_EXP_LIMIT) {
        return f32_pack(sign, exp, sig << (shift - F32_ROUND_BITS));      /* exact */
    }
    return f32_round_pack(sign, exp, sig << shift);
}

SOFTFP_INLINE uint64_t f64_round_pack(uint64_t sign, int exp, uint64_t sig)
{
    if ((unsigned)exp >= F64_EXP_LIMIT) {
        if (exp < 0) {
            sig = shr_jam64(sig, (unsigned)-exp);
            exp = 0;
        } else if (exp > F64_EXP_LIMIT || sig + F64_ROUND_HALF >= F64_SIGN) {
            return sign | F64_INF;
        }
    }
    const uint64_t round = sig & F64_ROUND_MASK;
    sig = (sig + F64_ROUND_HALF) >> F64_ROUND_BITS;
    if (round == F64_ROUND_HALF) {
        sig &= ~(uint64_t)1;
    }
    if (!sig) {
        exp = 0;
    }
    return f64_pack(sign, exp, sig);
}

SOFTFP_INLINE uint64_t f64_norm_round_pack(uint64_t sign, int exp, uint64_t sig)
{
    const int shift = (int)rv_clz64(sig) - 1;
    exp -= shift;
    if (shift >= F64_ROUND_BITS && (unsigned)exp < F64_EXP_LIMIT) {
        return f64_pack(sign, exp, sig << (shift - F64_ROUND_BITS));
    }
    return f64_round_pack(sign, exp, sig << shift);
}

/* ---- binary32 add / subtract ---- */
SOFTFP_INLINE uint32_t f32_add_mags(uint32_t a, uint32_t b)
{
    int ea = f32_exp(a), eb = f32_exp(b);
    uint32_t sa = f32_frac(a), sb = f32_frac(b);
    const uint32_t sign = a & F32_SIGN;
    const int diff = ea - eb;
    int ez;
    uint32_t sz;
    if (!diff) {
        if (!ea) {
            return a + sb;                                    /* subnormals: plain integer add */
        }
        if (ea == F32_EXP_MAX) {
            return (sa | sb) ? f32_nan2(a, b) : a;
        }
        ez = ea;
        sz = 2u * F32_IMPLICIT + sa + sb;
        if (!(sz & 1u) && ez < F32_EXP_MAX - 1) {
            return f32_pack(sign, ez, sz >> 1);               /* exact */
        }
        sz <<= F32_ROUND_BITS - 1;
    } else {
        sa <<= F32_ROUND_BITS - 1;
        sb <<= F32_ROUND_BITS - 1;
        if (diff < 0) {
            if (eb == F32_EXP_MAX) {
                return sb ? f32_nan2(a, b) : (sign | F32_INF);
            }
            ez = eb;
            sa += ea ? F32_SIG_TOP >> 1 : sa;
            sa = shr_jam32(sa, (unsigned)-diff);
        } else {
            if (ea == F32_EXP_MAX) {
                return sa ? f32_nan2(a, b) : a;
            }
            ez = ea;
            sb += eb ? F32_SIG_TOP >> 1 : sb;
            sb = shr_jam32(sb, (unsigned)diff);
        }
        sz = (F32_SIG_TOP >> 1) + sa + sb;
        if (sz < F32_SIG_TOP) {
            --ez;
            sz <<= 1;
        }
    }
    return f32_round_pack(sign, ez, sz);
}

SOFTFP_INLINE uint32_t f32_sub_mags(uint32_t a, uint32_t b)
{
    int ea = f32_exp(a), eb = f32_exp(b);
    uint32_t sa = f32_frac(a), sb = f32_frac(b);
    uint32_t sign = a & F32_SIGN;
    const int diff = ea - eb;
    if (!diff) {
        if (ea == F32_EXP_MAX) {
            return (sa | sb) ? f32_nan2(a, b) : F32_DEFAULT_NAN;          /* inf - inf */
        }
        int32_t d = (int32_t)(sa - sb);
        if (!d) {
            return 0u;                                        /* +0 in round-to-nearest */
        }
        if (ea) {
            --ea;
        }
        if (d < 0) {
            sign ^= F32_SIGN;
            d = -d;
        }
        int shift = (int)rv_clz32((uint32_t)d) - (I32_BITS - 1 - F32_FRAC_BITS);
        int ez = ea - shift;
        if (ez < 0) {
            shift = ea;
            ez = 0;
        }
        return f32_pack(sign, ez, (uint32_t)d << shift);
    }
    int ez;
    uint32_t x, y;
    sa <<= F32_ROUND_BITS;
    sb <<= F32_ROUND_BITS;
    if (diff < 0) {
        sign ^= F32_SIGN;
        if (eb == F32_EXP_MAX) {
            return sb ? f32_nan2(a, b) : (sign | F32_INF);
        }
        ez = eb - 1;
        x = sb | F32_SIG_TOP;
        y = sa + (ea ? F32_SIG_TOP : sa);
        y = shr_jam32(y, (unsigned)-diff);
    } else {
        if (ea == F32_EXP_MAX) {
            return sa ? f32_nan2(a, b) : a;
        }
        ez = ea - 1;
        x = sa | F32_SIG_TOP;
        y = sb + (eb ? F32_SIG_TOP : sb);
        y = shr_jam32(y, (unsigned)diff);
    }
    return f32_norm_round_pack(sign, ez, x - y);
}

/* the magnitude routines take the sign from a and ignore b's */
float SOFTFP_NAME(addsf3)(float fa, float fb)
{
    const uint32_t a = f32_bits(fa), b = f32_bits(fb);
    return f32_from(((a ^ b) & F32_SIGN) ? f32_sub_mags(a, b) : f32_add_mags(a, b));
}

float SOFTFP_NAME(subsf3)(float fa, float fb)
{
    const uint32_t a = f32_bits(fa), b = f32_bits(fb);
    return f32_from(((a ^ b) & F32_SIGN) ? f32_add_mags(a, b) : f32_sub_mags(a, b));
}

/* subnormal significand -> normalised, returns the matching exponent */
static inline int f32_norm_subnormal(uint32_t *sig)
{
    const int shift = (int)rv_clz32(*sig) - (I32_BITS - 1 - F32_FRAC_BITS);
    *sig <<= shift;
    return 1 - shift;
}

float SOFTFP_NAME(mulsf3)(float fa, float fb)
{
    const uint32_t a = f32_bits(fa), b = f32_bits(fb);
    int ea = f32_exp(a), eb = f32_exp(b);
    uint32_t sa = f32_frac(a), sb = f32_frac(b);
    const uint32_t sign = (a ^ b) & F32_SIGN;
    if (ea == F32_EXP_MAX || eb == F32_EXP_MAX) {
        if (f32_is_nan(a) || f32_is_nan(b)) {
            return f32_from(f32_nan2(a, b));
        }
        const uint32_t other = (ea == F32_EXP_MAX) ? (b & ~F32_SIGN) : (a & ~F32_SIGN);
        return f32_from(other ? (sign | F32_INF) : F32_DEFAULT_NAN);      /* inf x 0 */
    }
    if (!ea) {
        if (!sa) {
            return f32_from(sign);
        }
        ea = f32_norm_subnormal(&sa);
    }
    if (!eb) {
        if (!sb) {
            return f32_from(sign);
        }
        eb = f32_norm_subnormal(&sb);
    }
    int ez = ea + eb - F32_EXP_BIAS;
    sa = (sa | F32_IMPLICIT) << F32_ROUND_BITS;
    sb = (sb | F32_IMPLICIT) << (F32_ROUND_BITS + 1);
    const uint64_t p = (uint64_t)sa * sb;
    uint32_t sz = (uint32_t)(p >> I32_BITS) | (uint32_t)((uint32_t)p != 0u);
    if (sz < F32_SIG_TOP) {
        --ez;
        sz <<= 1;
    }
    return f32_from(f32_round_pack(sign, ez, sz));
}

float SOFTFP_NAME(divsf3)(float fa, float fb)
{
    const uint32_t a = f32_bits(fa), b = f32_bits(fb);
    int ea = f32_exp(a), eb = f32_exp(b);
    uint32_t sa = f32_frac(a), sb = f32_frac(b);
    const uint32_t sign = (a ^ b) & F32_SIGN;
    if (ea == F32_EXP_MAX) {
        if (sa || f32_is_nan(b)) {
            return f32_from(f32_nan2(a, b));
        }
        return f32_from(eb == F32_EXP_MAX ? F32_DEFAULT_NAN : (sign | F32_INF));   /* inf / inf */
    }
    if (eb == F32_EXP_MAX) {
        return f32_from(sb ? f32_nan2(a, b) : sign);
    }
    if (!eb) {
        if (!sb) {
            return f32_from((ea | (int)sa) ? (sign | F32_INF) : F32_DEFAULT_NAN);  /* x / 0, 0 / 0 */
        }
        eb = f32_norm_subnormal(&sb);
    }
    if (!ea) {
        if (!sa) {
            return f32_from(sign);
        }
        ea = f32_norm_subnormal(&sa);
    }
    int ez = ea - eb + F32_EXP_BIAS - 1;
    sa |= F32_IMPLICIT;
    sb |= F32_IMPLICIT;
    /* q = floor(sa 2^s / sb) with its leading one at bit 30: s = 31 if sa < sb, else 30.
     * x ~ 2^126 / (sb 2^40) = 2^86 / sb, so q ~ (sa << (s - 22)) x / 2^64, never above. */
    unsigned s = I32_BITS - 2;
    if (sa < sb) {
        --ez;
        s = I32_BITS - 1;
    }
    const uint64_t x = rv_recip((uint64_t)sb << (64 - 1 - F32_FRAC_BITS));
    uint64_t q = rv_mulhu((uint64_t)sa << (s - (RV_RECIP_FRAC_BITS - 64 + 1 + F32_FRAC_BITS)), x);
    uint64_t r = ((uint64_t)sa << s) - q * sb;
    while (r >= sb) {
        ++q;
        r -= sb;
    }
    return f32_from(f32_round_pack(sign, ez, (uint32_t)q | (uint32_t)(r != 0u)));
}

/* ---- binary64 add / subtract ---- */
SOFTFP_INLINE uint64_t f64_add_mags(uint64_t a, uint64_t b)
{
    int ea = f64_exp(a), eb = f64_exp(b);
    uint64_t sa = f64_frac(a), sb = f64_frac(b);
    const uint64_t sign = a & F64_SIGN;
    const int diff = ea - eb;
    int ez;
    uint64_t sz;
    if (!diff) {
        if (!ea) {
            return a + sb;
        }
        if (ea == F64_EXP_MAX) {
            return (sa | sb) ? f64_nan2(a, b) : a;
        }
        ez = ea;
        sz = (2u * F64_IMPLICIT + sa + sb) << (F64_ROUND_BITS - 1);
    } else {
        sa <<= F64_ROUND_BITS - 1;
        sb <<= F64_ROUND_BITS - 1;
        if (diff < 0) {
            if (eb == F64_EXP_MAX) {
                return sb ? f64_nan2(a, b) : (sign | F64_INF);
            }
            ez = eb;
            if (ea) {
                sa += F64_SIG_TOP >> 1;
            } else {
                sa <<= 1;
            }
            sa = shr_jam64(sa, (unsigned)-diff);
        } else {
            if (ea == F64_EXP_MAX) {
                return sa ? f64_nan2(a, b) : a;
            }
            ez = ea;
            if (eb) {
                sb += F64_SIG_TOP >> 1;
            } else {
                sb <<= 1;
            }
            sb = shr_jam64(sb, (unsigned)diff);
        }
        sz = (F64_SIG_TOP >> 1) + sa + sb;
        if (sz < F64_SIG_TOP) {
            --ez;
            sz <<= 1;
        }
    }
    return f64_round_pack(sign, ez, sz);
}

SOFTFP_INLINE uint64_t f64_sub_mags(uint64_t a, uint64_t b)
{
    int ea = f64_exp(a), eb = f64_exp(b);
    uint64_t sa = f64_frac(a), sb = f64_frac(b);
    uint64_t sign = a & F64_SIGN;
    const int diff = ea - eb;
    if (!diff) {
        if (ea == F64_EXP_MAX) {
            return (sa | sb) ? f64_nan2(a, b) : F64_DEFAULT_NAN;
        }
        int64_t d = (int64_t)(sa - sb);
        if (!d) {
            return 0u;
        }
        if (ea) {
            --ea;
        }
        if (d < 0) {
            sign ^= F64_SIGN;
            d = -d;
        }
        int shift = (int)rv_clz64((uint64_t)d) - (64 - 1 - F64_FRAC_BITS);
        int ez = ea - shift;
        if (ez < 0) {
            shift = ea;
            ez = 0;
        }
        return f64_pack(sign, ez, (uint64_t)d << shift);
    }
    int ez;
    uint64_t sz;
    sa <<= F64_ROUND_BITS;
    sb <<= F64_ROUND_BITS;
    if (diff < 0) {
        sign ^= F64_SIGN;
        if (eb == F64_EXP_MAX) {
            return sb ? f64_nan2(a, b) : (sign | F64_INF);
        }
        sa += ea ? F64_SIG_TOP : sa;
        sa = shr_jam64(sa, (unsigned)-diff);
        sb |= F64_SIG_TOP;
        ez = eb;
        sz = sb - sa;
    } else {
        if (ea == F64_EXP_MAX) {
            return sa ? f64_nan2(a, b) : a;
        }
        sb += eb ? F64_SIG_TOP : sb;
        sb = shr_jam64(sb, (unsigned)diff);
        sa |= F64_SIG_TOP;
        ez = ea;
        sz = sa - sb;
    }
    return f64_norm_round_pack(sign, ez - 1, sz);
}

double SOFTFP_NAME(adddf3)(double fa, double fb)
{
    const uint64_t a = f64_bits(fa), b = f64_bits(fb);
    return f64_from(((a ^ b) & F64_SIGN) ? f64_sub_mags(a, b) : f64_add_mags(a, b));
}

double SOFTFP_NAME(subdf3)(double fa, double fb)
{
    const uint64_t a = f64_bits(fa), b = f64_bits(fb);
    return f64_from(((a ^ b) & F64_SIGN) ? f64_add_mags(a, b) : f64_sub_mags(a, b));
}

static inline int f64_norm_subnormal(uint64_t *sig)
{
    const int shift = (int)rv_clz64(*sig) - (64 - 1 - F64_FRAC_BITS);
    *sig <<= shift;
    return 1 - shift;
}

double SOFTFP_NAME(muldf3)(double fa, double fb)
{
    const uint64_t a = f64_bits(fa), b = f64_bits(fb);
    int ea = f64_exp(a), eb = f64_exp(b);
    uint64_t sa = f64_frac(a), sb = f64_frac(b);
    const uint64_t sign = (a ^ b) & F64_SIGN;
    if (ea == F64_EXP_MAX || eb == F64_EXP_MAX) {
        if (f64_is_nan(a) || f64_is_nan(b)) {
            return f64_from(f64_nan2(a, b));
        }
        const uint64_t other = (ea == F64_EXP_MAX) ? (b & ~F64_SIGN) : (a & ~F64_SIGN);
        return f64_from(other ? (sign | F64_INF) : F64_DEFAULT_NAN);
    }
    if (!ea) {
        if (!sa) {
            return f64_from(sign);
        }
        ea = f64_norm_subnormal(&sa);
    }
    if (!eb) {
        if (!sb) {
            return f64_from(sign);
        }
        eb = f64_norm_subnormal(&sb);
    }
    int ez = ea + eb - F64_EXP_BIAS;
    sa = (sa | F64_IMPLICIT) << F64_ROUND_BITS;
    sb = (sb | F64_IMPLICIT) << (F64_ROUND_BITS + 1);
    const __uint128_t p = (__uint128_t)sa * sb;
    uint64_t sz = (uint64_t)(p >> 64) | (uint64_t)((uint64_t)p != 0u);
    if (sz < F64_SIG_TOP) {
        --ez;
        sz <<= 1;
    }
    return f64_from(f64_round_pack(sign, ez, sz));
}

#define F64_DIV_Q_TOP 54        /* leading one of the divide's quotient: 53 result bits + guard */

double SOFTFP_NAME(divdf3)(double fa, double fb)
{
    const uint64_t a = f64_bits(fa), b = f64_bits(fb);
    int ea = f64_exp(a), eb = f64_exp(b);
    uint64_t sa = f64_frac(a), sb = f64_frac(b);
    const uint64_t sign = (a ^ b) & F64_SIGN;
    if (ea == F64_EXP_MAX) {
        if (sa || f64_is_nan(b)) {
            return f64_from(f64_nan2(a, b));
        }
        return f64_from(eb == F64_EXP_MAX ? F64_DEFAULT_NAN : (sign | F64_INF));
    }
    if (eb == F64_EXP_MAX) {
        return f64_from(sb ? f64_nan2(a, b) : sign);
    }
    if (!eb) {
        if (!sb) {
            return f64_from((ea || sa) ? (sign | F64_INF) : F64_DEFAULT_NAN);
        }
        eb = f64_norm_subnormal(&sb);
    }
    if (!ea) {
        if (!sa) {
            return f64_from(sign);
        }
        ea = f64_norm_subnormal(&sa);
    }
    int ez = ea - eb + F64_EXP_BIAS - 1;
    sa |= F64_IMPLICIT;
    sb |= F64_IMPLICIT;
    /* q = floor(sa 2^s / sb) with its leading one at bit 54 (s = 55 if sa < sb, else 54): the 53
     * result bits, a guard bit and one more; everything below is the sticky bit r != 0.
     * x ~ 2^126 / (sb 2^11) = 2^115 / sb, so q ~ (sa << (s - 51)) x / 2^64, never above. The
     * remainder is computed modulo 2^64: its true value is below a few sb. */
    unsigned s = F64_DIV_Q_TOP;
    if (sa < sb) {
        --ez;
        s = F64_DIV_Q_TOP + 1;
    }
    const uint64_t x = rv_recip(sb << (64 - 1 - F64_FRAC_BITS));
    uint64_t q = rv_mulhu(sa << (s - (RV_RECIP_FRAC_BITS - 64 + 1 + F64_FRAC_BITS)), x);
    uint64_t r = (sa << s) - q * sb;
    while (r >= sb) {
        ++q;
        r -= sb;
    }
    return f64_from(f64_round_pack(sign, ez, (q << (62 - F64_DIV_Q_TOP)) | (uint64_t)(r != 0u)));
}

/* ---- comparisons: -1 / 0 / +1, `unordered` when either operand is a NaN ---- */
SOFTFP_INLINE int f32_cmp(float fa, float fb, int unordered)
{
    const uint32_t a = f32_bits(fa), b = f32_bits(fb);
    if (f32_is_nan(a) || f32_is_nan(b)) {
        return unordered;
    }
    if (a == b || !((a | b) & ~F32_SIGN)) {
        return 0;                                             /* equal, or +0 vs -0 */
    }
    const int neg_a = (a & F32_SIGN) != 0u, neg_b = (b & F32_SIGN) != 0u;
    if (neg_a != neg_b) {
        return neg_a ? -1 : 1;
    }
    return ((a < b) != neg_a) ? -1 : 1;
}

SOFTFP_INLINE int f64_cmp(double fa, double fb, int unordered)
{
    const uint64_t a = f64_bits(fa), b = f64_bits(fb);
    if (f64_is_nan(a) || f64_is_nan(b)) {
        return unordered;
    }
    if (a == b || !((a | b) & ~F64_SIGN)) {
        return 0;
    }
    const int neg_a = (a & F64_SIGN) != 0u, neg_b = (b & F64_SIGN) != 0u;
    if (neg_a != neg_b) {
        return neg_a ? -1 : 1;
    }
    return ((a < b) != neg_a) ? -1 : 1;
}

/* libgcc contract: the caller tests the result against 0 with the operator in the name, and a
 * NaN operand must make every test but != come out false. */
int SOFTFP_NAME(eqsf2)(float a, float b) { return f32_cmp(a, b, 1); }
int SOFTFP_NAME(nesf2)(float a, float b) { return f32_cmp(a, b, 1); }
int SOFTFP_NAME(ltsf2)(float a, float b) { return f32_cmp(a, b, 1); }
int SOFTFP_NAME(lesf2)(float a, float b) { return f32_cmp(a, b, 1); }
int SOFTFP_NAME(gtsf2)(float a, float b) { return f32_cmp(a, b, -1); }
int SOFTFP_NAME(gesf2)(float a, float b) { return f32_cmp(a, b, -1); }
int SOFTFP_NAME(unordsf2)(float a, float b) { return f32_is_nan(f32_bits(a)) || f32_is_nan(f32_bits(b)); }
int SOFTFP_NAME(eqdf2)(double a, double b) { return f64_cmp(a, b, 1); }
int SOFTFP_NAME(nedf2)(double a, double b) { return f64_cmp(a, b, 1); }
int SOFTFP_NAME(ltdf2)(double a, double b) { return f64_cmp(a, b, 1); }
int SOFTFP_NAME(ledf2)(double a, double b) { return f64_cmp(a, b, 1); }
int SOFTFP_NAME(gtdf2)(double a, double b) { return f64_cmp(a, b, -1); }
int SOFTFP_NAME(gedf2)(double a, double b) { return f64_cmp(a, b, -1); }
int SOFTFP_NAME(unorddf2)(double a, double b) { return f64_is_nan(f64_bits(a)) || f64_is_nan(f64_bits(b)); }

/* ---- conversions ---- */
double SOFTFP_NAME(extendsfdf2)(float fa)
{
    const uint32_t a = f32_bits(fa);
    int e = f32_exp(a);
    uint32_t s = f32_frac(a);
    const uint64_t sign = (uint64_t)(a & F32_SIGN) << I32_BITS;
    if (e == F32_EXP_MAX) {
        return f64_from(sign | F64_INF | (s ? ((uint64_t)s << (F64_FRAC_BITS - F32_FRAC_BITS)) | F64_QUIET : 0u));
    }
    if (!e) {
        if (!s) {
            return f64_from(sign);
        }
        e = f32_norm_subnormal(&s) - 1;
        s &= F32_FRAC_MASK;
        return f64_from(f64_pack(sign, e + 1 + (F64_EXP_BIAS - F32_EXP_BIAS), (uint64_t)s << (F64_FRAC_BITS - F32_FRAC_BITS)));
    }
    return f64_from(f64_pack(sign, e + (F64_EXP_BIAS - F32_EXP_BIAS), (uint64_t)s << (F64_FRAC_BITS - F32_FRAC_BITS)));
}

float SOFTFP_NAME(truncdfsf2)(double fa)
{
    const uint64_t a = f64_bits(fa);
    const int e = f64_exp(a);
    const uint64_t s = f64_frac(a);
    const uint32_t sign = (uint32_t)(a >> I32_BITS) & F32_SIGN;
    if (e == F64_EXP_MAX) {
        return f32_from(sign | F32_INF | (s ? (uint32_t)(s >> (F64_FRAC_BITS - F32_FRAC_BITS)) | F32_QUIET : 0u));
    }
    /* 52 fraction bits -> leading one at bit 30: keep 30, jam the other 22 */
    const unsigned drop = F64_FRAC_BITS - (I32_BITS - 2);
    const uint32_t s32 = (uint32_t)(s >> drop) | (uint32_t)((s & (((uint64_t)1 << drop) - 1u)) != 0u);
    if (!e && !s32) {
        return f32_from(sign);
    }
    return f32_from(f32_round_pack(sign, e - (F64_EXP_BIAS - F32_EXP_BIAS) - 1, s32 | F32_SIG_TOP));
}

float SOFTFP_NAME(floatunsisf)(uint32_t a)
{
    if (!a) {
        return f32_from(0u);
    }
    if (a & F32_SIGN) {
        return f32_from(f32_round_pack(0u, F32_INT_EXP - 1, (a >> 1) | (a & 1u)));
    }
    return f32_from(f32_norm_round_pack(0u, F32_INT_EXP - 2, a));
}

float SOFTFP_NAME(floatsisf)(int32_t a)
{
    const uint32_t sign = (uint32_t)a & F32_SIGN;
    if (!((uint32_t)a & ~F32_SIGN)) {
        return f32_from(sign ? f32_pack(F32_SIGN, F32_INT_EXP, 0u) : 0u);     /* 0 or -2^31 */
    }
    const uint32_t mag = sign ? 0u - (uint32_t)a : (uint32_t)a;
    return f32_from(f32_norm_round_pack(sign, F32_INT_EXP - 2, mag));
}

double SOFTFP_NAME(floatsidf)(int32_t a)
{
    if (!a) {
        return f64_from(0u);
    }
    const uint64_t sign = a < 0 ? F64_SIGN : 0u;
    const uint32_t mag = a < 0 ? 0u - (uint32_t)a : (uint32_t)a;
    const int shift = (int)rv_clz32(mag) + (F64_FRAC_BITS + 1 - I32_BITS);
    return f64_from(f64_pack(sign, F64_INT_EXP - 1 - shift, (uint64_t)mag << shift));
}

int32_t SOFTFP_NAME(fixsfsi)(float fa)
{
    const uint32_t a = f32_bits(fa);
    const int shift = F32_INT_EXP - f32_exp(a);
    if (shift >= I32_BITS) {
        return 0;
    }
    if (shift <= 0) {
        return (a & F32_SIGN) ? INT32_MIN : INT32_MAX;                        /* saturate (and NaN) */
    }
    const uint32_t mag = ((f32_frac(a) | F32_IMPLICIT) << (I32_BITS - 1 - F32_FRAC_BITS)) >> shift;
    return (a & F32_SIGN) ? -(int32_t)mag : (int32_t)mag;
}

uint32_t SOFTFP_NAME(fixunssfsi)(float fa)
{
    const uint32_t a = f32_bits(fa);
    const int shift = F32_INT_EXP - f32_exp(a);
    if (shift >= I32_BITS || (a & F32_SIGN)) {
        return 0u;
    }
    if (shift < 0) {
        return UINT32_MAX;
    }
    return ((f32_frac(a) | F32_IMPLICIT) << (I32_BITS - 1 - F32_FRAC_BITS)) >> shift;
}

int32_t SOFTFP_NAME(fixdfsi)(double fa)
{
    const uint64_t a = f64_bits(fa);
    const int shift = F64_INT_EXP - f64_exp(a);
    if (shift > F64_FRAC_BITS) {
        return 0;
    }
    if (shift <= F64_FRAC_BITS - (I32_BITS - 1)) {
        return (a & F64_SIGN) ? INT32_MIN : INT32_MAX;
    }
    const uint32_t mag = (uint32_t)((f64_frac(a) | F64_IMPLICIT) >> shift);
    return (a & F64_SIGN) ? -(int32_t)mag : (int32_t)mag;
}
