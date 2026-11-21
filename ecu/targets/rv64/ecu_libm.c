/* Portable, deterministic fabsf / expf / powf for the ECU slow loop (see include/math.h).
 *
 * Everything is evaluated in IEEE-754 binary64 with +, -, *, / only (no fused multiply-add, no
 * rounding-mode dependence, no tables), then rounded once to binary32. On RV64IM the operations
 * are libgcc soft-float; on the host they are SSE2: the results are bit-identical by construction.
 * Accuracy: the binary64 result is good to ~1e-15 relative, so the binary32 result is the correctly
 * rounded one except within ~1e-7 ulp of a rounding boundary (checked against glibc in
 * tests/test_m9_hil.py to <= 1 ulp).
 *
 * Domain: expf for any finite x; powf for x > 0 (the slow loop raises a positive speed ratio to a
 * calibration exponent). powf(x <= 0, y) returns 0 for y > 0 and +inf otherwise, documented, not IEEE. */
#include <stdint.h>

#include "math.h"

typedef union { double f; uint64_t u; } f64_bits_t;
typedef union { float f; uint32_t u; } f32_bits_t;

#define F32_SIGN_MASK 0x80000000u
#define F64_EXP_SHIFT 52
#define F64_EXP_MASK 0x7FFu
#define F64_EXP_BIAS 1023
#define F64_FRAC_MASK 0x000FFFFFFFFFFFFFull
#define F64_INF_BITS 0x7FF0000000000000ull

static const double LN2 = 0.693147180559945309417232121458;        /* ln 2 */
static const double INV_LN2 = 1.44269504088896340735992468100;     /* 1 / ln 2 */
static const double SQRT2 = 1.41421356237309504880168872421;
/* exp argument limits: beyond these the binary32 result is +inf / 0 anyway */
static const double EXP_ARG_MAX = 89.0;       /* > ln(FLT_MAX)  = 88.72 */
static const double EXP_ARG_MIN = -104.0;     /* < ln(FLT_TRUE_MIN / 2) = -103.97 */
/* Taylor orders: |r| <= ln2/2 so r^13/13! < 2e-16; |s| <= 0.1716 so s^22/23 < 1e-18 */
#define EXP_TAYLOR_ORDER 12
#define LOG_SERIES_TERMS 10

float ecu_fabsf(float x)
{
    f32_bits_t b;
    b.f = x;
    b.u &= ~F32_SIGN_MASK;
    return b.f;
}

static double pow2i(int k)
{
    f64_bits_t b;
    b.u = (uint64_t)(k + F64_EXP_BIAS) << F64_EXP_SHIFT;
    return b.f;
}

/* e^x for EXP_ARG_MIN <= x <= EXP_ARG_MAX: x = k ln2 + r, Horner Taylor of e^r, times 2^k. */
static double exp_core(double x)
{
    double kf = x * INV_LN2;
    int k = (int)(kf >= 0.0 ? kf + 0.5 : kf - 0.5);
    double r = x - (double)k * LN2;
    double p = 1.0;
    for (int n = EXP_TAYLOR_ORDER; n >= 1; --n) {
        p = 1.0 + p * r / (double)n;
    }
    return p * pow2i(k);
}

/* ln x for finite x > 0: x = m 2^e, m in [1/sqrt2, sqrt2); ln m = 2 atanh((m-1)/(m+1)). */
static double log_core(double x)
{
    f64_bits_t b;
    b.f = x;
    int e = (int)((b.u >> F64_EXP_SHIFT) & F64_EXP_MASK) - F64_EXP_BIAS;
    b.u = (b.u & F64_FRAC_MASK) | ((uint64_t)F64_EXP_BIAS << F64_EXP_SHIFT);
    double m = b.f;
    if (m > SQRT2) {
        m *= 0.5;
        e += 1;
    }
    double s = (m - 1.0) / (m + 1.0);
    double s2 = s * s;
    double acc = 0.0;
    for (int n = LOG_SERIES_TERMS; n >= 0; --n) {
        acc = acc * s2 + 1.0 / (double)(2 * n + 1);
    }
    return 2.0 * s * acc + (double)e * LN2;
}

static float exp_clamped(double x)
{
    if (x != x) {                 /* NaN propagates (and never reaches a float -> int conversion) */
        return (float)x;
    }
    if (x > EXP_ARG_MAX) {
        f64_bits_t inf;
        inf.u = F64_INF_BITS;
        return (float)inf.f;
    }
    if (x < EXP_ARG_MIN) {
        return 0.0f;
    }
    return (float)exp_core(x);
}

float ecu_expf(float x)
{
    return exp_clamped((double)x);
}

float ecu_powf(float x, float y)
{
    if (!(x > 0.0f)) {
        return y > 0.0f ? 0.0f : exp_clamped(2.0 * EXP_ARG_MAX);
    }
    return exp_clamped((double)y * log_core((double)x));
}
