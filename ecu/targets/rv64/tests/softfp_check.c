/* Host check of the RV64 target's arithmetic helpers (M9.1): rv_softfp.c and rv_div.c are
 * compiled natively under other names and every operation is compared with the host's own
 * IEEE-754 arithmetic (SSE2, correctly rounded) and integer division, bit for bit:
 *   - all pairs of a list of special values (zeros, subnormals, powers of two and their
 *     neighbours, the largest finite values, infinities, NaNs),
 *   - seeded random operands: uniform bit patterns, and pairs with nearby exponents (the cases
 *     where add / subtract cancel and round).
 * NaN results are compared as "is a NaN" only (x86's default NaN differs in the sign bit).
 * Float -> integer conversions are compared in range; out of range the library saturates.
 * Prints one line per operation and exits non-zero on any mismatch.
 *   usage: softfp_check [random vectors per operation]   (built and run by `make softfp-check`) */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define SOFTFP_NAME(x) t_##x
#define RV_DIV_NAME(x) t_##x
#include "../rv_softfp.c"
#include "../rv_div.c"

#define DEFAULT_RANDOM 2000000u
#define SEED 0x9E3779B97F4A7C15ull
#define NEAR_EXP_SPREAD 3u       /* exponent difference range of the "nearby" random pairs */
#define FIX32_LIMIT 2147483520.0f           /* largest float below 2^31 */
#define FIXU32_LIMIT 4294967040.0f          /* largest float below 2^32 */
#define FIX64_LIMIT 2147483647.0

static uint64_t g_rng = SEED;
static uint64_t rnd(void)        /* splitmix64 */
{
    uint64_t z = (g_rng += SEED);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

static const uint32_t S32[] = {
    0x00000000u, 0x80000000u, 0x00000001u, 0x80000001u, 0x00000002u, 0x007FFFFFu, 0x807FFFFFu, 0x00400000u,
    0x00800000u, 0x80800000u, 0x00800001u, 0x00FFFFFFu, 0x01000000u, 0x3F7FFFFFu, 0x3F800000u, 0xBF800000u,
    0x3F800001u, 0x3FFFFFFFu, 0x40000000u, 0x40400000u, 0x3F000000u, 0x3EFFFFFFu, 0x4B800000u, 0x4B7FFFFFu,
    0x4F000000u, 0xCF000000u, 0x4EFFFFFFu, 0x4F800000u, 0x4F7FFFFFu, 0x5F000000u, 0x7F7FFFFFu, 0xFF7FFFFFu,
    0x7F000000u, 0x7F800000u, 0xFF800000u, 0x7FC00000u, 0xFFC00000u, 0x7F800001u, 0x7FA55555u, 0x33800000u,
    0x34000000u, 0x3A83126Fu, 0x41200000u, 0x3DCCCCCDu, 0x0D800000u, 0x71800000u, 0x00FFFFFEu, 0x3FC00000u};
static const uint64_t S64[] = {
    0x0000000000000000ull, 0x8000000000000000ull, 0x0000000000000001ull, 0x8000000000000001ull,
    0x000FFFFFFFFFFFFFull, 0x800FFFFFFFFFFFFFull, 0x0008000000000000ull, 0x0010000000000000ull,
    0x8010000000000000ull, 0x0010000000000001ull, 0x001FFFFFFFFFFFFFull, 0x0020000000000000ull,
    0x3FEFFFFFFFFFFFFFull, 0x3FF0000000000000ull, 0xBFF0000000000000ull, 0x3FF0000000000001ull,
    0x3FFFFFFFFFFFFFFFull, 0x4000000000000000ull, 0x4008000000000000ull, 0x3FE0000000000000ull,
    0x3FDFFFFFFFFFFFFFull, 0x4330000000000000ull, 0x432FFFFFFFFFFFFFull, 0x41E0000000000000ull,
    0xC1E0000000000000ull, 0x41DFFFFFFFC00000ull, 0x7FEFFFFFFFFFFFFFull, 0xFFEFFFFFFFFFFFFFull,
    0x7FE0000000000000ull, 0x7FF0000000000000ull, 0xFFF0000000000000ull, 0x7FF8000000000000ull,
    0xFFF8000000000000ull, 0x7FF0000000000001ull, 0x3FE62E42FEFA39EFull, 0x3FF71547652B82FEull,
    0x3CA0000000000000ull, 0x3CB0000000000000ull, 0x36A0000000000000ull, 0x47EFFFFFE0000000ull,
    0x47EFFFFFF0000000ull, 0x380FFFFFE0000000ull, 0x36A0000000000001ull, 0x3690000000000000ull,
    0x3FB999999999999Aull, 0x4024000000000000ull, 0x1FF0000000000000ull, 0x5FF0000000000000ull};
#define N32 (sizeof(S32) / sizeof(S32[0]))
#define N64 (sizeof(S64) / sizeof(S64[0]))

static int g_fail;
static uint64_t g_bad, g_n;

static void begin(void) { g_bad = 0; g_n = 0; }
static void end(const char *name)
{
    printf("%-14s %10" PRIu64 " vectors, %" PRIu64 " mismatches\n", name, g_n, g_bad);
    g_fail |= g_bad != 0;
}

static void same32(float got, float want, uint32_t a, uint32_t b, const char *op)
{
    const uint32_t g = f32_bits(got), w = f32_bits(want);
    g_n++;
    if (g != w && !(f32_is_nan(g) && f32_is_nan(w))) {
        if (g_bad++ < 5) {
            printf("  %s(%08x, %08x) = %08x, host %08x\n", op, a, b, g, w);
        }
    }
}

static void same64(double got, double want, uint64_t a, uint64_t b, const char *op)
{
    const uint64_t g = f64_bits(got), w = f64_bits(want);
    g_n++;
    if (g != w && !(f64_is_nan(g) && f64_is_nan(w))) {
        if (g_bad++ < 5) {
            printf("  %s(%016" PRIx64 ", %016" PRIx64 ") = %016" PRIx64 ", host %016" PRIx64 "\n", op, a, b, g, w);
        }
    }
}

static void samei(int64_t got, int64_t want, uint64_t a, uint64_t b, const char *op)
{
    g_n++;
    if (got != want) {
        if (g_bad++ < 5) {
            printf("  %s(%016" PRIx64 ", %016" PRIx64 ") = %" PRId64 ", host %" PRId64 "\n", op, a, b, got, want);
        }
    }
}

/* random pair k: even = independent bit patterns, odd = b's exponent near a's */
static void pair32(uint64_t k, uint32_t *a, uint32_t *b)
{
    const uint64_t r = rnd();
    *a = (uint32_t)r;
    *b = (uint32_t)(r >> 32);
    if (k & 1u) {
        const uint32_t e = (uint32_t)f32_exp(*a) + (uint32_t)(rnd() % (2u * NEAR_EXP_SPREAD + 1u)) - NEAR_EXP_SPREAD;
        *b = (*b & ~F32_INF) | ((e & F32_EXP_MAX) << F32_FRAC_BITS);
    }
}

static void pair64(uint64_t k, uint64_t *a, uint64_t *b)
{
    *a = rnd();
    *b = rnd();
    if (k & 1u) {
        const uint64_t e = (uint64_t)f64_exp(*a) + rnd() % (2u * NEAR_EXP_SPREAD + 1u) - NEAR_EXP_SPREAD;
        *b = (*b & ~F64_INF) | ((e & F64_EXP_MAX) << F64_FRAC_BITS);
    }
}

typedef float (*op32_t)(float, float);
typedef double (*op64_t)(double, double);
static float h_add32(float a, float b) { return a + b; }
static float h_sub32(float a, float b) { return a - b; }
static float h_mul32(float a, float b) { return a * b; }
static float h_div32(float a, float b) { return a / b; }
static double h_add64(double a, double b) { return a + b; }
static double h_sub64(double a, double b) { return a - b; }
static double h_mul64(double a, double b) { return a * b; }
static double h_div64(double a, double b) { return a / b; }

static void bin32(const char *name, op32_t soft, op32_t host, uint64_t n)
{
    begin();
    for (size_t i = 0; i < N32; ++i) {
        for (size_t j = 0; j < N32; ++j) {
            same32(soft(f32_from(S32[i]), f32_from(S32[j])), host(f32_from(S32[i]), f32_from(S32[j])), S32[i], S32[j], name);
        }
    }
    for (uint64_t k = 0; k < n; ++k) {
        uint32_t a, b;
        pair32(k, &a, &b);
        same32(soft(f32_from(a), f32_from(b)), host(f32_from(a), f32_from(b)), a, b, name);
    }
    end(name);
}

static void bin64(const char *name, op64_t soft, op64_t host, uint64_t n)
{
    begin();
    for (size_t i = 0; i < N64; ++i) {
        for (size_t j = 0; j < N64; ++j) {
            same64(soft(f64_from(S64[i]), f64_from(S64[j])), host(f64_from(S64[i]), f64_from(S64[j])), S64[i], S64[j], name);
        }
    }
    for (uint64_t k = 0; k < n; ++k) {
        uint64_t a, b;
        pair64(k, &a, &b);
        same64(soft(f64_from(a), f64_from(b)), host(f64_from(a), f64_from(b)), a, b, name);
    }
    end(name);
}

static int sgn(int v) { return (v > 0) - (v < 0); }

static void cmp32_one(uint32_t ua, uint32_t ub)
{
    const float a = f32_from(ua), b = f32_from(ub);
    samei(t_eqsf2(a, b) == 0, a == b, ua, ub, "eqsf2");
    samei(t_nesf2(a, b) != 0, a != b, ua, ub, "nesf2");
    samei(t_ltsf2(a, b) < 0, a < b, ua, ub, "ltsf2");
    samei(t_lesf2(a, b) <= 0, a <= b, ua, ub, "lesf2");
    samei(t_gtsf2(a, b) > 0, a > b, ua, ub, "gtsf2");
    samei(t_gesf2(a, b) >= 0, a >= b, ua, ub, "gesf2");
    samei(t_unordsf2(a, b) != 0, a != a || b != b, ua, ub, "unordsf2");
    if (a == a && b == b) {
        samei(sgn(t_ltsf2(a, b)), (a > b) - (a < b), ua, ub, "cmpsf2 sign");
    }
}

static void cmp64_one(uint64_t ua, uint64_t ub)
{
    const double a = f64_from(ua), b = f64_from(ub);
    samei(t_eqdf2(a, b) == 0, a == b, ua, ub, "eqdf2");
    samei(t_nedf2(a, b) != 0, a != b, ua, ub, "nedf2");
    samei(t_ltdf2(a, b) < 0, a < b, ua, ub, "ltdf2");
    samei(t_ledf2(a, b) <= 0, a <= b, ua, ub, "ledf2");
    samei(t_gtdf2(a, b) > 0, a > b, ua, ub, "gtdf2");
    samei(t_gedf2(a, b) >= 0, a >= b, ua, ub, "gedf2");
    samei(t_unorddf2(a, b) != 0, a != a || b != b, ua, ub, "unorddf2");
}

static void conv32_one(uint32_t ua)
{
    const float a = f32_from(ua);
    same64(t_extendsfdf2(a), (double)a, ua, 0, "extendsfdf2");
    if (a >= -FIX32_LIMIT && a <= FIX32_LIMIT) {
        samei(t_fixsfsi(a), (int32_t)a, ua, 0, "fixsfsi");
    } else if (a == a) {
        samei(t_fixsfsi(a), a < 0.0f ? INT32_MIN : INT32_MAX, ua, 0, "fixsfsi saturate");
    }
    if (a > -1.0f && a <= FIXU32_LIMIT) {
        samei(t_fixunssfsi(a), (uint32_t)a, ua, 0, "fixunssfsi");
    } else if (a == a) {
        samei(t_fixunssfsi(a), a < 0.0f ? 0 : UINT32_MAX, ua, 0, "fixunssfsi saturate");
    }
    same32(t_floatsisf((int32_t)ua), (float)(int32_t)ua, ua, 0, "floatsisf");
    same32(t_floatunsisf(ua), (float)ua, ua, 0, "floatunsisf");
    same64(t_floatsidf((int32_t)ua), (double)(int32_t)ua, ua, 0, "floatsidf");
}

static void conv64_one(uint64_t ua)
{
    const double a = f64_from(ua);
    same32(t_truncdfsf2(a), (float)a, (uint32_t)(ua >> 32), (uint32_t)ua, "truncdfsf2");
    if (a > -FIX64_LIMIT - 1.0 && a < FIX64_LIMIT + 1.0) {
        samei(t_fixdfsi(a), (int32_t)a, ua, 0, "fixdfsi");
    } else if (a == a) {
        samei(t_fixdfsi(a), a < 0.0 ? INT32_MIN : INT32_MAX, ua, 0, "fixdfsi saturate");
    }
}

static void divs(uint64_t n)
{
    begin();
    static const uint64_t D[] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 255, 256, 384, 36571, 98304, 0x7FFFFFFFull, 0x80000000ull,
                                 0xFFFFFFFFull, 0x100000000ull, 0x7FFFFFFFFFFFFFFFull, 0x8000000000000000ull,
                                 0x8000000000000001ull, 0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull};
    const size_t nd = sizeof(D) / sizeof(D[0]);
    for (uint64_t k = 0; k < n + nd * nd; ++k) {
        uint64_t a, b;
        if (k < nd * nd) {
            a = D[k / nd];
            b = D[k % nd];
        } else {
            a = rnd() >> (rnd() % 64u);          /* every magnitude of dividend and divisor */
            b = rnd() >> (rnd() % 64u);
            if (k % 3u == 0u) {
                b = D[k % nd];                   /* repeated divisors exercise the reciprocal cache */
            }
        }
        samei((int64_t)t_udivdi3(a, b), (int64_t)(b ? a / b : ~(uint64_t)0), a, b, "udivdi3");
        int64_t want;
        if (b == 0u) {
            want = -1;
        } else if ((int64_t)a == INT64_MIN && (int64_t)b == -1) {
            want = INT64_MIN;
        } else {
            want = (int64_t)a / (int64_t)b;
        }
        samei(t_divdi3((int64_t)a, (int64_t)b), want, a, b, "divdi3");
    }
    end("udivdi3/divdi3");
}

int main(int argc, char **argv)
{
    const uint64_t n = argc > 1 ? strtoull(argv[1], NULL, 0) : DEFAULT_RANDOM;
    bin32("addsf3", t_addsf3, h_add32, n);
    bin32("subsf3", t_subsf3, h_sub32, n);
    bin32("mulsf3", t_mulsf3, h_mul32, n);
    bin32("divsf3", t_divsf3, h_div32, n);
    bin64("adddf3", t_adddf3, h_add64, n);
    bin64("subdf3", t_subdf3, h_sub64, n);
    bin64("muldf3", t_muldf3, h_mul64, n);
    bin64("divdf3", t_divdf3, h_div64, n);
    begin();
    for (size_t i = 0; i < N32; ++i) {
        for (size_t j = 0; j < N32; ++j) {
            cmp32_one(S32[i], S32[j]);
        }
    }
    for (size_t i = 0; i < N64; ++i) {
        for (size_t j = 0; j < N64; ++j) {
            cmp64_one(S64[i], S64[j]);
        }
    }
    for (uint64_t k = 0; k < n; ++k) {
        uint32_t a, b;
        uint64_t c, d;
        pair32(k, &a, &b);
        cmp32_one(a, (k % 7u) ? b : a);
        pair64(k, &c, &d);
        cmp64_one(c, (k % 7u) ? d : c);
    }
    end("comparisons");
    begin();
    for (size_t i = 0; i < N32; ++i) {
        conv32_one(S32[i]);
    }
    for (size_t i = 0; i < N64; ++i) {
        conv64_one(S64[i]);
    }
    for (uint64_t k = 0; k < n; ++k) {
        const uint64_t r = rnd();
        conv32_one((uint32_t)r);
        conv32_one((uint32_t)(r >> 32) >> (r % 32u));        /* small integers too */
        conv64_one(rnd());
        /* doubles in the binary32 range, where truncdfsf2 rounds instead of overflowing */
        conv64_one(f64_bits((double)f32_from((uint32_t)r)) ^ (rnd() & 0x3FFFFFFFull));
    }
    end("conversions");
    divs(n);
    printf(g_fail ? "FAIL\n" : "PASS: bit-identical to the host FPU and integer divide\n");
    return g_fail;
}
