/* int8 MLP inference with Xmlacc Tier-1 instructions ($VORTEX_V_ROOT/docs/XMLACC-SPEC-v1.md),
 * the same mapping as Vortex-V's NeuroFlight workload (verif/workloads/neuroflight):
 *     acc = sat_int32(b_q[j] + sum_i W_q[j,i] x_q[i])      mlacc.macc, 8 int8 lanes per instruction
 *     hidden:  y = clamp(rshift_round(acc * M[j], S[j]), -128, 127)      mul + mlacc.rqs
 *              a = max(y, 0)                                             ReLU
 *     last:    acc is the output
 * Integer only (this unit is compiled with the no-float check). */
#include "vortex.h"             /* xm_macc: .insn encodings maintained by Vortex-V */

#define KNOCK_DEFINE_TABLES
#include "knock_kernel.h"

#define XM_LANES 8u             /* int8 lanes per 64-bit register */
#define XM_OPCODE_CUSTOM0 0x0b
#define XM_FUNCT3_RQS 5         /* mlacc.rqs, funct7 = 0 (XMLACC-SPEC-v1 table) */

typedef uint64_t __attribute__((may_alias)) lanes_t;

static inline int64_t xm_rqs(int64_t v, uint64_t sh)
{
    int64_t r;
    __asm__ volatile(".insn r %c3, %c4, 0, %0, %1, %2"
                     : "=r"(r) : "r"(v), "r"(sh), "i"(XM_OPCODE_CUSTOM0), "i"(XM_FUNCT3_RQS));
    return r;
}

static inline int64_t sat_int32(int64_t acc)
{
    if (acc > INT32_MAX) {
        return INT32_MAX;
    }
    if (acc < INT32_MIN) {
        return INT32_MIN;
    }
    return acc;
}

int32_t knock_infer(const int8_t *x)
{
    static _Alignas(8) int8_t act[2][KNOCK_MAX_HIDDEN];
    const int8_t *in = x;
    for (uint32_t l = 0u; l < KNOCK_N_LAYERS; ++l) {
        const knock_layer_t *ly = &KNOCK_LAYERS[l];
        const lanes_t *w = (const lanes_t *)ly->w;
        const lanes_t *xin = (const lanes_t *)in;
        int8_t *out = act[l & 1u];
        for (uint32_t j = 0u; j < ly->n_out; ++j) {
            int64_t acc = ly->b[j];
            for (uint32_t k = 0u; k < ly->in_words; ++k) {
                acc = xm_macc(xin[k], *w++, acc);
            }
            acc = sat_int32(acc);
            if (ly->m == 0) {
                return (int32_t)acc;                 /* output layer: one unit */
            }
            int64_t y = xm_rqs(acc * ly->m[j], ly->s[j]);
            out[j] = (int8_t)(y < 0 ? 0 : y);
        }
        in = out;
    }
    return 0;
}
