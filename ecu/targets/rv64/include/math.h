/* Redline ECU, portable-libm shim (M9). The slow loop's only libm calls are fabsf, expf and powf.
 * Bare-metal RV64IM has no libm, and glibc's expf / powf are not bit-reproducible elsewhere, so the
 * RV64 target and the host *reference* build (libecu_host_ref.so) both resolve <math.h> to this
 * header (it is first on the include path) and link ecu_libm.c. Same C, same IEEE-754 binary32 /
 * binary64 operations, hence bit-identical slow-loop outputs on x86-64 SSE2 and on libgcc soft-float.
 * The core sources in ecu/src are unchanged. */
#ifndef ECU_TARGET_MATH_H
#define ECU_TARGET_MATH_H

float ecu_fabsf(float x);
float ecu_expf(float x);
float ecu_powf(float x, float y);

#define fabsf ecu_fabsf
#define expf ecu_expf
#define powf ecu_powf

#endif
