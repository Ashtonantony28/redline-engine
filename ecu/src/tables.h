/* Table interpolation for the slow loops (float). Clamped at the axis ends. */
#ifndef ECU_TABLES_H
#define ECU_TABLES_H

float ecu_interp1(const float *x, const float *y, int n, float xq);
float ecu_interp2(const float *xr, int nr, const float *xc, int nc, const float *z, float r, float c);

#endif
