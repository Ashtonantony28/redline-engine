#include "tables.h"

static int seg(const float *x, int n, float q, float *f)
{
    if (q <= x[0]) {
        *f = 0.0f;
        return 0;
    }
    if (q >= x[n - 1]) {
        *f = 1.0f;
        return n - 2;
    }
    int i = 0;
    while (i < n - 2 && q >= x[i + 1]) {
        ++i;
    }
    *f = (q - x[i]) / (x[i + 1] - x[i]);
    return i;
}

float ecu_interp1(const float *x, const float *y, int n, float xq)
{
    float f;
    int i = seg(x, n, xq, &f);
    return y[i] + f * (y[i + 1] - y[i]);
}

float ecu_interp2(const float *xr, int nr, const float *xc, int nc, const float *z, float r, float c)
{
    float fr, fc;
    int i = seg(xr, nr, r, &fr);
    int j = seg(xc, nc, c, &fc);
    const float *z0 = z + i * nc;
    const float *z1 = z0 + nc;
    float a = z0[j] + fc * (z0[j + 1] - z0[j]);
    float b = z1[j] + fc * (z1[j + 1] - z1[j]);
    return a + fr * (b - a);
}
