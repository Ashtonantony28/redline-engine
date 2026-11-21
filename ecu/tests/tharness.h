/* Tiny self-contained C test harness: no external dependencies. */
#ifndef THARNESS_H
#define THARNESS_H

#include <math.h>
#include <stdio.h>

static int th_fail, th_checks, th_tests;

#define CHECK(cond)                                                                       \
    do {                                                                                  \
        ++th_checks;                                                                      \
        if (!(cond)) {                                                                    \
            ++th_fail;                                                                    \
            printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);                      \
        }                                                                                 \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                             \
    do {                                                                                  \
        ++th_checks;                                                                      \
        double th_a = (double)(a), th_b = (double)(b);                                    \
        if (!(fabs(th_a - th_b) <= (tol))) {                                              \
            ++th_fail;                                                                    \
            printf("  FAIL %s:%d: %s = %.6g, expected %.6g +/- %.3g\n", __FILE__, __LINE__, \
                   #a, th_a, th_b, (double)(tol));                                        \
        }                                                                                 \
    } while (0)

#define RUN(fn)                                                                           \
    do {                                                                                  \
        int th_before = th_fail;                                                          \
        ++th_tests;                                                                       \
        fn();                                                                             \
        printf("%s %s\n", th_fail == th_before ? "ok  " : "FAIL", #fn);                   \
    } while (0)

#define TH_SUMMARY()                                                                      \
    (printf("%d tests, %d checks, %d failures\n", th_tests, th_checks, th_fail), th_fail ? 1 : 0)

#endif
