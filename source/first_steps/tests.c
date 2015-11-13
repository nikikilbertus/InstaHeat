#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "main.h"
#include "evolution_toolkit.h"

void run_all_tests() {
    test_fft_first_derivative();
    test_fft_second_derivative();
    test_fft_apply_filter();
}

void test_fft_first_derivative() {
    size_t N = 101;
    int count = 0;
    double *x = malloc(N * sizeof *x);
    double *f = malloc(N * sizeof *f);
    double *resn = malloc(N * sizeof *resn);
    double *resex = malloc(N * sizeof *resex);

    for (size_t i = 0; i < N; ++i)
    {
        x[i] = -PI + i * 2.0 * PI / N;
        f[i] = 3.0 * sin(x[i]) + cos(5 * x[i]);
        resex[i] = 3.0 * cos(x[i]) - 5.0 * sin(5 * x[i]);
    }

    fft_D1(f, resn, N);

    for (size_t i = 0; i < N; ++i)
    {
        count += equal(resn[i], resex[i]);
    }

    puts("test_fft_first_derivative: ");
    if (count == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
}

void test_fft_second_derivative() {
    size_t N = 101;
    int count = 0;
    double *x = malloc(N * sizeof *x);
    double *f = malloc(N * sizeof *f);
    double *resn = malloc(N * sizeof *resn);
    double *resex = malloc(N * sizeof *resex);

    for (size_t i = 0; i < N; ++i)
    {
        x[i] = -PI + i * 2.0 * PI / N;
        f[i] = 3.0 * sin(x[i]) + cos(5 * x[i]);
        resex[i] = -3.0 * sin(x[i]) - 25.0 * cos(5 * x[i]);
    }

    fft_D2(f, resn, N);

    for (size_t i = 0; i < N; ++i)
    {
        count += equal(resn[i], resex[i]);
    }

    puts("test_fft_second_derivative: ");
    if (count == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
}

void test_fft_apply_filter() {
    size_t N = 63;
    int count = 0;
    double *x = malloc(N * sizeof *x);
    double *f = malloc(N * sizeof *f);
    double *ext = malloc(N * sizeof *ext);

    for (size_t i = 0; i < N; ++i)
    {
        x[i] = -PI + i * 2.0 * PI / N;
        f[i] = sin(x[i]) + sin(16.0 * x[i]);
        ext[i] = sin(x[i]);
    }

    fft_apply_filter(f, N);

    for (size_t i = 0; i < N; ++i)
    {
        count += equal(ext[i], f[i]);
    }

    puts("test_fft_apply_filter: ");
    if (count == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
}

int equal(double a, double b) {
    return fabs(a - b) < 1e-8 ? 0 : -1;
}