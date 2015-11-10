#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "main.h"
#include "evolution_toolkit.h"

void run_all_tests() {
#ifdef USE_ACCELERATE_FRAMEWORK
    test_matrix_vector();
    test_matrix_matrix();
#endif
    test_fft_first_derivative();
    test_fft_second_derivative();
    test_fft_apply_filter();
}
#ifdef USE_ACCELERATE_FRAMEWORK
void test_matrix_matrix() {
    size_t N = 3;
    int count = 0;
    double *result = malloc(N*N * sizeof *result);
    double *A = create_test_matrix(N, 0.0, 0);
    double *B = create_test_matrix(N, 1.0, 0);

    double actual_result[] = {24.0, 30.0, 36.0, 51.0, 66.0,
                                81.0, 78.0, 102.0, 126.0};

    matrix_matrix(A, B, result, N);

    for (size_t i = 0; i < N*N; ++i)
    {
        count += equal(result[i], actual_result[i]);
    }

    puts("test_matrix_times_matrix: ");
    if (count == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
}

void test_matrix_vector() {
    size_t N = 3;
    double *result = malloc(N * sizeof *result);
    double *A = create_test_matrix(N, 0.0, 0);
    double *x = create_test_vector(N, 0.0, 0);

    matrix_vector(A, x, result, N);

    puts("test_matrix_times_vector: ");
    if (equal(result[0], 15.0) == 0 &&
        equal(result[1], 18.0) == 0 &&
        equal(result[2], 21.0) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
}
#endif

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

double * create_test_matrix(size_t N, double shift, int alternate) {
    size_t N2 = N*N;
    double *matrix = malloc(N2 * sizeof *matrix);
    for (size_t i = 0; i < N2; ++i)
    {
        matrix[i] = pow(-1.0, alternate * i) * i + shift;
    }
    return matrix;
}

double * create_test_vector(size_t N, double shift, int alternate) {
    double *vector = malloc(N * sizeof *vector);
    for (size_t i = 0; i < N; ++i)
    {
        vector[i] = pow(-1.0, alternate * i) * i + shift;
    }
    return vector;
}

int equal(double a, double b) {
    return fabs(a - b) < 1e-8 ? 0 : -1;
}