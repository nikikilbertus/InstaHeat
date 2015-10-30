#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include "tests.h"
#include "setup.h"

void run_all_tests() {
    test_matrix_vector();
    test_matrix_matrix();
}

void test_matrix_matrix() {
    size_t N = 3, count = 0;
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

    printf("\n\ntest_matrix_times_matrix: ");
    if (count == 0)
    {
        printf("passed\n\n");
    }
    else
    {
        printf("failed\n\n");
    }
}

void test_matrix_vector() {
    size_t N = 3;
    double *result = malloc(N * sizeof *result);
    double *A = create_test_matrix(N, 0.0, 0);
    double *x = create_test_vector(N, 0.0, 0);

    matrix_vector(A, x, result, N);

    printf("\n\ntest_matrix_times_vector: ");
    if (equal(result[0], 15.0) == 0 &&
        equal(result[1], 18.0) == 0 &&
        equal(result[2], 21.0) == 0)
    {
        printf("passed\n\n");
    }
    else
    {
        printf("failed\n\n");
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
    return fabs(a - b) < 1e-10 ? 0 : -1;
}