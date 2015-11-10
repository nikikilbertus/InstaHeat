#ifndef __tests__
#define __tests__

#include <stddef.h>

void run_all_tests(void);
#ifdef USE_ACCELERATE_FRAMEWORK
void test_matrix_matrix(void);
void test_matrix_vector(void);
#endif
void test_fft_first_derivative();
void test_fft_second_derivative();
void test_fft_apply_filter();
double * create_test_matrix(size_t N, double shift, int alternate);
double * create_test_vector(size_t N, double shift, int alternate);
int equal(double a, double b);

#endif
