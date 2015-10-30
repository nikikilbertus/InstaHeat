#ifndef __tests__
#define __tests__

#include <stddef.h>

void run_all_tests(void);
void test_matrix_matrix(void);
void test_matrix_vector(void);
double * create_test_matrix(size_t N, double shift, int alternate);
double * create_test_vector(size_t N, double shift, int alternate);
int equal(double a, double b);

#endif
