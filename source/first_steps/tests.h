#ifndef __tests__
#define __tests__

#include <stddef.h>

void run_all_tests(void);
void test_fft_first_derivative();
void test_fft_second_derivative();
void test_fft_apply_filter();
int equal(double a, double b);

#endif
