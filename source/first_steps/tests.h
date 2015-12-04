#ifndef __tests__
#define __tests__

#include "main.h"

void run_all_tests(parameters_t *pars);
void test_mk_gradient_squared_and_laplacian(parameters_t *pars);
void test_fft_apply_filter(parameters_t *pars);
double test_func(double x, double y, double z);
double test_func_gradsq(double x, double y, double z);
double test_func_lap(double x, double y, double z);
double test_func_Dx(double x, double y, double z);
double test_func_Dy(double x, double y, double z);
double test_func_Dz(double x, double y, double z);
double test_func_D2(double x, double y, double z);
void fill_field(double *f, double (*func)(double, double, double),
                parameters_t *pars);
int are_fields_equal(double *f, double *g, parameters_t *pars);
int equal(double a, double b);

#endif
