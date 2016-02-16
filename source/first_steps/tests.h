#ifndef __tests__
#define __tests__

#include "main.h"

void run_all_tests();
void test_mk_gradient_squared_and_laplacian();
void test_full_multigrid();
double test_func_gradsq(const double x, const double y, const double z);
double test_func_lap(const double x, const double y, const double z);
double test_func(const double x, const double y, const double z);
double test_func_Dx(const double x, const double y, const double z);
double test_func_Dy(const double x, const double y, const double z);
double test_func_Dz(const double x, const double y, const double z);
double test_func_D2x(const double x, const double y, const double z);
double test_func_D2y(const double x, const double y, const double z);
double test_func_D2z(const double x, const double y, const double z);
void fill_field(double *f, double (*func)(const double, const double,
            const double));
int are_fields_equal(const double *f, const double *g);
int equal(const double a, const double b);

#endif
