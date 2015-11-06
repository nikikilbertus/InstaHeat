#ifndef __RK4_STEPPER__
#define __RK4_STEPPER__

#include <stddef.h>
#include "main.h"

void run_RK4_stepper(parameters_t *pars);
void mk_field_velocity(double *f, size_t nt, double *result, size_t N);
inline double potential(double f);
inline double potential_prime_term(double field_value);
inline double get_a_velocity(size_t nt);
inline double get_hubble(size_t nt);
double mk_rho(double *f, double a, size_t N);

#endif