#ifndef __RK4_STEPPER__
#define __RK4_STEPPER__

#include <stddef.h>
#include "main.h"

void run_RK4_stepper(parameters_t *pars);
double mk_velocities(double *f, double a, double *result, parameters_t *pars);
extern double potential(double f);
extern double potential_prime_term(double field_value);
double mk_rho(double *f, double a, size_t N);

#endif