#ifndef __RK4_STEPPER__
#define __RK4_STEPPER__

#include <stddef.h>
#include "main.h"

void run_RK4_stepper(parameters_t *pars);
void get_field_velocity(double* f, double *result, size_t N);
inline double potential(double f);
inline double potential_prime_term(double field_value);
double get_a_velocity(double *f, size_t nt, size_t N);
double get_total_energy(double *f, size_t N);

#endif