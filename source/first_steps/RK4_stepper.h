#ifndef __RK4_STEPPER__
#define __RK4_STEPPER__

#include <stddef.h>
#include "main.h"

void run_RK4_stepper(parameters_t *pars);
void get_field_velocity(double* f, double *result, size_t N);
void spectral_op_D2(double *f, double *result, size_t N);
void fft_D2(double *f, double *result, size_t N);

#endif