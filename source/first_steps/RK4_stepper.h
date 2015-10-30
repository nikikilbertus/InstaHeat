#ifndef __RK4_STEPPER__
#define __RK4_STEPPER__

#include <stddef.h>

void run_RK4_stepper(double dt, size_t N);
void get_field_velocity(double* f, double *df, double *result, size_t N);

#endif