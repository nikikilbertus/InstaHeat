#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>
#include <fftw3.h>
#include "main.h"

typedef struct {
    uint8_t filter;
    uint8_t compute_pow_spec;
}evolution_flags_t;

extern evolution_flags_t evo_flags;

void mk_gradient_squared_and_laplacian(double *in, double *grad2, double *lap);
void apply_filter_real(double *inout);
void apply_filter_fourier(fftw_complex *inout);
void mk_velocities(double t, double *f, double *result);
extern double potential(double f);
extern double potential_prime(double f);
double mk_rho(double *f);
void mk_power_spectrum(fftw_complex *in);
extern double filter_window_function(double x);

#endif
