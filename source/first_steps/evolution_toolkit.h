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

void mk_gradient_squared_and_laplacian(double *in, double *grad2,
							double *laplacian, parameters_t *pars);
void fft_apply_filter(fftw_complex *inout, parameters_t *pars);
double mk_velocities(double t, double *f, double a, double *result, parameters_t *pars);
extern double potential(double f);
extern double potential_prime(double f);
double mk_rho(double *f, double a, parameters_t *pars);
void mk_power_spectrum(fftw_complex *in, parameters_t *pars);
void mk_filter_window(double *out, size_t cutoffindex, size_t windowlength);
extern double filter_window_function(double x);

#endif