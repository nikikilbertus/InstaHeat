#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>
#include "main.h"

void mk_gradient_squared_and_laplacian(double *in, double *grad2,
							double *laplacian, parameters_t *pars);
void set_adaptive_cutoff_fraction(size_t nc);
void fft_apply_filter(double *f, parameters_t *pars);
void mk_filter_window(double *out, size_t cutoffindex, size_t windowlength);
extern double filter_window_function(double x);

#endif