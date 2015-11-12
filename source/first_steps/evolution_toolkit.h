#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>

void fft_D1(double *in, double *out, int direction, parameters_t *pars);
void fft_D2(double *in, double *out, int direction, parameters_t *pars);
void set_adaptive_cutoff_fraction(size_t nc);
void fft_apply_filter(double *f, size_t N);
void mk_filter_window(double *window, size_t nmax, size_t nc);
extern double filter_window_function(double x);

#endif