#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>
#include <fftw3.h>
#include "main.h"

// some flags to pass information about the progress of the stepping to the
// evolution_toolkit
typedef struct {
    uint8_t filter; // apply a filter in this step? (currently unused)
    uint8_t compute_pow_spec; // compute the power spectrum in this step?
}evolution_flags_t;

extern evolution_flags_t evo_flags;

void mk_rhs(const double t, double *f, double *result);
void mk_rho(double *f);
void mk_gradient_squared_and_laplacian(double *in);
void mk_grad_phi_times_grad_psi();
extern double potential(const double f);
extern double potential_prime(const double f);
void mk_psi_and_dpsi(double *f);
void mk_power_spectrum(const fftw_complex *in);
void apply_filter_real(double *inout);
void apply_filter_fourier(fftw_complex *inout, fftw_complex *dinout);
extern double filter_window_function(const double x);
void prepare_and_save_timeslice();
double mean(const double *in, size_t N);

#endif
