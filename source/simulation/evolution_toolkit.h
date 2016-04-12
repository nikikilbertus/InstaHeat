#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>
#include <fftw3.h>
#include "main.h"

// some flags to pass information about the progress of the stepping to the
// evolution_toolkit
struct evolution_flags
{
    uint8_t filter; // apply a filter in this step? (currently unused)
    uint8_t compute_pow_spec; // compute the power spectrum in this step?
};

extern struct evolution_flags evo_flags;

void mk_rhs(const double t, double *f, double *result);
void mk_gradient_squared_and_laplacian(double *in);
void assemble_gradient_squared();
void mk_rho(const double *f);
extern double potential(const double f);
extern double potential_prime(const double f);
void mk_psi(double *f);
void mk_power_spectrum(const fftw_complex *in, struct output out);
void apply_filter_real(double *inout);
void apply_filter_fourier(fftw_complex *phi_io, fftw_complex *dphi_io,
        fftw_complex *psi_io, fftw_complex *dpsi_io);
void prepare_and_save_timeslice();
void center(double *f, const size_t N);
void mk_means_and_variances();
extern double mean(const double *f, const size_t N);
void mean_var_min_max(const double *f, double *smry);
double variance(const double mean, const double *f, const size_t N);
void contains_nan(const double *f, const size_t N);
void contains_nanc(const complex *f, const size_t N);

#endif
