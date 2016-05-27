#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include "main.h"

/**
 * @file toolbox.h
 *
 * @brief Struct definitions and function declarations for `toolbox.c`
 */

/**
 * @brief Flags to pass program flow information and the status of the
 * integration to the routines in `toolbox.c`.
 */
struct evolution_flags
{
    uint8_t filter; ///< If set, apply filtering (unused)
    /**
     * @brief If set, compute power spectra via `mk_power_spectrum(const
     * fftw_complex *in, struct output out)` when calling
     * `mk_gradient_squared_and_laplacian(double *in)`
     */
    uint8_t compute_pow_spec;
    /**
     * @brief If set, compute the constraints via `mk_constraints()` when
     * calling `mk_rhs(const double t, double *f, double *result)`
     */
    uint8_t compute_cstr;
};

extern struct evolution_flags evo_flags; ///< Instance of `evolution_flags`

void mk_rhs(const double t, double *f, double *result);
void mk_gradient_squared_and_laplacian(double *in);
void mk_rho_and_p(const double *f);
void mk_psi(double *f);
#ifdef ENABLE_FFT_FILTER
void apply_filter(double *inout);
#endif
void mk_summary();
void prepare_and_save_timeslice();

#endif
