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
     * @brief Set if output at current timeslice is desired. This tells methods
     * in toolbox to compute additional quantities.
     */
    uint8_t output;
};

extern struct evolution_flags evo_flags; ///< Instance of `evolution_flags`

void mk_rhs(const double t, double *f, double *result);
void mk_gradient_squared_and_laplacian(double *in);
void mk_rho(const double *f);
void mk_psi(double *f);
void mk_summary();
void prepare_and_save_timeslice();

#endif
