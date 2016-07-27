#include "main.h"
#ifdef GSL_STEPPER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "gsl_stepper.h"
#include "toolbox.h"
#include "io.h"

/**
 * @file gsl_stepper.c
 *
 * @brief Contains various explicit ODE Runge Kutta stepper as implemented by
 * GNU Scientific Library (GSL).
 *
 * Available methods are: Runge Kutta Felberg 4(5) (RKF45), Runge Kutta
 * Cash-Karp 4(5) (RKCK), Runge Kutta Dormand Prince 8(9) (DOPRI89).  Only the
 * function run_gsl_stepper() is called from outside this file. It is the only
 * one that needs to be visible.
 *
 * @see <a href="https://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html">GSL ODE Manual</a>
 */

static void write_start_info();
static int mk_rhs_wrap(double t, const double f[], double res[], void *params);

/**
 * @brief Run one of the supported GSL integration routines to evolve the
 * fields.
 *
 * This is the main routine of `gsl_stepper.c`. It integrates the field
 * (assuming initial values are already present) and writes the evolution to
 * disk as specified by calling integration routines implemented by GNU GSL.
 */
void run_gsl_stepper()
{
    write_start_info();
    prepare_and_save_timeslice();
    gsl_odeiv2_system sys = {mk_rhs_wrap, NULL, pars.Ntot, NULL};
    gsl_odeiv2_driver *d;

    //TODO: use gsl_odeiv2_driver_alloc_scaled_new for abs err vector
    #if INTEGRATION_METHOD == RKF45
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, pars.t.dt,
            ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    #elif INTEGRATION_METHOD == RKCK
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkck, pars.t.dt,
            ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    #elif INTEGRATION_METHOD == DOPRI89
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, pars.t.dt,
            ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    #endif
    gsl_odeiv2_driver_set_hmin(d, MINIMAL_DELTA_T);

    for (size_t i = 1; i <= GSL_OUTPUT_NUMBER; ++i) {
        double ti = i * (pars.t.tf - pars.t.ti) / GSL_OUTPUT_NUMBER;
        int stat = gsl_odeiv2_driver_apply(d, &pars.t.t, ti, field);
        if (stat != GSL_SUCCESS) {
          printf("error in GSL integration: %d\n", stat);
          break;
        }
        prepare_and_save_timeslice();
    }
    gsl_odeiv2_driver_free(d);
}

/**
 * @brief A wrapper for `mk_rhs(const double t, double *f, double *result)` in
 * `toolbox.c` to suit the GSL integration routines.
 *
 * @param[in] t The current time.
 * @param[in] f The fields.
 * @param[out] res The right hand side of the partial differential equation,
 * i.e. the first temporal derivatives of the fields in @p f.
 * @param[in] params Arbitrary parameters.
 *
 * @see `mk_rhs(const double t, double *f, double *result)` in `toolbox.c`
 */
static int mk_rhs_wrap(double t, const double f[], double res[], void *params)
{
    mk_rhs(t, (double*) f, res);
    return GSL_SUCCESS;
}

/**
 * @brief Print information about the integration to stdout.
 */
static void write_start_info()
{
    #if INTEGRATION_METHOD == RKF45
    INFO(puts("Starting rkf45 integration with:"));
    #elif INTEGRATION_METHOD == RKCK
    INFO(puts("Starting rkck integration with:"));
    #elif INTEGRATION_METHOD == DOPRI89
    INFO(puts("Starting dopri89 integration with:"));
    #endif
    INFO(printf("  Initial time: %g\n", pars.t.ti));
    INFO(printf("  Final time: %g\n", pars.t.tf));
    INFO(printf("  Initial time step dt: %g\n", pars.t.dt));
    INFO(printf("  Minimal time step dt: %g\n", MINIMAL_DELTA_T));
    INFO(printf("  Max number of steps: %g\n", MAX_STEPS));
    // TODO: change this once I have more definitions
    INFO(printf("  Relative tolerance: %g\n", RELATIVE_TOLERANCE));
    INFO(printf("  Absolute tolerance: %g\n\n", ABSOLUTE_TOLERANCE));
}
#endif
