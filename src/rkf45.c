#include "main.h"
#if INTEGRATION_METHOD == RKF45

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "rkf45.h"
#include "toolbox.h"
#include "io.h"

static int mk_rhs_wrapper(double t, const double f[], double result[],
        void *params);

void run_rkf45()
{
    INFO(puts("Starting rkf45 integration with:"));
    INFO(printf("  Initial time: %.17f\n", pars.t.ti));
    INFO(printf("  Final time: %.17f\n", pars.t.tf));
    INFO(printf("  Initial time step dt: %.17f\n", pars.t.dt));
    INFO(printf("  Minimal time step dt: %.17f\n", MINIMAL_DELTA_T));
    INFO(printf("  Max number of steps: %zu\n", MAX_STEPS));
    // TODO: change this once I have more definitions
    INFO(printf("  Relative tolerance: %.17f\n", RELATIVE_TOLERANCE));
    INFO(printf("  Absolute tolerance: %.17f\n\n", ABSOLUTE_TOLERANCE));

    prepare_and_save_timeslice();

    gsl_odeiv2_system sys = {mk_rhs_wrapper, NULL, pars.Ntot, NULL};

    //TODO: use gsl_odeiv2_driver_alloc_scaled_new for abs err vector
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
            gsl_odeiv2_step_rkf45, pars.t.dt,
            ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    gsl_odeiv2_driver_set_hmin(d, MINIMAL_DELTA_T);

    for (size_t i = 1; i <= RKF45_OUTPUT_NUMBER; ++i) {
        double ti = i * (pars.t.tf - pars.t.ti) / RKF45_OUTPUT_NUMBER;
        int status = gsl_odeiv2_driver_apply(d, &pars.t.t, ti, field);
        if (status != GSL_SUCCESS) {
          printf ("error in GSL integration: return value=%d\n", status);
          break;
        }
        prepare_and_save_timeslice();
    }
    gsl_odeiv2_driver_free(d);
}

static int mk_rhs_wrapper(double t, const double f[], double result[],
        void *params)
{
    mk_rhs(t, (double*) f, result);
    return GSL_SUCCESS;
}

#endif
