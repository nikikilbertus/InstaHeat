#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "rkf45.h"
#include "main.h"
#include "evolution_toolkit.h"

void run_rkf45()
{
    INFO(puts("Starting rkf45 integration with:"));
    INFO(printf("initial time: %f\n", pars.t.ti));
    INFO(printf("final time: %f\n", pars.t.tf));
    INFO(printf("initial time step dt: %f\n", pars.t.dt));
    INFO(printf("minimal time step dt: %f\n", MINIMAL_DELTA_T));
    INFO(printf("max number of steps: %zu\n", MAX_STEPS));
    INFO(printf("relative tolerance: %.15f\n", RELATIVE_TOLERANCE));
    INFO(printf("absolute tolerance: %.15f\n\n", ABSOLUTE_TOLERANCE));

    #ifdef OUTPUT_PHI_PS
    evo_flags.compute_pow_spec = 1;
    #endif
    #ifdef OUTPUT_CONSTRAINTS
    evo_flags.compute_cstr = 1;
    #endif
    mk_rhs(pars.t.t, field, dfield);
    evo_flags.compute_pow_spec = 0;
    evo_flags.compute_cstr = 0;
    mk_summary();
    save();

    gsl_odeiv2_system sys = {mk_rhs_wrapper, NULL, pars.Nall, NULL};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
            gsl_odeiv2_step_rkf45, pars.t.dt, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);



    for (size_t i = 1; i <= 100; i++) {
        double ti = i * pars.t.tf / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &pars.t.t, ti, field);
        if (status != GSL_SUCCESS) {
          printf ("error, return value=%d\n", status);
          break;
        }
        printf("t = %f\n", ti);
    }
    gsl_odeiv2_driver_free (d);
}

int mk_rhs_wrapper(double t, const double f[], double result[], void * params)
{
    mk_rhs(t, f, result);
    return GSL_SUCCESS;
}
