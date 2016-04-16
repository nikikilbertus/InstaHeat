#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

void runrkf45()
{
    double mu = 10;
    gsl_odeiv2_system sys = {mk_rhs_wrapper, NULL, pars.Nall, NULL};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                    ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE, 0.0);
    for (i = 1; i <= 100; i++) {
        double ti = i * pars.t.tf / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &pars.t.t, ti, field);
        if (status != GSL_SUCCESS) {
          printf ("error, return value=%d\n", status);
          break;
        }
    }
    gsl_odeiv2_driver_free (d);
    return 0;
}

int mk_rhs_wrapper(double t, const double f[], double result[], void * params)
{
    mk_rhs(t, f, result);
    return GSL_SUCCESS;
}
