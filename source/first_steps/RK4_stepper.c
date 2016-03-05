#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

void run_rk4() {
    size_t Ntot = pars.Ntot;
    size_t Nall = pars.Nall;
    size_t Nt = pars.t.Nt;
    double dt = pars.t.dt;
    double t = pars.t.t;

    double *k1, *k2, *k3, *k4, *tmp_k;
    k1    = fftw_malloc(Ntot * sizeof *k1);
    k2    = fftw_malloc(Ntot * sizeof *k2);
    k3    = fftw_malloc(Ntot * sizeof *k3);
    k4    = fftw_malloc(Ntot * sizeof *k4);
    tmp_k = fftw_malloc(Nall * sizeof *tmp_k);

    if (!(k1 && k2 && k3 && k4 && tmp_k))
    {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }

    RUNTIME_INFO(puts("Starting RK4 time evolution with:"));
    RUNTIME_INFO(printf("initial time: %f\n", pars.t.ti));
    RUNTIME_INFO(printf("final time: %f\n", pars.t.tf));
    RUNTIME_INFO(printf("time step dt: %f\n", dt));
    RUNTIME_INFO(printf("number of steps: %zu\n", Nt));

    #ifdef SHOW_TIMING_INFO
    double secs = -get_wall_time();
    #endif

    for (size_t nt = 0; t < pars.t.tf; ++nt)
    {
        #ifdef ENABLE_FFT_FILTER
        apply_filter_real(field);
        #endif

        // to precisely reach final time in the last step, change dt
        if (t + dt * 1.0001 > pars.t.tf)
        {
            dt = pars.t.tf - t;
            pars.t.dt = dt;
            RUNTIME_INFO(printf("overshoot, new dt = %f\n", dt));
        }

        // step 1 (and write out data if required)
        if (nt % pars.file.skip == 0)
        {
            evo_flags.compute_pow_spec = 1;
            mk_rhs(t, field, k1);
            evo_flags.compute_pow_spec = 0;
            mk_means_and_variances();
            save();
        }
        else
        {
            mk_rhs(t, field, k1);
        }

        // step 2
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k1[i] / 2.0;
        }
        mk_rhs(t + dt / 2.0, tmp_k, k2);

        // step 3
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k2[i] / 2.0;
        }
        mk_rhs(t + dt / 2.0, tmp_k, k3);

        // step 4
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k3[i];
        }
        mk_rhs(t + dt, tmp_k, k4);

        // perform time step
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            field[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }

        t += dt;
        pars.t.t = t;
    }

    // make sure to write out last time slice
    prepare_and_save_timeslice();

    // info about last timeslice
    if (fabs(pars.t.tf - pars.t.t) > 1e-10)
    {
        RUNTIME_INFO(puts("The time of the last step does not coincide "
                          "with the specified final time."));
    }

    RUNTIME_INFO(puts("Finished rk4"));
    #ifdef SHOW_TIMING_INFO
    secs += get_wall_time();
    RUNTIME_INFO(printf("time: %f seconds\n\n", secs));
    RUNTIME_INFO(puts("Writing simulation meta data to disk\n"));
    double val[1] = {secs};
    h5_write_parameter(H5_RUNTIME_STEPPER_NAME, val, 1);
    #endif

    fftw_free(k1);
    fftw_free(k2);
    fftw_free(k3);
    fftw_free(k4);
    fftw_free(tmp_k);
}
