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
    size_t N = pars.N;
    size_t N2 = 2 * N;
    size_t Ntot = N2 + 1;
    size_t Nt = pars.t.Nt;
    double dt = pars.t.dt;
    double t = pars.t.t;

    double *k1, *k2, *k3, *k4, *tmp_k;
    k1    = fftw_malloc(Ntot * sizeof *k1);
    k2    = fftw_malloc(Ntot * sizeof *k2);
    k3    = fftw_malloc(Ntot * sizeof *k3);
    k4    = fftw_malloc(Ntot * sizeof *k4);
    tmp_k = fftw_malloc(Ntot * sizeof *tmp_k);

    if (!(k1 && k2 && k3 && k4 && tmp_k))
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    RUNTIME_INFO(puts("Starting RK4 time evolution with:"));
    RUNTIME_INFO(printf("initial time: %f\n", pars.t.ti));
    RUNTIME_INFO(printf("final time: %f\n", pars.t.tf));
    RUNTIME_INFO(printf("time step dt: %f\n", dt));
    RUNTIME_INFO(printf("number of steps: %zu\n", Nt));
    RUNTIME_INFO(puts("Using DFT (fftw3) for spatial derivatives."));

    #ifdef ENABLE_FFT_FILTER
        RUNTIME_INFO(puts("Frequency cutoff filtering enabled."));
    #else
        RUNTIME_INFO(puts("Filtering disabled."));
    #endif

    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif

    for (size_t nt = 0; nt < Nt - 1; ++nt)
    {
        #ifdef ENABLE_FFT_FILTER
            evo_flags.filter = 1;
        #endif
        if (nt % pars.file.skip == 0)
        {
            evo_flags.compute_pow_spec = 1;
        }

        // k1 & a1
        mk_velocities(t, field, k1);
        #ifdef ENABLE_FFT_FILTER
            evo_flags.filter = 0;
        #endif
        evo_flags.compute_pow_spec = 0;

        // k2 & a2
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k1[i] / 2.0;
        }
        mk_velocities(t + dt / 2.0, tmp_k, k2);

        // k3 & a3
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k2[i] / 2.0;
        }
        mk_velocities(t + dt / 2.0, tmp_k, k3);

        // k4 & a4
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            tmp_k[i] = field[i] + dt * k3[i];
        }
        mk_velocities(t + dt, tmp_k, k4);

        rho = mk_rho(field);

        if (nt % pars.file.skip == 0)
        {
            save();
        }

        // perform one time step for the field and a
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i)
        {
            field[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }

        t += dt;
        pars.t.t = t;
    }

    // make sure to write out last time slice
    rho = mk_rho(field);
    save();

    // info about last timeslice
    if (fabs(pars.t.tf - pars.t.t) > 1e-10)
    {
        RUNTIME_INFO(fputs("The time of the last step does not coincide "
                            "with the specified final time.", stderr));
    }

    RUNTIME_INFO(puts("Finished rk4"));
    #ifdef SHOW_TIMING_INFO
    double end = get_wall_time();
    double secs = end - start;
    RUNTIME_INFO(printf("time: %f seconds\n\n", secs));
    #endif

    fftw_free(k1);
    fftw_free(k2);
    fftw_free(k3);
    fftw_free(k4);
    fftw_free(tmp_k);
}
