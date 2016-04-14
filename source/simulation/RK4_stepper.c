#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

/**
 * @file RK4_stepper.c
 * @brief A simple standard 4th order explicit Runge Kutta integrator with fixed
 * stepsize.
 */

/**
 * @brief Standard 4th order explicit Runge Kutta integrator with fixed
 * stepsize.
 *
 * @note The stepsize might be changed for the very last time step in order to
 * exactly end at the given final time.
 */
void run_rk4()
{
    const size_t Ntot = pars.Ntot;
    const size_t Nall = pars.Nall;
    const size_t Nt = pars.t.Nt;
    double dt = pars.t.dt;
    double t = pars.t.t;

    double *k1, *k2, *k3, *k4, *tmp_k;
    k1    = fftw_malloc(Ntot * sizeof *k1);
    k2    = fftw_malloc(Ntot * sizeof *k2);
    k3    = fftw_malloc(Ntot * sizeof *k3);
    k4    = fftw_malloc(Ntot * sizeof *k4);
    tmp_k = fftw_malloc(Nall * sizeof *tmp_k);

    if (!(k1 && k2 && k3 && k4 && tmp_k)) {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }

    INFO(puts("Starting RK4 time evolution with:"));
    INFO(printf("initial time: %f\n", pars.t.ti));
    INFO(printf("final time: %f\n", pars.t.tf));
    INFO(printf("time step dt: %f\n", dt));
    INFO(printf("number of steps: %zu\n", Nt));

    #ifdef SHOW_TIMING_INFO
    double secs = -get_wall_time();
    #endif

    for (size_t nt = 0; t < pars.t.tf; ++nt) {
        #ifdef ENABLE_FFT_FILTER
        apply_filter_real(field);
        #endif

        // to precisely reach final time in the last step, change dt
        if (t + dt * 1.0001 > pars.t.tf) {
            dt = pars.t.tf - t;
            pars.t.dt = dt;
            INFO(printf("overshoot, new dt = %f\n", dt));
        }

        // step 1 (and write out data if required)
        if (nt % pars.file.skip == 0) {
            #ifdef OUTPUT_PHI_PS
            evo_flags.compute_pow_spec = 1;
            #endif
            #ifdef OUTPUT_CONSTRAINTS
            evo_flags.compute_cstr = 1;
            #endif
            mk_rhs(t, field, k1);
            evo_flags.compute_pow_spec = 0;
            evo_flags.compute_cstr = 0;
            mk_summary();
            save();
        } else {
            mk_rhs(t, field, k1);
        }

        // step 2
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i) {
            tmp_k[i] = field[i] + dt * k1[i] / 2.0;
        }
        mk_rhs(t + dt / 2.0, tmp_k, k2);

        // step 3
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i) {
            tmp_k[i] = field[i] + dt * k2[i] / 2.0;
        }
        mk_rhs(t + dt / 2.0, tmp_k, k3);

        // step 4
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i) {
            tmp_k[i] = field[i] + dt * k3[i];
        }
        mk_rhs(t + dt, tmp_k, k4);

        // perform time step
        #pragma omp parallel for
        for (size_t i = 0; i < Ntot; ++i) {
            field[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }

        t += dt;
        pars.t.t = t;
        t_out.tmp[0] = t;
    }

    // make sure to write out last time slice
    prepare_and_save_timeslice();

    // info about last timeslice
    if (fabs(pars.t.tf - pars.t.t) > 1e-14) {
        INFO(puts("The time of the last step does not coincide "
                          "with the specified final time."));
    }

    INFO(puts("Finished rk4"));
    #ifdef SHOW_TIMING_INFO
    secs += get_wall_time();
    INFO(printf("time: %f seconds\n\n", secs));
    INFO(puts("Writing simulation meta data to disk\n"));
    double val[1] = {secs};
    h5_write_parameter(H5_RUNTIME_STEPPER_NAME, val, 1);
    #endif

    fftw_free(k1);
    fftw_free(k2);
    fftw_free(k3);
    fftw_free(k4);
    fftw_free(tmp_k);
}
