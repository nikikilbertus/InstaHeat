#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "dopri853_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

dopri853_control_t dp;
dopri853_values_t dpv;

void initialize_dopri853() {
    dp.t = pars.t.ti;
    dp.t_old = pars.t.ti;
    dp.ti = pars.t.ti;
    dp.tf = pars.t.tf;
    dp.dt = pars.t.dt;
    dp.dt_did = 0.0;
    dp.dt_next = pars.t.dt;
    dp.dt_min = MINIMAL_DELTA_T;
    dp.max_steps = MAX_STEPS;
    dp.n_stp = 0;
    dp.n_ok = 0;
    dp.n_bad = 0;
    dp.beta  = BETA;
    dp.alpha = 1.0/8.0 - dp.beta * 0.2;
    dp.safe = SAFE;
    dp.minscale = SMALLEST_SCALING;
    dp.maxscale = LARGEST_SCALING;
    dp.a_tol = ABSOLUTE_TOLERANCE;
    dp.r_tol = RELATIVE_TOLERANCE;
    dp.err_old = 1.0e-4;
    dp.reject = 0;
    dp.eps = DBL_EPSILON;
    RUNTIME_INFO(puts("Initialized dopri853 parameters.\n"));
}

void run_dopri853() {
    initialize_dopri853();
    allocate_dopri853_values();
    RUNTIME_INFO(puts("Starting dopri853 integration with:"));
    RUNTIME_INFO(printf("initial time: %f\n", dp.ti));
    RUNTIME_INFO(printf("final time: %f\n", dp.tf));
    RUNTIME_INFO(printf("initial time step dt: %f\n", dp.dt));
    RUNTIME_INFO(printf("minimal time step dt: %f\n", dp.dt_min));
    RUNTIME_INFO(printf("max number of steps: %zu\n", dp.max_steps));
    RUNTIME_INFO(printf("relative tolerance: %.15f\n", dp.r_tol));
    RUNTIME_INFO(printf("absolute tolerance: %.15f\n", dp.a_tol));
    RUNTIME_INFO(puts("Using DFT (fftw3) for spatial derivatives."));

    #ifdef ENABLE_FFT_FILTER
    RUNTIME_INFO(puts("Frequency cutoff filtering enabled.\n"));
    #else
    RUNTIME_INFO(puts("Filtering disabled.\n"));
    #endif

    evo_flags.compute_pow_spec = 1;
    mk_rhs(dp.t, field, dfield);
    evo_flags.compute_pow_spec = 0;
    mk_means_and_variances();
    save();

    #ifdef SHOW_TIMING_INFO
    double secs = -get_wall_time();
    #endif

    for (dp.n_stp = 0; dp.n_stp < dp.max_steps; ++dp.n_stp)
    {
        if (dp.t + dp.dt * 1.0001 > dp.tf)
        {
            dp.dt = dp.tf - dp.t;
            RUNTIME_INFO(printf("overshoot, new dt = %f\n", dp.dt));
        }
        if (perform_step(dp.dt))
        {
            break;
        }
        if (dp.dt_did == dp.dt)
        {
            ++dp.n_ok;
        }
        else
        {
            ++dp.n_bad;
        }
        #ifdef DEBUG
        RUNTIME_INFO(printf("did step: %d with dt: %f\n", dp.n_stp, dp.dt_did));
        #endif

        if ((dp.n_stp + 1) % pars.file.skip == 0)
        {
            mk_means_and_variances();
            save();
        }
        if (dp.t >= dp.tf)
        {
            break;
        }
        if (fabs(dp.dt_next) <= dp.dt_min)
        {
            fputs("Step size too small. Still exiting gracefully.\n", stderr);
            break;
        }
        dp.dt = dp.dt_next;
    }

    size_t index = pars.file.index;
    if (index != 0 && time_buf[index - 1] < dp.t)
    {
        prepare_and_save_timeslice();
    }

    #ifdef SHOW_TIMING_INFO
    secs += get_wall_time();
    #endif

    free_dopri853_values();

    RUNTIME_INFO(puts("Finished dopri853"));
    #ifdef SHOW_TIMING_INFO
    RUNTIME_INFO(printf("time: %f seconds\n", secs));
    #endif
    RUNTIME_INFO(printf("steps: %d\n", dp.n_stp + 1));
    RUNTIME_INFO(printf("good steps: %d\n", dp.n_ok));
    RUNTIME_INFO(printf("bad steps: %d\n\n", dp.n_bad));
}

int perform_step(const double dt_try) {
    size_t Ntot = 3 * pars.N + 1;
    double dt = dt_try;
    for ( ; ; )
    {
        try_step(dt);
        double err = error(dt);
        if (success(err, &dt))
        {
            break;
        }
        if (fabs(dt) <= fabs(dp.t) * dp.eps)
        {
            fputs("Stepsize underflow. Still exiting gracefully.\n", stderr);
            return 1;
        }
    }
    #ifdef ENABLE_FFT_FILTER
    apply_filter_real(field_new);
    #endif
    if ((dp.n_stp + 1) % pars.file.skip == 0)
    {
        evo_flags.compute_pow_spec = 1;
    }
    mk_rhs(dp.t + dt, field_new, dfield_new);
    evo_flags.compute_pow_spec = 0;

    #pragma omp parallel for
    for (size_t i = 0; i < Ntot; ++i)
    {
        field[i] = field_new[i];
        dfield[i] = dfield_new[i];
    }
    dp.t_old = dp.t;
    dp.t += (dp.dt_did = dt);
    pars.t.t = dp.t;
    return 0;
}

void try_step(const double dt) {
    size_t i;
    double t = dp.t;
    size_t Ntot = 3 * pars.N + 1;
    // ------------ 1 ------------
    // is already done in perform_step

    // ------------ 2 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt * dpc.a21 * dfield[i];
    }
    mk_rhs(t + dpc.c2 * dt, dpv.k_tmp, dpv.k2);

    // ------------ 3 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a31 * dfield[i] + dpc.a32 * dpv.k2[i]);
    }
    mk_rhs(t + dpc.c3 * dt, dpv.k_tmp, dpv.k3);

    // ------------ 4 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a41 * dfield[i] + dpc.a43 * dpv.k3[i]);
    }
    mk_rhs(t + dpc.c4 * dt, dpv.k_tmp, dpv.k4);

    // ------------ 5 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a51 * dfield[i] + dpc.a53 * dpv.k3[i] + dpc.a54 * dpv.k4[i]);
    }
    mk_rhs(t + dpc.c5 * dt, dpv.k_tmp, dpv.k5);

    // ------------ 6 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a61 * dfield[i] + dpc.a64 * dpv.k4[i] + dpc.a65 * dpv.k5[i]);
    }
    mk_rhs(t + dpc.c6 * dt, dpv.k_tmp, dpv.k6);

    // ------------ 7 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a71 * dfield[i] + dpc.a74 * dpv.k4[i] + dpc.a75 * dpv.k5[i] +
             dpc.a76 * dpv.k6[i]);
    }
    mk_rhs(t + dpc.c7 * dt, dpv.k_tmp, dpv.k7);

    // ------------ 8 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a81 * dfield[i] + dpc.a84 * dpv.k4[i] + dpc.a85 * dpv.k5[i] +
             dpc.a86 * dpv.k6[i] + dpc.a87 * dpv.k7[i]);
    }
    mk_rhs(t + dpc.c8 * dt, dpv.k_tmp, dpv.k8);

    // ------------ 9 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a91 * dfield[i] + dpc.a94 * dpv.k4[i] + dpc.a95 * dpv.k5[i] +
             dpc.a96 * dpv.k6[i] + dpc.a97 * dpv.k7[i] + dpc.a98 * dpv.k8[i]);
    }
    mk_rhs(t + dpc.c9 * dt, dpv.k_tmp, dpv.k9);

    // ------------ 10 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a101 * dfield[i] + dpc.a104 * dpv.k4[i] + dpc.a105 * dpv.k5[i] +
             dpc.a106 * dpv.k6[i] + dpc.a107 * dpv.k7[i] + dpc.a108 * dpv.k8[i] +
             dpc.a109 * dpv.k9[i]);
    }
    mk_rhs(t + dpc.c10 * dt, dpv.k_tmp, dpv.k10);

    // ------------ 11 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a111 * dfield[i] + dpc.a114 * dpv.k4[i] + dpc.a115 * dpv.k5[i] +
             dpc.a116 * dpv.k6[i] + dpc.a117 * dpv.k7[i] + dpc.a118 * dpv.k8[i] +
             dpc.a119 * dpv.k9[i] + dpc.a1110 * dpv.k10[i]);
    }
    mk_rhs(t + dpc.c11 * dt, dpv.k_tmp, dpv.k2);

    // ------------ new dt ------------
    double tpdt = t + dt;

    // ------------ 12 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a121 * dfield[i] + dpc.a124 * dpv.k4[i] + dpc.a125 * dpv.k5[i] +
             dpc.a126 * dpv.k6[i] + dpc.a127 * dpv.k7[i] + dpc.a128 * dpv.k8[i] +
             dpc.a129 * dpv.k9[i] + dpc.a1210 * dpv.k10[i] + dpc.a1211 * dpv.k2[i]);
    }
    mk_rhs(tpdt, dpv.k_tmp, dpv.k3);

    // ------------ step ahead ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.k4[i] = dpc.b1 * dfield[i] + dpc.b6 * dpv.k6[i] + dpc.b7 * dpv.k7[i] +
            dpc.b8 * dpv.k8[i] + dpc.b9 * dpv.k9[i] + dpc.b10 * dpv.k10[i] +
            dpc.b11 * dpv.k2[i] + dpc.b12 * dpv.k3[i];

        field_new[i] = field[i] + dt * dpv.k4[i];
    }

    // ------------ error estimates ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i)
    {
        dpv.yerr[i]  = dpv.k4[i] - dpc.bhh1 * dfield[i] - dpc.bhh2 * dpv.k9[i] -
                        dpc.bhh3 * dpv.k3[i];
        dpv.yerr2[i] = dpc.er1 * dfield[i] + dpc.er6 * dpv.k6[i] +
                       dpc.er7 * dpv.k7[i] + dpc.er8 * dpv.k8[i] +
                       dpc.er9 * dpv.k9[i] + dpc.er10 * dpv.k10[i] +
                       dpc.er11 * dpv.k2[i] + dpc.er12 * dpv.k3[i];
    }
}

double error(const double dt) {
    size_t Ntot = 3 * pars.N + 1;
    double err = 0.0, err2 = 0.0, sk, deno;

    #pragma omp parallel for private(sk) reduction(+:err,err2)
    for (size_t i = 0; i < Ntot; ++i)
    {
        sk = dp.a_tol + dp.r_tol * MAX(fabs(field[i]), fabs(field_new[i]));
        err2 += (dpv.yerr[i] / sk) * (dpv.yerr[i] / sk);
        err  += (dpv.yerr2[i] / sk) * (dpv.yerr2[i] / sk);
    }
    deno = err + 0.01 * err2;
    if (deno <= 0.0)
    {
        deno = 1.0;
    }
    return dt * err * sqrt(1.0 / (Ntot * deno));
}

int success(const double err, double *dt) {
    double beta  = dp.beta;
    double alpha = dp.alpha;
    double safe  = dp.safe;
    double minscale = dp.minscale;
    double maxscale = dp.maxscale;
    double scale;

    if (err <= 1.0)
    {
        if (err == 0.0)
        {
            scale = maxscale;
        }
        else
        {
            scale = safe * pow(err, -alpha) * pow(dp.err_old, beta);
            if (scale < minscale)
            {
                scale = minscale;
            }
            if (scale > maxscale)
            {
                scale = maxscale;
            }
        }
        if (dp.reject)
        {
            dp.dt_next = (*dt) * MIN(scale, 1.0);
        }
        else
        {
            dp.dt_next = (*dt) * scale;
        }
        dp.err_old = MAX(err, 1.0e-4);
        dp.reject = 0;
        return 1;
    }
    else
    {
        scale = MAX(safe * pow(err, -alpha), minscale);
        (*dt) *= scale;
        dp.reject = 1;
        return 0;
    }
}

void allocate_dopri853_values() {
    size_t Ntot = 3 * pars.N + 1;

    dpv.k2    = fftw_malloc(Ntot * sizeof *dpv.k2);
    dpv.k3    = fftw_malloc(Ntot * sizeof *dpv.k3);
    dpv.k4    = fftw_malloc(Ntot * sizeof *dpv.k4);
    dpv.k5    = fftw_malloc(Ntot * sizeof *dpv.k5);
    dpv.k6    = fftw_malloc(Ntot * sizeof *dpv.k6);
    dpv.k7    = fftw_malloc(Ntot * sizeof *dpv.k7);
    dpv.k8    = fftw_malloc(Ntot * sizeof *dpv.k8);
    dpv.k9    = fftw_malloc(Ntot * sizeof *dpv.k9);
    dpv.k10   = fftw_malloc(Ntot * sizeof *dpv.k10);
    dpv.k_tmp = fftw_malloc(Ntot * sizeof *dpv.k_tmp);

    dpv.yerr  = fftw_malloc(Ntot * sizeof *dpv.yerr);
    dpv.yerr2 = fftw_malloc(Ntot * sizeof *dpv.yerr2);

    if (!(dpv.k2 && dpv.k3 && dpv.k4 && dpv.k5 && dpv.k6 && dpv.k7 && dpv.k8 &&
          dpv.k9 && dpv.k10 && dpv.k_tmp && dpv.yerr && dpv.yerr2))
    {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    RUNTIME_INFO(puts("Allocated memory for dopri853 variables.\n"));
}

void free_dopri853_values() {
    free(dpv.k2);
    free(dpv.k3);
    free(dpv.k4);
    free(dpv.k5);
    free(dpv.k6);
    free(dpv.k7);
    free(dpv.k8);
    free(dpv.k9);
    free(dpv.k10);
    free(dpv.k_tmp);

    free(dpv.yerr);
    free(dpv.yerr2);
    RUNTIME_INFO(puts("Freed memory of dopri853 variables.\n"));
}
