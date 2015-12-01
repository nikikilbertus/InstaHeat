#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "dopri853_stepper.h"
#include "main.h"
#include "RK4_stepper.h"
#include "evolution_toolkit.h"

dopri853_control_t dp;
dopri853_values_t dpv;

void integrate(parameters_t *pars) {
	initialize_dopri853(pars);
	allocate_dopri853_values();
	RUNTIME_INFO(puts("Starting dopri853 integration with:"));
	RUNTIME_INFO(printf("initial time: %f\n", dp.ti));
	RUNTIME_INFO(printf("final time: %f\n", dp.tf));
	RUNTIME_INFO(printf("initial time step dt: %f\n", dp.dt));
	RUNTIME_INFO(printf("minimal time step dt: %f\n", dp.dt_min));
	RUNTIME_INFO(printf("max number of steps: %zu\n", dp.MAX_STEPS));
	RUNTIME_INFO(printf("relative tolerance: %f, "
						"absolute tolerance: %f\n", dp.r_tol, dp.a_tol));
	RUNTIME_INFO(puts("Using DFT (fftw3) for spatial derivatives."));
	#ifdef ENABLE_FFT_FILTER
		RUNTIME_INFO(puts("Frequency cutoff filtering enabled.\n"));
	#else
		RUNTIME_INFO(puts("Filtering disabled.\n"));
	#endif

	df_a = mk_velocities_new(dp.t, field, f_a, dfield, pars);
	if (dp.dense)
	{
		//out.out(-1,x,y,s,h)
	}
	else
	{
		//out.save(x,y)
	}
	for (dp.n_stp = 0; dp.n_stp < dp.MAX_STEPS; ++dp.n_stp)
	{
		if (dp.t + dp.dt * 1.0001 > dp.tf)
		{
			dp.dt = dp.tf - dp.t;
			RUNTIME_INFO(printf("overshoot, new dt = %f\n", dp.dt));
		}
		perform_step(dp.dt, pars);
		if (dp.dt_did == dp.dt)
		{
			++dp.n_ok;
		}
		else
		{
			++dp.n_bad;
		}
		if (dp.dense)
		{
			// out.out(dp.n_stp,x,y,s,s.hdid)
		}
		else
		{
			// out.save(x,y)
		}
		if (dp.t >= dp.tf)
		{
			// for (size_t i = 0; i < dp.Ntot2; ++i)
			// {
			// 	ystart[i] = y[i];
			// }
			// if (...)
			// {
			// 	out.save(x,y)
			// }
			break;
		}
		if (fabs(dp.dt_next) <= dp.dt_min)
		{
			fputs("Step size too small.", stderr);
			exit(EXIT_FAILURE);
		}
		dp.dt = dp.dt_next;
		RUNTIME_INFO(printf("step: %d, next dt: %f\n", dp.n_stp, dp.dt));
	}
	RUNTIME_INFO(puts("finished dopri853 with:"));
	RUNTIME_INFO(printf("steps: %d\n", dp.n_stp + 1));
	RUNTIME_INFO(printf("good steps: %d\n", dp.n_ok));
	RUNTIME_INFO(printf("bad steps: %d\n\n", dp.n_bad));
	destroy_dopri853_values();
}

void initialize_dopri853(parameters_t *pars) {
	dp.t = pars->t.ti;
    dp.t_old = pars->t.ti;
    dp.ti = pars->t.ti;
    dp.tf = pars->t.tf;
    dp.dt = pars->t.dt;
    dp.dt_did = 0.0;
    dp.dt_next = pars->t.dt;
    dp.dt_min = 1.0e-5;
    dp.Ntot2 = 2 * (pars->x.N * pars->y.N * pars->z.N);
    dp.MAX_STEPS = 50000;
    dp.n_stp = 0;
    dp.n_ok = 0;
    dp.n_bad = 0;
    dp.beta  = 0.0;
 	dp.alpha = 1.0/8.0 - dp.beta * 0.2;
 	dp.safe = 0.9;
 	dp.minscale = 0.333;
 	dp.maxscale = 6.0;
	dp.a_tol = 1.0e-5;
	dp.r_tol = dp.a_tol;
	dp.err_old = 1.0e-4;
	dp.reject = 0;
	dp.eps = DBL_EPSILON;
    //dp.dense ;
    RUNTIME_INFO(puts("Initialized dopri853 parameters.\n"));
}

void allocate_dopri853_values() {
	size_t N = dp.Ntot2;

	dpv.k2    = fftw_malloc(N * sizeof *dpv.k2);
	dpv.k3    = fftw_malloc(N * sizeof *dpv.k3);
	dpv.k4    = fftw_malloc(N * sizeof *dpv.k4);
	dpv.k5    = fftw_malloc(N * sizeof *dpv.k5);
	dpv.k6    = fftw_malloc(N * sizeof *dpv.k6);
	dpv.k7    = fftw_malloc(N * sizeof *dpv.k7);
	dpv.k8    = fftw_malloc(N * sizeof *dpv.k8);
	dpv.k9    = fftw_malloc(N * sizeof *dpv.k9);
	dpv.k10   = fftw_malloc(N * sizeof *dpv.k10);
    dpv.k_tmp = fftw_malloc(N * sizeof *dpv.k_tmp);

    dpv.yerr  = fftw_malloc(N * sizeof *dpv.yerr);
    dpv.yerr2 = fftw_malloc(N * sizeof *dpv.yerr2);

    dpv.rcont1 = fftw_malloc(N * sizeof *dpv.rcont1);
    dpv.rcont2 = fftw_malloc(N * sizeof *dpv.rcont2);
    dpv.rcont3 = fftw_malloc(N * sizeof *dpv.rcont3);
    dpv.rcont4 = fftw_malloc(N * sizeof *dpv.rcont4);
    dpv.rcont5 = fftw_malloc(N * sizeof *dpv.rcont5);
    dpv.rcont6 = fftw_malloc(N * sizeof *dpv.rcont6);
    dpv.rcont7 = fftw_malloc(N * sizeof *dpv.rcont7);
    dpv.rcont8 = fftw_malloc(N * sizeof *dpv.rcont8);

    if (!(dpv.k2 && dpv.k3 && dpv.k4 && dpv.k5 && dpv.k6 && dpv.k7 && dpv.k8 &&
    	  dpv.k9 && dpv.k10 && dpv.k_tmp && dpv.yerr && dpv.yerr2 &&
    	  dpv.rcont1 && dpv.rcont2 && dpv.rcont3 && dpv.rcont4 && dpv.rcont5 &&
    	  dpv.rcont6 && dpv.rcont7 && dpv.rcont8))
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }
    RUNTIME_INFO(puts("Allocated memory for dopri853 variables.\n"));
}

void destroy_dopri853_values() {
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

    free(dpv.rcont1);
    free(dpv.rcont2);
    free(dpv.rcont3);
    free(dpv.rcont4);
    free(dpv.rcont5);
    free(dpv.rcont6);
    free(dpv.rcont7);
    free(dpv.rcont8);
    RUNTIME_INFO(puts("Freed memory of dopri853 variables.\n"));
}

void perform_step(const double dt_try, parameters_t *pars) {
	double dt = dt_try;
	for ( ; ; )
	{
		try_step(dt, pars);
		double err = error(dt);
		if (success(err, &dt))
		{
			break;
		}
		if (fabs(dt) <= fabs(dp.t) * dp.eps)
		{
			fputs("Stepsize underflow", stderr);
			exit(EXIT_FAILURE);
		}
	}
	#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 1;
	#endif
	df_a_new = mk_velocities_new(dp.t + dt, field_new, f_a_new, dfield_new, pars);
	#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 0;
	#endif
	if (dp.dense)
	{
		// prepare_dense_output(dt);
	}
	#pragma omp parallel for
	for (size_t i = 0; i < dp.Ntot2; ++i)
	{
		dfield[i] = dfield_new[i];
		field[i] = field_new[i];
	}
	df_a = df_a_new;
	f_a = f_a_new;
	dp.t_old = dp.t;
	dp.t += (dp.dt_did = dt);
	//dp.dt_next = con.dt_next;
}

double error(const double dt) {
	size_t N = dp.Ntot2;

	double err = 0.0, err2 = 0.0, sk, deno;
	// TODO[performance] parallelize with reduction to err2, err
	for (size_t i = 0; i < N; ++i)
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
	return dt * err * sqrt(1.0 / (N * deno));
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

void try_step(const double dt, parameters_t *pars) {
	size_t i;
	double t = dp.t;
	size_t Ntot2 = dp.Ntot2;
	// ------------ 1 ------------
	// dpv.a1 = mk_velocities_new(t, field, f_a, dpv.k1, pars);

	// ------------ 2 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt * dpc.a21 * dfield[i];
	}
	dpv.a_tmp = f_a + dt * dpc.a21 * df_a;
	dpv.a2 = mk_velocities_new(t + dpc.c2 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k2, pars);

	// ------------ 3 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a31 * dfield[i] + dpc.a32 * dpv.k2[i]);
	}
	dpv.a_tmp = f_a + dt * (dpc.a31 * df_a + dpc.a32 * dpv.a2);
	dpv.a3 = mk_velocities_new(t + dpc.c3 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k3, pars);

	// ------------ 4 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a41 * dfield[i] + dpc.a43 * dpv.k3[i]);
	}
	dpv.a_tmp = f_a + dt * (dpc.a41 * df_a + dpc.a43 * dpv.a3);
	dpv.a4 = mk_velocities_new(t + dpc.c4 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k4, pars);

	// ------------ 5 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a51 * dfield[i] + dpc.a53 * dpv.k3[i] + dpc.a54 * dpv.k4[i]);
	}
	dpv.a_tmp = f_a + dt * (dpc.a51 * df_a + dpc.a53 * dpv.a3 + dpc.a54 * dpv.a4);
	dpv.a5 = mk_velocities_new(t + dpc.c5 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k5, pars);

	// ------------ 6 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a61 * dfield[i] + dpc.a64 * dpv.k4[i] + dpc.a65 * dpv.k5[i]);
	}
	dpv.a_tmp = f_a + dt * (dpc.a61 * df_a + dpc.a64 * dpv.a4 + dpc.a65 * dpv.a5);
	dpv.a6 = mk_velocities_new(t + dpc.c6 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k6, pars);

	// ------------ 7 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a71 * dfield[i] + dpc.a74 * dpv.k4[i] + dpc.a75 * dpv.k5[i] +
			 dpc.a76 * dpv.k6[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a71 * df_a + dpc.a74 * dpv.a4 + dpc.a75 * dpv.a5 +
			 dpc.a76 * dpv.a6);
	dpv.a7 = mk_velocities_new(t + dpc.c7 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k7, pars);

	// ------------ 8 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a81 * dfield[i] + dpc.a84 * dpv.k4[i] + dpc.a85 * dpv.k5[i] +
			 dpc.a86 * dpv.k6[i] + dpc.a87 * dpv.k7[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a81 * df_a + dpc.a84 * dpv.a4 + dpc.a85 * dpv.a5 +
			 dpc.a86 * dpv.a6 + dpc.a87 * dpv.a7);
	dpv.a8 = mk_velocities_new(t + dpc.c8 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k8, pars);

	// ------------ 9 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a91 * dfield[i] + dpc.a94 * dpv.k4[i] + dpc.a95 * dpv.k5[i] +
			 dpc.a96 * dpv.k6[i] + dpc.a97 * dpv.k7[i] + dpc.a98 * dpv.k8[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a91 * df_a + dpc.a94 * dpv.a4 + dpc.a95 * dpv.a5 +
			 dpc.a96 * dpv.a6 + dpc.a97 * dpv.a7 + dpc.a98 * dpv.a8);
	dpv.a9 = mk_velocities_new(t + dpc.c9 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k9, pars);

	// ------------ 10 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a101 * dfield[i] + dpc.a104 * dpv.k4[i] + dpc.a105 * dpv.k5[i] +
			 dpc.a106 * dpv.k6[i] + dpc.a107 * dpv.k7[i] + dpc.a108 * dpv.k8[i] +
			 dpc.a109 * dpv.k9[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a101 * df_a + dpc.a104 * dpv.a4 + dpc.a105 * dpv.a5 +
			 dpc.a106 * dpv.a6 + dpc.a107 * dpv.a7 + dpc.a108 * dpv.a8 +
			 dpc.a109 * dpv.a9);
	dpv.a10 = mk_velocities_new(t + dpc.c10 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k10, pars);

	// ------------ 11 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a111 * dfield[i] + dpc.a114 * dpv.k4[i] + dpc.a115 * dpv.k5[i] +
			 dpc.a116 * dpv.k6[i] + dpc.a117 * dpv.k7[i] + dpc.a118 * dpv.k8[i] +
			 dpc.a119 * dpv.k9[i] + dpc.a1110 * dpv.k10[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a111 * df_a + dpc.a114 * dpv.a4 + dpc.a115 * dpv.a5 +
			 dpc.a116 * dpv.a6 + dpc.a117 * dpv.a7 + dpc.a118 * dpv.a8 +
			 dpc.a119 * dpv.a9 + dpc.a1110 * dpv.a10);
	dpv.a2 = mk_velocities_new(t + dpc.c11 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k2, pars);

	// ------------ new dt ------------
	double tpdt = t + dt;

	// ------------ 12 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a121 * dfield[i] + dpc.a124 * dpv.k4[i] + dpc.a125 * dpv.k5[i] +
			 dpc.a126 * dpv.k6[i] + dpc.a127 * dpv.k7[i] + dpc.a128 * dpv.k8[i] +
			 dpc.a129 * dpv.k9[i] + dpc.a1210 * dpv.k10[i] + dpc.a1211 * dpv.k2[i]);
	}
	dpv.a_tmp = f_a + dt *
			(dpc.a121 * df_a + dpc.a124 * dpv.a4 + dpc.a125 * dpv.a5 +
			 dpc.a126 * dpv.a6 + dpc.a127 * dpv.a7 + dpc.a128 * dpv.a8 +
			 dpc.a129 * dpv.a9 + dpc.a1210 * dpv.a10 + dpc.a1211 * dpv.a2);
	dpv.a3 = mk_velocities_new(tpdt, dpv.k_tmp, dpv.a_tmp, dpv.k3, pars);

	// ------------ step ahead ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k4[i] = dpc.b1 * dfield[i] + dpc.b6 * dpv.k6[i] + dpc.b7 * dpv.k7[i] +
					dpc.b8 * dpv.k8[i] + dpc.b9 * dpv.k9[i] + dpc.b10 * dpv.k10[i] +
					dpc.b11 * dpv.k2[i] + dpc.b12 * dpv.k3[i];

		field_new[i] = field[i] + dt * dpv.k4[i];
	}
	dpv.a4 = dpc.b1 * df_a + dpc.b6 * dpv.a6 + dpc.b7 * dpv.a7 +
			 dpc.b8 * dpv.a8 + dpc.b9 * dpv.a9 + dpc.b10 * dpv.a10 +
			 dpc.b11 * dpv.a2 + dpc.b12 * dpv.a3;
	f_a_new = f_a + dt * dpv.a4;

	// ------------ error estimates ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.yerr[i]  = dpv.k4[i] - dpc.bhh1 * dfield[i] - dpc.bhh2 * dpv.k9[i] -
					   dpc.bhh3 * dpv.k3[i];
		dpv.yerr2[i] = dpc.er1 * dfield[i] + dpc.er6 * dpv.k6[i] +
					   dpc.er7 * dpv.k7[i] + dpc.er8 * dpv.k8[i] +
					   dpc.er9 * dpv.k9[i] + dpc.er10 * dpv.k10[i] +
					   dpc.er11 * dpv.k2[i] + dpc.er12 * dpv.k3[i];
	}
}

/*
compute the right hand side of the pde, ie the first order temporal derivatives
*/
double mk_velocities_new(double t, double *f, double a, double *result, parameters_t *pars) {
	size_t Nx  = pars->x.N;
	size_t Ny  = pars->y.N;
	size_t Nz  = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t Ntot2 = 2 * Ntot;

	double current_rho = mk_rho(f, a, pars);
	double hubble = sqrt(current_rho / 3.0);

	#pragma omp parallel for
	for (size_t i = 0; i < Ntot; ++i)
	{
		result[i] = f[Ntot + i];
	}

	#pragma omp parallel for
	for (size_t i = Ntot; i < Ntot2; ++i)
	{
		result[i] = dtmp_lap[i - Ntot] / (a * a);
		result[i] -= ( 3.0 * hubble * f[i]
						+ potential_prime(f[i - Ntot]) );
	}
	return a * hubble;
}