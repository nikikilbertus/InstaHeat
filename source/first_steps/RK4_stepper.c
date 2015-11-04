#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <Accelerate/Accelerate.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "setup.h"

void run_RK4_stepper(parameters_t *pars) {
	size_t N  = pars->Nx;
	size_t N2 = 2 * N;
	size_t Nt = pars->Nt;
	double dt = pars->dt;

	double *k1, *k2, *k3, *k4, *tmp;
	k1  = malloc(N2 * sizeof *k1);
	k2  = malloc(N2 * sizeof *k2);
	k3  = malloc(N2 * sizeof *k3);
	k4  = malloc(N2 * sizeof *k4);
	tmp = malloc(N2 * sizeof *tmp);

	//double a1, a2, a3, a4;

	if (!(k1 && k2 && k3 && k4 && tmp))
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}
	size_t os, new_os;

	DEBUG(puts("Starting RK4 time evolution with:\n"));
	DEBUG(printf("initial time: %f\n", pars->ti));
	DEBUG(printf("final time: %f\n", pars->tf));
	DEBUG(printf("time step dt: %f\n", dt));
	DEBUG(printf("number of steps: %zu\n", Nt));
	#ifdef SPECTRAL_OPERATOR_DERIVATIVE
		DEBUG(puts("Using multiplication (lapack) by spectral"
				   "operators for spatial derivatives.\n"));
	#endif
	#ifdef FFT_DERIVATIVE
		DEBUG(puts("Using DFT (fftw3) for spatial derivatives.\n"));
	#endif

	clock_t start = clock();

	for (size_t nt = 0; nt < Nt; ++nt)
	{
		os = nt * N2;

		// k1 & a1
		get_field_velocity(field + os, k1, N);
		// a1 = get_a_velocity(field + os, a[nt], N);

		// k2 & k2
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k1[i] / 2.0;
		}
		get_field_velocity(tmp, k2, N);
		// a2 = get_a_velocity(tmp, a[nt], N);

		// k3 & a3
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k2[i] / 2.0;
		}
		get_field_velocity(tmp, k3, N);
		// a3 = get_a_velocity(tmp, a[nt], N);

		// k4 & a4
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k3[i];
		}
		get_field_velocity(tmp, k4, N);
		// a4 = get_a_velocity(tmp, a[nt], N);

		// perform one time step for the field and a
		new_os = os + N2;
		for (size_t i = 0; i < N2; ++i)
		{
			field[new_os+i] = field[os+i] +
						dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}
		// a[nt+1] = a[nt] + dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;

		#ifdef CHECK_FOR_NAN
			for (size_t i = 0; i < N2; ++i)
			{
				if (isnan(field[new_os+i]))
				{
					fprintf(stderr,
						"A nan value was discovered in timestep: %zu \n", nt);
				}
			}
		#endif
	}

	clock_t end = clock();

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(tmp);

	double secs = (double)(end - start) / CLOCKS_PER_SEC;
	DEBUG(printf("Finished RK4 time evolution in: %f seconds.\n\n", secs));
}

void get_field_velocity(double *f, double *result, size_t N) {
	for (size_t i = 0; i < N; ++i)
	{
		result[i] = f[N + i];
	}
	#ifdef SPECTRAL_OPERATOR_DERIVATIVE
		spectral_op_D2(f, result + N, N);
	#endif
	#ifdef FFT_DERIVATIVE
		fft_D2(f, result + N, N);
	#endif

	for (size_t i = 0; i < N; ++i)
	{
		result[i + N] += potential_term(f[i]);
	}
}

void spectral_op_D2(double *f, double *result, size_t N) {
	matrix_vector(D2, f, result, N);
}

void fft_D2(double *f, double *result, size_t N) {
	size_t nc = N / 2 + 1;
	fftw_complex *tmp;
	fftw_plan p_fw;
	fftw_plan p_bw;
	tmp = fftw_malloc(nc * sizeof *tmp);
	if (!tmp)
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}
	p_fw = fftw_plan_dft_r2c_1d(N, f, tmp, FFTW_ESTIMATE);
	fftw_execute(p_fw);
	double L = pars.b - pars.a;
	double factor = - 4 * PI * PI / (L*L);
	for (size_t i = 0; i < nc; ++i)
	{
		tmp[i] *= (factor * i * i) / N;
	}
	p_bw = fftw_plan_dft_c2r_1d(N, tmp, result, FFTW_ESTIMATE);
	fftw_execute(p_bw);
	fftw_destroy_plan(p_fw);
	fftw_destroy_plan(p_bw);
	fftw_free(tmp);
}

inline double potential_term(double f) {
	// return -MASS*MASS*f - COUPLING * f*f*f / 6.0;
	// return - 20.0 * tanh(pow(f, 50));
	return 0.0;
}

double get_a_velocity(double *f, double a, size_t N) {
	return a * get_hubble(f, N);
}

double get_hubble(double *f, size_t N) {
	/*
	TODO: compute average energy at current time:
		* what exactly do I need
		* how do I average
		* ^0.5 in the end?
		* what about units
	*/
	return 0.0;
}