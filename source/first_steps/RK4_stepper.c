#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"

void run_RK4_stepper(parameters_t *pars) {
	size_t N  = pars->Nx;
	size_t N2 = 2 * N;
	size_t Nt = pars->Nt;
	double dt = pars->dt;

	double a1, a2, a3, a4;
	double *k1, *k2, *k3, *k4, *tmp;
	k1  = malloc(N2 * sizeof *k1);
	k2  = malloc(N2 * sizeof *k2);
	k3  = malloc(N2 * sizeof *k3);
	k4  = malloc(N2 * sizeof *k4);
	tmp = malloc(N2 * sizeof *tmp);

	if (!(k1 && k2 && k3 && k4 && tmp))
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}
	size_t os, new_os;

	DEBUG(puts("Starting RK4 time evolution with:"));
	DEBUG(printf("initial time: %f\n", pars->ti));
	DEBUG(printf("final time: %f\n", pars->tf));
	DEBUG(printf("time step dt: %f\n", dt));
	DEBUG(printf("number of steps: %zu\n", Nt));
	#ifdef SPECTRAL_OPERATOR_DERIVATIVE
		DEBUG(puts("Using multiplication (lapack) by spectral"
				   "operators for spatial derivatives."));
	#endif
	#ifdef FFT_DERIVATIVE
		DEBUG(puts("Using DFT (fftw3) for spatial derivatives."));
	#endif
	#ifdef ENABLE_FFT_FILTER
		DEBUG(puts("Frequency cutoff filtering enabled."));
	#else
		DEBUG(puts("Filtering disabled."));
	#endif

	clock_t start = clock();

	for (size_t nt = 0; nt < Nt - 1; ++nt)
	{
		os = nt * N2;

		// k1 & a1
		mk_field_velocity(field + os, nt, k1, N);
		a1 = get_a_velocity(nt);

		// k2 & k2
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k1[i] / 2.0;
		}
		mk_field_velocity(tmp, nt, k2, N);
		a2 = get_a_velocity(nt);

		// k3 & a3
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k2[i] / 2.0;
		}
		mk_field_velocity(tmp, nt, k3, N);
		a3 = get_a_velocity(nt);

		// k4 & a4
		for (size_t i = 0; i < N2; ++i)
		{
			tmp[i] = field[os+i] + dt * k3[i];
		}
		mk_field_velocity(tmp, nt, k4, N);
		a4 = get_a_velocity(nt);

		// perform one time step for the field and a
		new_os = os + N2;
		for (size_t i = 0; i < N2; ++i)
		{
			field[new_os+i] = field[os+i] +
						dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}

		frw_a[nt + 1] = frw_a[nt] + dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;

		rho[nt] = mk_rho(field + os, frw_a[nt], N);

		// apply filter
		#ifdef ENABLE_FFT_FILTER
		fft_apply_filter(field + new_os, N);
		#endif

		#ifdef CHECK_FOR_NAN
			for (size_t i = 0; i < N2; ++i)
			{
				if (isnan(field[new_os+i]))
				{
					fprintf(stderr,
						"A nan value was discovered in timestep: %zu \n", nt);
				}
			}
			if (isnan(a[nt]))
			{
				fprintf(stderr,
						"A nan value was discovered in timestep: %zu \n", nt);
			}
		#endif
	}

	// compute the final 00 component of the stress energy
	rho[Nt - 1] = mk_rho(field + 2 * N * (Nt - 1), frw_a[Nt - 1], N);

	clock_t end = clock();

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(tmp);

	double secs = (double)(end - start) / CLOCKS_PER_SEC;
	DEBUG(printf("Finished RK4 time evolution in: %f seconds.\n\n", secs));
}

void mk_field_velocity(double *f, size_t nt, double *result, size_t N) {
	size_t N2 = 2 * N;
	double a = frw_a[nt];
	rho[nt] = mk_rho(f, a, N);

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

	for (size_t i = N; i < N2; ++i)
	{
		result[i] /= (a * a);
		result[i] -= ( 3.0 * get_hubble(nt) * f[i]
						+ potential_prime_term(f[i - N]) );
	}
}

/*
A small selection of potentials one can try, make sure to set the corresponding
potential_prime_term, the derivative is not computed automatically yet
TODO: change that?
*/
inline double potential(double f){
	// return MASS * MASS * f * f / 2.0;
	// return MASS * MASS * f * f / 2.0 + COUPLING * f * f * f * f / 24.0;
	return 0.0;
}

inline double potential_prime_term(double f) {
	// return MASS * MASS * f;
	// return MASS * MASS * f + COUPLING * f * f * f / 6.0;
	// return 20.0 * tanh(pow(f, 50));
	return 0.0;
}

inline double get_a_velocity(size_t nt) {
	return frw_a[nt] * get_hubble(nt);
}

inline double get_hubble(size_t nt) {
	return sqrt(rho[nt] / 3.0);
}

/*
TODO: compute average 00 component of stress energy at current time:
	* what about units
	* is my formula correct (sign)
*/
double mk_rho(double *f, double a, size_t N) {
	double T00 = 0.0;

	double *fD1 = malloc(N * sizeof *fD1);
	fft_D1(f, fD1, N);

	double ft, fxa;
	for (size_t i = 0; i < N; ++i)
	{
		ft = f[N + i];
		fxa = fD1[i] / a;
		T00 += (ft * ft + fxa * fxa) / 2. + potential(f[i]);
	}

	free(fD1);
	return T00 / N;
}