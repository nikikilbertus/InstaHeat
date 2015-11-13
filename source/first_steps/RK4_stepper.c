#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"

void run_RK4_stepper(parameters_t *pars) {
	size_t Nx  = pars->x.N;
	size_t Ny  = pars->y.N;
	size_t Nz  = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t Ntot2 = 2 * Ntot;
	size_t Nt = pars->Nt;
	double dt = pars->dt;

	double a1, a2, a3, a4, tmp_a;
	double *k1, *k2, *k3, *k4, *tmp_k;

	k1    = malloc(N2 * sizeof *k1);
	k2    = malloc(N2 * sizeof *k2);
	k3    = malloc(N2 * sizeof *k3);
	k4    = malloc(N2 * sizeof *k4);
	tmp_k = malloc(N2 * sizeof *tmp_k);

	if (!(k1 && k2 && k3 && k4 && tmp_k))
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}
	size_t os, new_os;

	RUNTIME_INFO(puts("Starting RK4 time evolution with:"));
	RUNTIME_INFO(printf("initial time: %f\n", pars->ti));
	RUNTIME_INFO(printf("final time: %f\n", pars->tf));
	RUNTIME_INFO(printf("time step dt: %f\n", dt));
	RUNTIME_INFO(printf("number of steps: %zu\n", Nt));
	RUNTIME_INFO(puts("Using DFT (fftw3) for spatial derivatives."));

#ifdef ENABLE_FFT_FILTER
		RUNTIME_INFO(puts("Frequency cutoff filtering enabled."));
#else
		RUNTIME_INFO(puts("Filtering disabled."));
#endif

	clock_t start = clock();

	for (size_t nt = 0; nt < Nt - 1; ++nt)
	{
		os = nt * Ntot2;

		// k1 & a1
		a1 = mk_velocities(field + os, frw_a[nt], k1, pars);

		// k2 & k2
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[os+i] + dt * k1[i] / 2.0;
		}
		tmp_a = frw_a[nt] + dt * a1 / 2.0;
		a2 = mk_velocities(tmp_k, tmp_a, k2, pars);

		// k3 & a3
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[os+i] + dt * k2[i] / 2.0;
		}
		tmp_a = frw_a[nt] + dt * a2 / 2.0;
		a3 = mk_velocities(tmp_k, tmp_a, k3, pars);

		// k4 & a4
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[os+i] + dt * k3[i];
		}
		tmp_a = frw_a[nt] + dt * a3;
		a4 = mk_velocities(tmp_k, tmp_a, k4, pars);

		// perform one time step for the field and a
		new_os = os + Ntot2;
		for (size_t i = 0; i < Ntot2; ++i)
		{
			field[new_os+i] = field[os+i] +
						dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}

		frw_a[nt + 1] = frw_a[nt] + dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;

		rho[nt] = mk_rho(field + os, frw_a[nt], pars);

		// apply filter
#ifdef ENABLE_FFT_FILTER
		fft_apply_filter(field + new_os, pars);
#endif

#ifdef CHECK_FOR_NAN
			for (size_t i = 0; i < Ntot2; ++i)
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
	rho[Nt - 1] = mk_rho(field + Ntot2 * (Nt - 1), frw_a[Nt - 1], pars);

	clock_t end = clock();

	fftw_free(k1);
	fftw_free(k2);
	fftw_free(k3);
	fftw_free(k4);
	fftw_free(tmp_k);

	double secs = (double)(end - start) / CLOCKS_PER_SEC;
	RUNTIME_INFO(printf("Finished time evolution in: %f seconds.\n\n", secs));
}

/*
compute the right hand side of the pde (first order in time)
*/
double mk_velocities(double *f, double a, double *result, parameters_t *pars) {
	size_t Nx  = pars->x.N;
	size_t Ny  = pars->y.N;
	size_t Nz  = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t Ntot2 = 2 * Ntot;

	double current_rho = mk_rho(f, a, pars);
	double hubble = sqrt(current_rho / 3.0);

	for (size_t i = 0; i < Ntot; ++i)
	{
		result[i] = f[Ntot + i];
	}

	fft_D2(f, result + Ntot, pars);

	for (size_t i = Ntot; i < Ntot2; ++i)
	{
		result[i] /= (a * a);
		result[i] -= ( 3.0 * hubble * f[i]
						+ potential_prime_term(f[i - Ntot]) );
	}
	return a * hubble;
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

/*
compute average 00 component of stress energy
*/
double mk_rho(double *f, double a, parameters_t *pars) {
	size_t Nx  = pars->x.N;
	size_t Ny  = pars->y.N;
	size_t Nz  = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;

	double T00 = 0.0;

	double *f_grad2 = malloc(Ntot * sizeof *f_grad2);
	if (!f_grad2)
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}

	mk_gradient_squared(f, f_grad2, pars);

	double ft, f_grad_a;
	for (size_t i = 0; i < Ntot; ++i)
	{
		ft = f[Ntot + i];
		f_grad_a = f_grad2[i] / a;
		T00 += (ft * ft + f_grad_a * f_grad_a) / 2. + potential(f[i]);
	}

	free(f_grad2);
	return T00 / Ntot;
}