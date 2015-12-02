#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

void run_RK4_stepper(parameters_t *pars) {
	size_t Ntot = pars->Ntot;
	size_t Ntot2 = 2 * Ntot;
	size_t Nt = pars->t.Nt;
	double dt = pars->t.dt;

	double a1, a2, a3, a4, tmp_a;
	double *k1, *k2, *k3, *k4, *tmp_k;
	k1    = fftw_malloc(Ntot2 * sizeof *k1);
	k2    = fftw_malloc(Ntot2 * sizeof *k2);
	k3    = fftw_malloc(Ntot2 * sizeof *k3);
	k4    = fftw_malloc(Ntot2 * sizeof *k4);
	tmp_k = fftw_malloc(Ntot2 * sizeof *tmp_k);

	if (!(k1 && k2 && k3 && k4 && tmp_k))
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}

	RUNTIME_INFO(puts("Starting RK4 time evolution with:"));
	RUNTIME_INFO(printf("initial time: %f\n", pars->t.ti));
	RUNTIME_INFO(printf("final time: %f\n", pars->t.tf));
	RUNTIME_INFO(printf("time step dt: %f\n", dt));
	RUNTIME_INFO(printf("number of steps: %zu\n", Nt));
	RUNTIME_INFO(puts("Using DFT (fftw3) for spatial derivatives."));

	#ifdef ENABLE_FFT_FILTER
		RUNTIME_INFO(puts("Frequency cutoff filtering enabled."));
	#else
		RUNTIME_INFO(puts("Filtering disabled."));
	#endif

	double start = get_wall_time();

	size_t write_count_field = 0, write_count_powspec = 0;

	for (size_t nt = 0; nt < Nt - 1; ++nt)
	{
		#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 1;
		#endif
		if (nt % pars->file.skip_powspec == 0 &&
			write_count_powspec < pars->file.num_powspec)
		{
			evo_flags.write_pow_spec = 1;
			write_count_powspec += 1;
		}
		if (nt % pars->file.skip_field == 0 &&
			write_count_field < pars->file.num_field)
		{
			file_append_by_name_1d(field, Ntot, 1, pars->file.name_field);
			write_count_field += 1;
			RUNTIME_INFO(printf("Writing field to disc at t = %f \n", nt * dt));
		}

		// k1 & a1
		a1 = mk_velocities(field, frw_a[nt], k1, pars);
		#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 0;
		#endif
		evo_flags.write_pow_spec = 0;

		// k2 & a2
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k1[i] / 2.0;
		}
		tmp_a = frw_a[nt] + dt * a1 / 2.0;
		a2 = mk_velocities(tmp_k, tmp_a, k2, pars);

		// k3 & a3
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k2[i] / 2.0;
		}
		tmp_a = frw_a[nt] + dt * a2 / 2.0;
		a3 = mk_velocities(tmp_k, tmp_a, k3, pars);

		// k4 & a4
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k3[i];
		}
		tmp_a = frw_a[nt] + dt * a3;
		a4 = mk_velocities(tmp_k, tmp_a, k4, pars);

		rho[nt] = mk_rho(field, frw_a[nt], pars);

		// perform one time step for the field and a
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			field[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}

		frw_a[nt + 1] = frw_a[nt] + dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;

		#ifdef CHECK_FOR_NAN
			#pragma omp parallel for
			for (size_t i = 0; i < Ntot2; ++i)
			{
				if (isnan(field[i]))
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
	rho[Nt - 1] = mk_rho(field, frw_a[Nt - 1], pars);

	// write out last time slice
	if (pars->file.mode_field == 1)
	{
		if (fabs(pars->t.tf - (Nt - 1) * dt) > 1e-10)
		{
			RUNTIME_INFO(fputs("The time of the last step does not coincide "
								"with the specified final time.", stderr));
		}
		RUNTIME_INFO(printf("Writing field at tf = %f \n", (Nt - 1) * dt));
		file_append_by_name_1d(field, Ntot, 1, pars->file.name_field);
	}

	double end = get_wall_time();

	fftw_free(k1);
	fftw_free(k2);
	fftw_free(k3);
	fftw_free(k4);
	fftw_free(tmp_k);

	double secs = end - start;
	RUNTIME_INFO(printf("Finished time evolution in: %f seconds.\n\n", secs));
}

/*
compute the right hand side of the pde, ie the first order temporal derivatives
*/
double mk_velocities(double *f, double a, double *result, parameters_t *pars) {
	size_t Ntot = pars->Ntot;
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

/*
A selection of potentials one can try, make sure to set the corresponding
potential_prime, the derivative is not computed automatically yet
TODO: change that?
*/
inline double potential(double f){
	double lambda = 100.0;
	return LAMBDA / (1.0 + exp(-lambda * f));

	// double theta, dtheta;
	// if (f > 5.0)
	// {
	// 	theta = 1.0;
	// }
	// else if (f < -5.0)
	// {
	// 	theta = 0.0;
	// }
	// else
	// {
	// 	theta = 1.0 / (1.0 + exp(- 100.0 * f));
	// }

	// return MASS * MASS * f * f / 2.0;

	// return MASS * MASS * f * f / 2.0 + COUPLING * f * f * f * f / 24.0;

	// return 0.0;
}

inline double potential_prime(double f) {
	double lambda = 100.0;
	double tmp = exp(lambda * f);
	return LAMBDA * lambda * tmp / ((1.0 + tmp) * (1.0 + tmp));

	// return MASS * MASS * f;

	// return MASS * MASS * f + COUPLING * f * f * f / 6.0;

	// return 20.0 * tanh(pow(f, 50));

	// return 0.0;
}

/*
compute average 00 component of stress energy
*/
double mk_rho(double *f, double a, parameters_t *pars) {
	size_t Ntot = pars->Ntot;

	double T00 = 0.0;

	mk_gradient_squared_and_laplacian(f, dtmp_grad2, dtmp_lap, pars);

	double ft, grad2_a;
	#pragma omp parallel for default(shared) private(ft, grad2_a) reduction(+:T00)
	for (size_t i = 0; i < Ntot; ++i)
	{
		ft = f[Ntot + i];
		grad2_a = dtmp_grad2[i] / (a * a);
		T00 += (ft * ft + grad2_a) / 2. + potential(f[i]);
	}
	return T00 / Ntot;
}