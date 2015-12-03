#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"

void run_rk4(parameters_t *pars) {
	size_t Ntot = pars->Ntot;
	size_t Ntot2 = 2 * Ntot;
	size_t Nt = pars->t.Nt;
	double dt = pars->t.dt;
	double t = pars->t.t;

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

	for (size_t nt = 0; nt < Nt - 1; ++nt)
	{
		#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 1;
		#endif
		if (nt % pars->file.skip == 0)
		{
			evo_flags.compute_pow_spec = 1;
		}

		// k1 & a1
		a1 = mk_velocities(t, field, f_a, k1, pars);
		#ifdef ENABLE_FFT_FILTER
			evo_flags.filter = 0;
		#endif
		evo_flags.compute_pow_spec = 0;

		// k2 & a2
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k1[i] / 2.0;
		}
		tmp_a = f_a + dt * a1 / 2.0;
		a2 = mk_velocities(t + dt / 2.0, tmp_k, tmp_a, k2, pars);

		// k3 & a3
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k2[i] / 2.0;
		}
		tmp_a = f_a + dt * a2 / 2.0;
		a3 = mk_velocities(t + dt / 2.0, tmp_k, tmp_a, k3, pars);

		// k4 & a4
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			tmp_k[i] = field[i] + dt * k3[i];
		}
		tmp_a = f_a + dt * a3;
		a4 = mk_velocities(t + dt, tmp_k, tmp_a, k4, pars);

		rho = mk_rho(field, f_a, pars);

		if (nt % pars->file.skip == 0)
		{
			save();
		}

		// perform one time step for the field and a
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			field[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}
		f_a += dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0;

		t += dt;
		pars->t.t = t;

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
			if (isnan(f_a))
			{
				fprintf(stderr,
						"A nan value was discovered in timestep: %zu \n", nt);
			}
		#endif
	}

	// make sure to write out last time slice
	rho = mk_rho(field, f_a, pars);
	save();

	// info about last timeslice
	if (fabs(pars->t.tf - pars->t.t) > 1e-10)
	{
		RUNTIME_INFO(fputs("The time of the last step does not coincide "
							"with the specified final time.", stderr));
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