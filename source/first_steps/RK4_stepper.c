#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include "RK4_stepper.h"
#include "main.h"
#include "setup.h"

void run_RK4_stepper(double dt, size_t N) {
	size_t M = 2 * N;
	double *k1, *k2, *k3, *k4, *tmp;
	k1 = malloc(M * sizeof *k1);
	k2 = malloc(M * sizeof *k2);
	k3 = malloc(M * sizeof *k3);
	k4 = malloc(M * sizeof *k4);
	tmp = malloc(M * sizeof *tmp);
	size_t os, new_os;

	for (nt = 0; nt < TS; ++nt)
	{
		os = nt * N;

		// k1
		get_field_velocity(phi + os, dphi + os, k1, N);

		// k2
		for (size_t i = 0; i < N; ++i)
		{
			tmp[i] = dphi[os+i] + dt * k1[i] / 2.0;
		}
		for (size_t i = N; i < M; ++i)
		{
			tmp[i] = phi[os+i] + dt * k1[i] / 2.0;
		}
		get_field_velocity(tmp, tmp + N, k2, N);

		// k3
		for (size_t i = 0; i < N; ++i)
		{
			tmp[i] = dphi[os+i] + dt * k2[i] / 2.0;
		}
		for (size_t i = N; i < M; ++i)
		{
			tmp[i] = phi[os+i] + dt * k2[i] / 2.0;
		}
		get_field_velocity(tmp, tmp + N, k3, N);

		// k4
		for (size_t i = 0; i < N; ++i)
		{
			tmp[i] = dphi[os+i] + dt * k3[i];
		}
		for (size_t i = N; i < M; ++i)
		{
			tmp[i] = phi[os+i] + dt * k3[i];
		}
		get_field_velocity(tmp, tmp + N, k4, N);

		// perform one time step
		new_os = os + N;
		for (size_t i = 0; i < N; ++i)
		{
			phi[new_os+i] = phi[os+i] +
						dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}
		for (size_t i = N; i < M; ++i)
		{
			dphi[os+i] = dphi[os+i-N] +
						dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}
	}

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(tmp);
}

void get_field_velocity(double* f, double *df, double *result, size_t N) {
	for (int i = 0; i < N; ++i)
	{
		result[i] = df[i];
	}
	matrix_vector(D2, f, result + N, N);
}