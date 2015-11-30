#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "RK4_stepper.h"
#include "main.h"
#include "evolution_toolkit.h"
#include "filehandling.h"
#include "dopri853_stepper.h"

dopri853_control_t dp;
dopri853_values_t dpv;

void allocate_dopri853_values() {

}

void perform_step(const double dt_try) {
	double dt = dt_try;
	for ( ; ; )
	{
		try_step();
		double err = error(dt);
		if (success(err, dt))
		{
			break;
		}
		if (fabs(dt) <= fabs(dp.t) * dp.eps)
		{
			fputs("Stepsize underflow", stderr);
			exit(EXIT_FAILURE);
		}
		dpv.a1 = mk_velocities(t + dt, field_new, frw_a_new, dpv.k1);
		if (dense)
		{
			prepare_dense_output(dt);
		}
		#pragma omp parallel for
		for (size_t i = 0; i < Ntot2; ++i)
		{
			field[i] = field_new[i];
		}
		dp.t_old = dp.t;
		dp.t += (dp.dt_did = dp.dt);
	}
}

double error(const double dt) {
	double err = 0.0, err2 = 0.0, sk, deno;
	for (size_t i = 0; i < Ntot2; ++i)
	{
		sk = pd.a_tol + pd.r_tol * MAX(fabs(field[i]), fabs(field_new[i]));
		err2 += sqrt(yerr[i] / sk);
		err  += sqrt(yerr2[i] / sk);
	}
	deno = err + 0.01 * err2;
	if (deno <= 0.0)
	{
		deno = 1.0;
	}
	return dt * err * sqrt(1.0 / (Ntot2 * deno));
}

void try_step() {
	size_t i;
	double t = dp.t;
	double dt = pars->t.dt;
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t Ntot2 = Ntot2;
	// ------------ 1 ------------
	dpv.a1 = mk_velocities(t, field, frw_a, dpv.k1, pars);

	// ------------ 2 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt * dpc.a21 * dpv.k1[i];
	}
	a_tmp = frw_a + dt * dpc.a21 * dpv.a1;
	dpv.a2 = mk_velocities(t + dpc.c2 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k2, pars);

	// ------------ 3 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a31 * dpv.k1[i] + dpc.a32 * dpv.k2[i]);
	}
	a_tmp = frw_a + dt * (dpc.a31 * dpv.a1 + dpc.a32 * dpv.a2);
	dpv.a3 = mk_velocities(t + dpc.c3 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k3, pars);

	// ------------ 4 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a41 * dpv.k1[i] + dpc.a43 * dpv.k3[i]);
	}
	a_tmp = frw_a + dt * (dpc.a41 * dpv.a1 + dpc.a43 * dpv.a3);
	dpv.a4 = mk_velocities(t + dpc.c4 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k4, pars);

	// ------------ 5 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a51 * dpv.k1[i] + dpc.a53 * dpv.k3[i] + dpc.a54 * dpv.k4[i]);
	}
	a_tmp = frw_a + dt * (dpc.a51 * dpv.a1 + dpc.a53 * dpv.a3 + dpc.a54 * dpv.a4);
	dpv.a5 = mk_velocities(t + dpc.c5 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k5, pars);

	// ------------ 6 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a61 * dpv.k1[i] + dpc.a64 * dpv.k4[i] + dpc.a65 * dpv.k5[i]);
	}
	a_tmp = frw_a + dt * (dpc.a61 * dpv.a1 + dpc.a64 * dpv.a4 + dpc.a65 * dpv.a5);
	dpv.a6 = mk_velocities(t + dpc.c6 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k6, pars);

	// ------------ 7 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a71 * dpv.k1[i] + dpc.a74 * dpv.k4[i] + dpc.a75 * dpv.k5[i] +
			 dpc.a76 * dpv.k6[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a71 * dpv.a1 + dpc.a74 * dpv.a4 + dpc.a75 * dpv.a5 +
			 dpc.a76 * dpv.a6);
	dpv.a7 = mk_velocities(t + dpc.c7 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k7, pars);

	// ------------ 8 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a81 * dpv.k1[i] + dpc.a84 * dpv.k4[i] + dpc.a85 * dpv.k5[i] +
			 dpc.a86 * dpv.k6[i] + dpc.a87 * dpv.k7[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a81 * dpv.a1 + dpc.a84 * dpv.a4 + dpc.a85 * dpv.a5 +
			 dpc.a86 * dpv.a6 + dpc.a87 * dpv.a7);
	dpv.a8 = mk_velocities(t + dpc.c8 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k8, pars);

	// ------------ 9 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a91 * dpv.k1[i] + dpc.a94 * dpv.k4[i] + dpc.a95 * dpv.k5[i] +
			 dpc.a96 * dpv.k6[i] + dpc.a97 * dpv.k7[i] + dpc.a98 * dpv.k8[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a91 * dpv.a1 + dpc.a94 * dpv.a4 + dpc.a95 * dpv.a5 +
			 dpc.a96 * dpv.a6 + dpc.a97 * dpv.a7 + dpc.a98 * dpv.a8);
	dpv.a9 = mk_velocities(t + dpc.c9 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k9, pars);

	// ------------ 10 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a101 * dpv.k1[i] + dpc.a104 * dpv.k4[i] + dpc.a105 * dpv.k5[i] +
			 dpc.a106 * dpv.k6[i] + dpc.a107 * dpv.k7[i] + dpc.a108 * dpv.k8[i] +
			 dpc.a109 * dpv.k9[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a101 * dpv.a1 + dpc.a104 * dpv.a4 + dpc.a105 * dpv.a5 +
			 dpc.a106 * dpv.a6 + dpc.a107 * dpv.a7 + dpc.a108 * dpv.a8 +
			 dpc.a109 * dpv.a9);
	dpv.a10 = mk_velocities(t + dpc.c10 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k10, pars);

	// ------------ 11 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a111 * dpv.k1[i] + dpc.a114 * dpv.k4[i] + dpc.a115 * dpv.k5[i] +
			 dpc.a116 * dpv.k6[i] + dpc.a117 * dpv.k7[i] + dpc.a118 * dpv.k8[i] +
			 dpc.a119 * dpv.k9[i] + dpc.a1110 * dpv.k10[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a111 * dpv.a1 + dpc.a114 * dpv.a4 + dpc.a115 * dpv.a5 +
			 dpc.a116 * dpv.a6 + dpc.a117 * dpv.a7 + dpc.a118 * dpv.a8 +
			 dpc.a119 * dpv.a9 + dpc.a1110 * dpv.a10);
	dpv.a2 = mk_velocities(t + dpc.c11 * dt, dpv.k_tmp, dpv.a_tmp, dpv.k2, pars);

	// ------------ new dt ------------
	double tpdt = t + dt;

	// ------------ 12 ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k_tmp[i] = field[i] + dt *
			(dpc.a121 * dpv.k1[i] + dpc.a124 * dpv.k4[i] + dpc.a125 * dpv.k5[i] +
			 dpc.a126 * dpv.k6[i] + dpc.a127 * dpv.k7[i] + dpc.a128 * dpv.k8[i] +
			 dpc.a129 * dpv.k9[i] + dpc.a1210 * dpv.k10[i] + dpc.a1211 * dpv.k2[i]);
	}
	a_tmp = frw_a + dt *
			(dpc.a121 * dpv.a1 + dpc.a124 * dpv.a4 + dpc.a125 * dpv.a5 +
			 dpc.a126 * dpv.a6 + dpc.a127 * dpv.a7 + dpc.a128 * dpv.a8 +
			 dpc.a129 * dpv.a9 + dpc.a1210 * dpv.a10 + dpc.a1211 * dpv.a2);
	dpv.a3 = mk_velocities(tpdt, dpv.k_tmp, dpv.a_tmp, dpv.k3, pars);

	// ------------ step ahead ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k4[i] = dpc.b1 * dpv.k1[i] + dpc.b6 * dpv.k6[i] + dpc.b7 * dpv.k7[i] +
					dpc.b8 * dpv.k8[i] + dpc.b9 * dpv.k9[i] + dpc.b10 * dpv.k10[i] +
					dpc.b11 * dpv.k2[i] + dpc.b12 * dpv.k3[i];

		field_new = field[i] + dt * dpv.k4[i];
	}
	dpv.a4 = dpc.b1 * dpv.a1 + dpc.b6 * dpv.a6 + dpc.b7 * dpv.a7 +
			 dpc.b8 * dpv.a8 + dpc.b9 * dpv.a9 + dpc.b10 * dpv.a10 +
			 dpc.b11 * dpv.a2 + dpc.b12 * dpv.a3;
	frw_a_new = frw_a + dt * dpv.a4;

	// ------------ step ahead ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.k4[i] = dpc.b1 * dpv.k1[i] + dpc.b6 * dpv.k6[i] + dpc.b7 * dpv.k7[i] +
					dpc.b8 * dpv.k8[i] + dpc.b9 * dpv.k9[i] + dpc.b10 * dpv.k10[i] +
					dpc.b11 * dpv.k2[i] + dpc.b12 * dpv.k3[i];

		field_new = field[i] + dt * dpv.k4[i];
	}
	dpv.a4 = dpc.b1 * dpv.a1 + dpc.b6 * dpv.a6 + dpc.b7 * dpv.a7 +
			 dpc.b8 * dpv.a8 + dpc.b9 * dpv.a9 + dpc.b10 * dpv.a10 +
			 dpc.b11 * dpv.a2 + dpc.b12 * dpv.a3;
	frw_a_new = frw_a + dt * dpv.a4;

	// ------------ error estimates ------------
	#pragma omp parallel for
	for (i = 0; i < Ntot2; ++i)
	{
		dpv.yerr[i] = dpv.k4[i] - dpc.bhh1 * dpv.k1[i] - dpc.bhh2 * dpv.k9[i] -
					  dpc.bhh3 * dpv.k3[i];
		dpv.yerr2[i] = pdc.er1 * dpv.k1[i] + pdc.er6 * dpv.k6[i] +
					   pdc.er7 * dpv.k7[i] + pdc.er8 * dpv.k8[i] +
					   pdc.er9 * dpv.k9[i] + pdc.er10 * dpv.k10[i] +
					   pdc.er11 * dpv.k2[i] + pdc.er12 * dpv.k3[i];
	}
}