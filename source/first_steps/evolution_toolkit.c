#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"

void mk_gradient_squared_and_laplacian(double *in, double *grad2,
							double *laplacian, parameters_t *pars) {
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t ncz = Nz / 2 + 1;

	clock_t start = clock();
	fftw_execute_dft_r2c(p_fw_3d, in, cfftw_tmp);
	clock_t end = clock();
	fftw_time_exe += (double)(end - start) / CLOCKS_PER_SEC;

	double Lx = pars->x.b - pars->x.a;
	double Ly = pars->y.b - pars->y.a;
	double Lz = pars->z.b - pars->z.a;

	complex prefac = 2. * PI * I;
	complex factor_x = prefac / (Lx * Ntot);
	complex factor_y = prefac / (Ly * Ntot);
	complex factor_z = prefac / (Lz * Ntot);

	double prefac2 = -4. * PI *PI;
	double factor_x2 = prefac2 / (Lx * Lx * Ntot);
	double factor_y2 = prefac2 / (Ly * Ly * Ntot);
	double factor_z2 = prefac2 / (Lz * Lz * Ntot);
	double k_sq;

	size_t osx, osy, id;
	for (size_t i = 0; i < Nx; ++i)
	{
		osx = i * Ny * ncz;
		for (size_t j = 0; j < Ny; ++j)
		{
			osy = osx + j * ncz;
			for (size_t k = 0; k < ncz; ++k)
			{
				id = osy + k;

				// for lagrangian
				k_sq = factor_z2 * k * k;

				// x derivative
				if (i > Nx / 2)
				{
					cfftw_tmp_x[id] = cfftw_tmp[id] * factor_x
														* ((int)i - (int)Nx);
					k_sq += factor_x2 * (Nx - i) * (Nx - i);
				}
				else if (2 * i == Nx)
				{
					cfftw_tmp_x[id] = 0.0;
					k_sq += factor_x2 * i * i;
				}
				else
				{
					cfftw_tmp_x[id] = cfftw_tmp[id] * factor_x * i;
					k_sq += factor_x2 * i * i;
				}

				// y derivative
				if (j > Ny / 2)
				{
					cfftw_tmp_y[id] = cfftw_tmp[id] * factor_y
														* ((int)j - (int)Ny);
					k_sq += factor_y2 * (Ny - j) * (Ny - j);
				}
				else if (2 * j == Ny)
				{
					cfftw_tmp_y[id] = 0.0;
					k_sq += factor_y2 * j * j;
				}
				else
				{
					cfftw_tmp_y[id] = cfftw_tmp[id] * factor_y * j;
					k_sq += factor_y2 * j * j;
				}

				// z derivative
				if (k == ncz - 1)
				{
					cfftw_tmp_z[id] = 0.0;
				}
				else
				{
					cfftw_tmp_z[id] = cfftw_tmp[id] * factor_z * k;
				}

				// lagrangian
				cfftw_tmp[id] *= k_sq;
			}
		}
	}

	start = clock();
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_x, dtmp_x);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_y, dtmp_y);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_z, dtmp_z);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp,  laplacian);
	end = clock();
	fftw_time_exe += (double)(end - start) / CLOCKS_PER_SEC;

	double gx, gy, gz;
	for (size_t i = 0; i < Ntot; ++i)
	{
		gx = dtmp_x[i];
		gy = dtmp_y[i];
		gz = dtmp_z[i];
		grad2[i] = gx * gx + gy * gy + gz * gz;
	}
}

void set_adaptive_cutoff_fraction(size_t nc) {
	//TODO
	// double *power_spec = dtmp_x;
	// double absval;
	// double cpsd = CPSD_FRACTION, threshold;
	// double total = 0.0, fraction = 0.0;

	// for (size_t i = 0; i < nc; ++i)
	// {
	// 	absval = cabs(cfftw_tmp[i]);
	// 	power_spec[i] = absval * absval;
	// 	total += power_spec[i];
	// }

	// threshold = cpsd * total;
	// for (size_t i = 0; i < nc; ++i)
	// {
	// 	fraction += power_spec[i];
	// 	if (fraction >= threshold)
	// 	{
	// 		pars.cutoff_fraction =
	// 					fmax(0.0, (double)(nc - i - 3) / (double)nc);
	// 		#ifdef DEBUG
	// 			printf("cutoff was set to: %f\n", pars.cutoff_fraction);
	// 		#endif
	// 		break;
	// 	}
	// }
}

void fft_apply_filter(double *in, parameters_t *pars) {
	//TODO
	// double cutoff_fraction = pars.cutoff_fraction;
	// if (cutoff_fraction < 0.0 || cutoff_fraction > 1.0)
	// {
	// 	fputs("cutoff_fraction must be between 0 and 1.", stderr);
 //    	exit(EXIT_FAILURE);
	// }
	// size_t nc = N / 2 + 1;
	// size_t nmax = (size_t) fmin(nc, ceil(nc * (1.0 - cutoff_fraction)));

	// // warning if there is nothing to cut off
	// if (nc == nmax)
	// {
	// 	fputs("Warning: filtering is enabled (ENABLE_FFT_FILTER "
	// 		  "but the cutoff fraction is 0. FFT is performed"
	// 		  "but nothing will be filtered -> Overhead!\n"
	// 		  , stderr);
	// }

	// double *window = dtmp_x;
	// mk_filter_window(window, nmax, nc);

	// fftw_plan p_fw;
	// fftw_plan p_bw;

	// // we filter the field and its temporal derivative
	// for (size_t i = 0; i <= N; i += N)
	// {
	// 	p_fw = fftw_plan_dft_r2c_1d(N, in + i, cfftw_tmp, FFTW_DEFAULT_FLAG);
	// 	p_bw = fftw_plan_dft_c2r_1d(N, cfftw_tmp, in + i, FFTW_DEFAULT_FLAG);

	// 	fftw_execute(p_fw);

	// 	for (size_t j = 0; j < nc; ++j)
	// 	{
	// 		cfftw_tmp[j] *= window[j] / N;
	// 	}

	// 	fftw_execute(p_bw);
	// }

	// fftw_destroy_plan(p_fw);
	// fftw_destroy_plan(p_bw);
}

void mk_filter_window(double *out, size_t cutoffindex, size_t windowlength) {

	for (size_t i = 0; i < cutoffindex; ++i)
	{
		out[i] = filter_window_function((double)i / cutoffindex);
	}
	for (size_t i = cutoffindex; i < windowlength; ++i)
	{
		out[i] = 0.0;
	}
}

inline double filter_window_function(double x) {
	return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. );
	// return exp(1. + 1. / ( pow(x, 8) - 1. ));
	// return 0.5 * ( 1. + cos( pow(x, 8) * PI ) );
	// return 0.0;
}