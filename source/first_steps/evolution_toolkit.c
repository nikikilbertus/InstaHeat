#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"
#include "filehandling.h"

evolution_flags_t evo_flags = {.filter = 0, .write_pow_spec = 0};

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

	#ifdef ENABLE_FFT_FILTER
		if (evo_flags.filter == 1)
		{
			fft_apply_filter(cfftw_tmp, pars);
		}
	#endif

	#ifdef WRITE_OUT_POWER_SPECTRUM
		if (evo_flags.write_pow_spec == 1)
		{
			mk_and_write_power_spectrum(cfftw_tmp, pars);
		}
	#endif

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
	// 		pars->cutoff_fraction =
	// 					fmax(0.0, (double)(nc - i - 3) / (double)nc);
	// 		#ifdef DEBUG
	// 			printf("cutoff was set to: %f\n", pars->cutoff_fraction);
	// 		#endif
	// 		break;
	// 	}
	// }
}

void mk_and_write_power_spectrum(fftw_complex *in, parameters_t *pars) {
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t ncz = Nz / 2 + 1;

	double Lx = pars->x.b - pars->x.a;
	double Ly = pars->y.b - pars->y.a;
	double Lz = pars->z.b - pars->z.a;

	double prefac2 = 4. * PI *PI;
	double k_x2 = prefac2 / (Lx * Lx);
	double k_y2 = prefac2 / (Ly * Ly);
	double k_z2 = prefac2 / (Lz * Lz);

	double k_min2 = 0.0;
	double k_max2 = k_x2 * (Nx - 1)  * (Nx - 1) +
					k_y2 * (Ny - 1)  * (Ny - 1) +
					k_z2 * (ncz - 1) * (ncz - 1);

	double dk2 = (k_max2 - k_min2) / pars->pow_spec_shells;
	double k2_tmp = 0.0;
	double pow_tmp = 0.0;

	for (size_t i = 0; i < pars->pow_spec_shells; ++i)
	{
		pow_spec[i] = 0.0;
	}

	size_t osx, osy, idx;
	for (size_t i = 0; i < Nx; ++i)
	{
		osx = i * Ny * ncz;
		for (size_t j = 0; j < Ny; ++j)
		{
			osy = osx + j * ncz;
			for (size_t k = 0; k < ncz; ++k)
			{
				k2_tmp = k_z2 * k * k;
				if (i > Nx / 2)
				{
					k2_tmp += k_x2 * (Nx - i) * (Nx - i);
				}
				else
				{
					k2_tmp += k_x2 * i * i;
				}

				if (j > Ny / 2)
				{
					k2_tmp += k_x2 * (Ny - j) * (Ny - j);
				}
				else
				{
					k2_tmp += k_y2 * j * j;
				}
				pow_tmp = cabs(in[osy + k]);
				idx = (int)(k2_tmp / dk2 - 1e-10);
				pow_spec[idx] += pow_tmp * pow_tmp / (Ntot * Ntot);
			}
		}
	}
	file_append_by_name_1d(pow_spec, pars->pow_spec_shells, 1, POW_SPEC_NAME, 0);
}

void fft_apply_filter(fftw_complex *inout, parameters_t *pars) {
	double cutoff_fraction = pars->cutoff_fraction;
	if (cutoff_fraction < 0.0 || cutoff_fraction > 1.0)
	{
		fputs("cutoff_fraction must be between 0 and 1.", stderr);
    	exit(EXIT_FAILURE);
	}
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;
	size_t ncz = Nz / 2 + 1;
	size_t nxmax = (size_t) fmin(Nx, ceil(Nx * (1.0 - cutoff_fraction)));
	size_t nymax = (size_t) fmin(Ny, ceil(Ny * (1.0 - cutoff_fraction)));
	size_t nzmax = (size_t) fmin(ncz, ceil(ncz * (1.0 - cutoff_fraction)));

	if (Nx == nxmax && Ny == nymax && ncz == nzmax)
	{
		fputs("Warning: filtering is enabled (ENABLE_FFT_FILTER "
			  "but nothing will be cut off -> Overhead!\n", stderr);
	}

	size_t osx, osy;
	for (size_t i = nxmax; i < Nx; ++i)
	{
		osx = i * Ny * ncz;
		for (size_t j = nymax; j < Ny; ++j)
		{
			osy = osx + j * ncz;
			for (size_t k = nzmax; k < ncz; ++k)
			{
				// might want to do filter window later
				inout[osy + k] = 0.0;
			}
		}
	}
}

void mk_filter_window(double *inout, size_t cutoffindex, size_t windowlength) {

	for (size_t i = 0; i < cutoffindex; ++i)
	{
		inout[i] = filter_window_function((double)i / cutoffindex);
	}
	for (size_t i = cutoffindex; i < windowlength; ++i)
	{
		inout[i] = 0.0;
	}
}

inline double filter_window_function(double x) {
	return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. );
	// return exp(1. + 1. / ( pow(x, 8) - 1. ));
	// return 0.5 * ( 1. + cos( pow(x, 8) * PI ) );
	// return 0.0;
}