#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
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

	double start =  get_wall_time();
	fftw_execute_dft_r2c(p_fw_3d, in, cfftw_tmp);
	double end =  get_wall_time();
	fftw_time_exe += end - start;


	if (evo_flags.write_pow_spec == 1)
	{
		mk_and_write_power_spectrum(cfftw_tmp, pars);
	}

	#ifdef ENABLE_FFT_FILTER
		if (evo_flags.filter == 1)
		{
			fft_apply_filter(cfftw_tmp, pars);
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
	#pragma omp parallel for private(osx, osy, id, k_sq)
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

	start =  get_wall_time();
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_x, dtmp_x);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_y, dtmp_y);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_z, dtmp_z);
	fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp,  laplacian);
	end =  get_wall_time();
	fftw_time_exe += end - start;

	double gx, gy, gz;
	#pragma omp parallel for private(gx, gy, gz)
	for (size_t i = 0; i < Ntot; ++i)
	{
		gx = dtmp_x[i];
		gy = dtmp_y[i];
		gz = dtmp_z[i];
		grad2[i] = gx * gx + gy * gy + gz * gz;
	}
}

void mk_and_write_power_spectrum(fftw_complex *in, parameters_t *pars) {
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;
	size_t Ntot = Nx * Ny * Nz;
	size_t ncx = Nx / 2 + 1;
	size_t ncy = Ny / 2 + 1;
	size_t ncz = Nz / 2 + 1;
	size_t bins = pars->file.bins_powspec;

	// todo[performance]: precompute bins only once and reuse
	double Lx = pars->x.b - pars->x.a;
	double Ly = pars->y.b - pars->y.a;
	double Lz = pars->z.b - pars->z.a;

	double prefac2 = 4. * PI *PI;
	double k_x2 = prefac2 / (Lx * Lx);
	double k_y2 = prefac2 / (Ly * Ly);
	double k_z2 = prefac2 / (Lz * Lz);

	double k_min2 = 0.0;
	double k_max2 = k_x2 * (Nx/2) * (Nx/2) +
					k_y2 * (Ny/2) * (Ny/2) +
					k_z2 * (Nz/2) * (Nz/2);

	double dk2 = (k_max2 - k_min2) / bins;
	double k2_tmp = 0.0;
	double pow2_tmp = 0.0;

	#pragma omp parallel for
	for (size_t i = 0; i < bins; ++i)
	{
		pow_spec[i] = 0.0;
	}

	size_t osx, osy, idx;
	// todo parallelize?
	for (size_t i = 0; i < Nx; ++i)
	{
		osx = i * Ny * ncz;
		for (size_t j = 0; j < Ny; ++j)
		{
			osy = osx + j * ncz;
			for (size_t k = 0; k < ncz; ++k)
			{
				if (k == 0 || 2 * k == Nz)
				{
					pow2_tmp = in[osy + k] * conj(in[osy + k]);
				}
				else
				{
					pow2_tmp = 2.0 * in[osy + k] * conj(in[osy + k]);
				}
				if (pow2_tmp > 0.0)
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
						k2_tmp += k_y2 * (Ny - j) * (Ny - j);
					}
					else
					{
						k2_tmp += k_y2 * j * j;
					}
					idx = (int)(k2_tmp / dk2 - 1e-7);
					pow_spec[idx] += pow2_tmp / Ntot;
				}
			}
		}
	}
	file_append_by_name_1d(pow_spec, bins, 1, pars->file.name_powspec);
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
	#pragma omp parallel for private(osx, osy)
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

	#pragma omp parallel for
	for (size_t i = 0; i < cutoffindex; ++i)
	{
		inout[i] = filter_window_function((double)i / cutoffindex);
	}
	#pragma omp parallel for
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