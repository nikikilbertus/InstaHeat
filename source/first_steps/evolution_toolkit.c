#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"

void fft_D1(double *in, double *out, int direction, parameters_t *pars) {
	size_t Nx = pars->x.N;
	size_t Ny = pars->y.N;
	size_t Nz = pars->z.N;

	size_t Neffx, Neffy, Neffz, nc;
	double a, b;

	fftw_complex *fftw_tmp;
	int N, rank, howmany, idist, odist, istride, ostride;

	switch (direction) {
		case 1:
			N = Nx;
			nc = N / 2 + 1;
			Neffx = nc;
			Neffy = Ny;
			Neffz = Nz;
			howmany = Ny * Nz;
			idist   = 1;
			odist   = 1;
			istride = Ny * Nz;
			ostride = istride;
			fftw_tmp = cfftw_tmp_x;
			a = pars->x.a;
			b = pars->x.b;
			break;
		case 2:
			N = Ny;
			nc = N / 2 + 1;
			Neffx = Nx;
			Neffy = nc;
			Neffz = Nz;
			howmany = Nx * Nz;
			idist   = 1; // TODO: that one's a bitch, see notes
			odist   = 1; // TODO
			istride = Nz;
			ostride = istride;
			fftw_tmp = cfftw_tmp_y;
			a = pars->y.a;
			b = pars->y.b;
			break;
		case 3:
			N = Nz;
			nc = N / 2 + 1;
			Neffx = Nx;
			Neffy = Ny;
			Neffz = nc;
			howmany = Nx * Ny;
			idist   = Nz;
			odist   = nc;
			istride = 1;
			ostride = istride;
			fftw_tmp = cfftw_tmp_z;
			a = pars->z.a;
			b = pars->z.b;
			break;
		default:
			fputs("Direction must be 1 (x), 2 (y) or 3 (z).", stderr);
			exit(EXIT_FAILURE);
	}

	rank = 1;
	int n[] = {N};
	int *inembed = n;
	int *onembed = n;

	fftw_plan p_fw;
	fftw_plan p_bw;

	p_fw = fftw_plan_many_dft_r2c(rank, n, howmany, in, inembed, istride, idist,
								fftw_tmp, onembed, ostride, odist, FFTW_ESTIMATE);
	fftw_execute(p_fw);

	// todo correct differentiation
	double L = b - a;
	double complex factor = 2. * PI * I / (L * N);
	size_t osx, osy;
	for (size_t i = 0; i < Neffx; ++i)
	{
		osx = i * Neffy * Neffz;
		for (size_t j = 0; j < Neffy; ++j)
		{
			osy = osx + j * Neffz;
			for (size_t k = 0; k < Neffz; ++k)
			{
				switch (direction) {
					case 1:
						fftw_tmp[osy + k] *= (factor * i);
						break;
					case 2:
						fftw_tmp[osy + k] *= (factor * j);
						break;
					case 3:
						fftw_tmp[osy + k] *= (factor * k);
						break;
				}
			}
		}
	}

	p_bw = fftw_plan_many_dft_c2r(rank, n, howmany, fftw_tmp, inembed, istride,
							idist, out, onembed, ostride, odist, FFTW_ESTIMATE);
	fftw_execute(p_bw);

	fftw_destroy_plan(p_fw);
	fftw_destroy_plan(p_bw);
}

void fft_D2(double *in, double *out, int direction, parameters_t *pars) {
	//TODO
// 	size_t nc = N / 2 + 1;
// 	fftw_plan p_fw;
// 	fftw_plan p_bw;

// 	p_fw = fftw_plan_dft_r2c_1d(N, in, cfftw_tmp, FFTW_ESTIMATE);
// 	fftw_execute(p_fw);

// #if defined(ENABLE_FFT_FILTER) && defined(ENABLE_ADAPTIVE_FILTER)
// 	set_adaptive_cutoff_fraction(nc);
// #endif

// 	double L = pars.b - pars.a;
// 	double factor = - 4. * PI * PI / (L*L);
// 	for (size_t i = 0; i < nc; ++i)
// 	{
// 		cfftw_tmp[i] *= (factor * i * i) / N;
// 	}

// 	p_bw = fftw_plan_dft_c2r_1d(N, cfftw_tmp, out, FFTW_ESTIMATE);
// 	fftw_execute(p_bw);

// 	fftw_destroy_plan(p_fw);
// 	fftw_destroy_plan(p_bw);
}

void mk_gradient_squared(double *in, double *out, parameters_t *pars) {
	//TODO
}

void mk_laplacian(double *in, double *out, parameters_t *pars) {
	//TODO
}

void set_adaptive_cutoff_fraction(size_t nc) {
	//TODO
	// double *power_spec = dmisc_tmp;
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

	// double *window = dmisc_tmp;
	// mk_filter_window(window, nmax, nc);

	// fftw_plan p_fw;
	// fftw_plan p_bw;

	// // we filter the field and its temporal derivative
	// for (size_t i = 0; i <= N; i += N)
	// {
	// 	p_fw = fftw_plan_dft_r2c_1d(N, in + i, cfftw_tmp, FFTW_ESTIMATE);
	// 	p_bw = fftw_plan_dft_c2r_1d(N, cfftw_tmp, in + i, FFTW_ESTIMATE);

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

void mk_filter_window(double *window, size_t nmax, size_t nc) {

	for (size_t i = 0; i < nmax; ++i)
	{
		window[i] = filter_window_function((double)i / nmax);
	}
	for (size_t i = nmax; i < nc; ++i)
	{
		window[i] = 0.0;
	}
}

inline double filter_window_function(double x) {
	return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. );
	// return exp(1. + 1. / ( pow(x, 8) - 1. ));
	// return 0.5 * ( 1. + cos( pow(x, 8) * PI ) );
	// return 0.0;
}