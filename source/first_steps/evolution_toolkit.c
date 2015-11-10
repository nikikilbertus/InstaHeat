#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
// #include <Accelerate/Accelerate.h>
#include "evolution_toolkit.h"
#include "main.h"

void fft_D1(double *in, double *out, size_t N) {
	size_t nc = N / 2 + 1;
	fftw_complex *tmp;
	fftw_plan p_fw;
	fftw_plan p_bw;

	tmp = fftw_malloc(nc * sizeof *tmp);
	if (!tmp)
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}

	p_fw = fftw_plan_dft_r2c_1d(N, in, tmp, FFTW_ESTIMATE);
	fftw_execute(p_fw);

	double L = pars.b - pars.a;
	double complex factor = 2. * PI * I / L;
	for (size_t i = 0; i < nc; ++i)
	{
		tmp[i] *= factor * i / N;
	}

	p_bw = fftw_plan_dft_c2r_1d(N, tmp, out, FFTW_ESTIMATE);
	fftw_execute(p_bw);

	fftw_destroy_plan(p_fw);
	fftw_destroy_plan(p_bw);
	fftw_free(tmp);
}

void fft_D2(double *in, double *out, size_t N) {
	size_t nc = N / 2 + 1;
	fftw_complex *tmp;
	fftw_plan p_fw;
	fftw_plan p_bw;

	tmp = fftw_malloc(nc * sizeof *tmp);
	if (!tmp)
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}

	p_fw = fftw_plan_dft_r2c_1d(N, in, tmp, FFTW_ESTIMATE);
	fftw_execute(p_fw);

	#if defined(ENABLE_FFT_FILTER) && defined(ENABLE_ADAPTIVE_FILTER)
	{
		double *power_spec = malloc(nc * sizeof *power_spec);
		double absval;
		double cpsd = CPSD_FRACTION, threshold;
		double total = 0.0, fraction = 0.0;

		for (size_t i = 0; i < nc; ++i)
		{
			absval = cabs(tmp[i]);
			power_spec[i] = absval * absval;
			total += power_spec[i];
			// printf("%zu: %f\n", i, absval);
		}
		// puts("\n\n");

		threshold = cpsd * total;
		for (size_t i = 0; i < nc; ++i)
		{
			fraction += power_spec[i];
			if (fraction >= threshold)
			{
				pars.cutoff_fraction =
							fmax(0.0, (double)(nc - i - 3) / (double)nc);
				printf("cutoff was set to: %f\n", pars.cutoff_fraction);
				break;
			}
		}
		free(power_spec);
	}
	#endif

	double L = pars.b - pars.a;
	double factor = - 4. * PI * PI / (L*L);
	for (size_t i = 0; i < nc; ++i)
	{
		tmp[i] *= (factor * i * i) / N;
	}

	p_bw = fftw_plan_dft_c2r_1d(N, tmp, out, FFTW_ESTIMATE);
	fftw_execute(p_bw);

	fftw_destroy_plan(p_fw);
	fftw_destroy_plan(p_bw);
	fftw_free(tmp);
}

void fft_apply_filter(double *in, size_t N) {
	double cutoff_fraction = pars.cutoff_fraction;
	if (cutoff_fraction < 0.0 || cutoff_fraction > 1.0)
	{
		fputs("cutoff_fraction must be between 0 and 1.", stderr);
    	exit(EXIT_FAILURE);
	}
	size_t nc = N / 2 + 1;
	size_t nmax = (size_t) fmin(nc, ceil(nc * (1.0 - cutoff_fraction)));

	// if there is nothing to cut off
	if (nc == nmax)
	{
		fputs("Warning: filtering is enabled (ENABLE_FFT_FILTER "
			  "but the cutoff fraction is 0. FFT is performed"
			  "but nothing will be filtered -> Overhead!\n"
			  , stderr);
	}

	double *window;
	fftw_complex *tmp;
	fftw_plan p_fw;
	fftw_plan p_bw;

	window = malloc(nc * sizeof *window);
	tmp = fftw_malloc(nc * sizeof *tmp);
	if (!tmp || !window)
	{
		fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
	}

	mk_filter_window(window, nmax, nc);

	// we filter the field and its temporal derivative
	for (size_t i = 0; i <= N; i += N)
	{
		p_fw = fftw_plan_dft_r2c_1d(N, in + i, tmp, FFTW_ESTIMATE);
		p_bw = fftw_plan_dft_c2r_1d(N, tmp, in + i, FFTW_ESTIMATE);

		fftw_execute(p_fw);

		for (size_t j = 0; j < nc; ++j)
		{
			tmp[j] *= window[j] / N;
		}

		fftw_execute(p_bw);
	}

	fftw_destroy_plan(p_fw);
	fftw_destroy_plan(p_bw);
	fftw_free(tmp);
	free(window);
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

#ifdef USE_ACCELERATE_FRAMEWORK
#ifdef SPECTRAL_OPERATOR_DERIVATIVE
void spectral_op_D2(double *in, double *out, size_t N) {
	matrix_vector(D2, in, out, N);
}
#endif

//****************************** matrix/vector operations
/*
matrix-matrix and matrix-vector function wrapper
*/
void matrix_matrix(double *matrixA, double *matrixB, double *result, size_t N) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0,
    				matrixA, N, matrixB, N, 0.0, result, N);
}

void sq_matrix(double *matrix, double *result, size_t N) {
    matrix_matrix(matrix, matrix, result, N);
}

void matrix_vector(double *matrix, double *vector, double *result, size_t N) {
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0,
    				matrix, N, vector, 1, 0.0, result, 1);
}
#endif

//****************************** printing functions
/*
for debugging mostly, watch out for colomn_major vs row_major in lapack routines
*/
void print_matrix(double *matrix, size_t N){
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            printf("%f\t", matrix[j*N+i]);
        }
        printf("\n");
    }
}

void print_vector(double *vector, size_t N) {
    for (size_t i = 0; i < N; i++)
    {
        printf("%f\n", vector[i]);
    }
}

/*
conversion from 2D to 1D indices
*/
size_t idx(size_t N, size_t row, size_t col) {
	return row * N + col;
}