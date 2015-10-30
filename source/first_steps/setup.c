#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include "setup.h"
#include "main.h"

void allocate_and_initialize_all() {
    allocate_external(NOGP);
    mk_fourier_spectral_operators(NOGP, LOW_BND, UP_BND);
    mk_initial_conditions(NOGP, phi_init, dphi_init);
}

/*
allocate memory and initialize all external variables
*/
void allocate_external(size_t N) {
    size_t N2 = N * N;

    //running parameter for time evolution and number of timesteps
    nt = 0;
    TS = ceil( ((double)(TF - TI)) / DT );

    //grid points and spectral operators
    x  = malloc(N  * sizeof *x);
    D1 = malloc(N2 * sizeof *D1);
    D2 = malloc(N2 * sizeof *D2);

    //solutions for phi (we are saving each timestep: N * TS space)
    phi   = calloc(N * TS, sizeof *phi);
    phiD1 = calloc(N * TS, sizeof *phiD1);
    phiD2 = calloc(N * TS, sizeof *phiD2);

    //solution for dphi
    dphi   = calloc(N * TS, sizeof *dphi);
    dphiD1 = calloc(N * TS, sizeof *dphiD1);
    dphiD2 = calloc(N * TS, sizeof *dphiD2);
}

/*
	make the N fourier spectral gridpoints for the interval [low_bnd, up_bnd],
	and the spectral operators (N \times N matrices) for first and second order
	derivatives; x, D1 and D2 need to be initialized and memory needs to be
	allocated before calling this function
	see mathematica notebook spectral_operators for further info
*/
void mk_fourier_spectral_operators(size_t N, double low_bnd, double up_bnd) {
	if (N <= 0 || N%2 == 0)
	{
		fputs("Need an odd number of gridpoints for a fourier grid", stderr);
	}

	// Grid points
	for (size_t i = 0; i < N; ++i)
	{
		x[i] = low_bnd + (up_bnd - low_bnd)*i/N;
	}

	//D1
	for (size_t i = 0; i < N; ++i)
	{
		size_t offset = i*N;
		double tmpD1 = 0.0;
		for (size_t j = 0; j < N; ++j)
		{
			if (i != j)
			{
				int diff = (int)j - (int)i;
				tmpD1 = PI / (sin((diff*PI)/(int)N) * (up_bnd-low_bnd));
				if ( (i+j)%2 != 0 )
				{
					tmpD1 *= -1.0;
				}
			}
			else
			{
				tmpD1 = 0.0;
			}
			D1[offset+j] = tmpD1;
		}
	}

	//D2
	sq_matrix(D1, D2, N);

	// Console output for debugging
	#ifdef PRINT_SPECTRAL_OPERATORS
		printf("\nx\n");
        print_vector(x, N);
        printf("\n\n\nD1\n");
        print_matrix(D1, N);
        printf("\n\n\nD2\n");
        print_matrix(D2, N);
	#endif
}

/*
setup initial conditions for phi
*/
void mk_initial_conditions(size_t N, double (*phi_init)(double),
								double (*dphi_init)(double)) {
    for (size_t i = 0; i < N; ++i)
    {
        phi[i] = phi_init(x[i]);
        dphi[i] = dphi_init(x[i]);
    }
}

/*
example test functions for initial conditions
when using fourier grid points, those need to be periodic in the spatial domain
specified by LOW_BND and UP_BND in main.h
*/
double phi_init(double x) {
	return sin(x);
}

double dphi_init(double x) {
	return -sin(x);
}

/*
free all allocated memory
*/
void free_all_external() {
	free(x);
	free(D1);
	free(D2);

	free(phi);
	free(phiD1);
	free(phiD2);
	free(dphi);
	free(dphiD1);
	free(dphiD2);
}

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