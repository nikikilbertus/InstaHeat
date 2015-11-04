#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include "setup.h"
#include "main.h"

void allocate_and_initialize_all(parameters_t *pars) {
	initialize_parameters(pars);
	DEBUG(puts("Initialized parameters.\n"));
	size_t Nx = pars->Nx;
	size_t Nt = pars->Nt;
    double a = pars->a;
	double b = pars->b;
    allocate_external(Nx, Nt);
    DEBUG(puts("Allocated memory for all external variables.\n"));
    mk_fourier_spectral_operators(Nx, a, b);
    DEBUG(puts("Constructed fourier gridpoints (and spectral operators).\n"));
    mk_initial_conditions(Nx, phi_init, dphi_init);
    DEBUG(puts("Initialized the field and its temporal derivative.\n"));
}

void initialize_parameters(parameters_t *pars) {
	pars->Nx = GRIDPOINTS_SPATIAL;
	pars->dt = DELTA_T > 0 ? DELTA_T : pars->dt;
	pars->ti = INITIAL_TIME;
	pars->tf = FINAL_TIME;
	pars->a  = SPATIAL_LOWER_BOUND;
	pars->b  = SPATIAL_UPPER_BOUND;
	pars->Nt = ceil((pars->tf - pars->ti) / pars->dt);
	pars->nt = 0;
	pars->t  = 0.0;
}

/*
allocate memory and initialize all external variables
*/
void allocate_external(size_t Nx, size_t Nt) {
    size_t N2 = Nx * Nx;

    //grid points and spectral operators
    x  = malloc(Nx * sizeof *x);
    if (!x)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }
    #ifdef SPECTRAL_OPERATOR_DERIVATIVE
    D1 = malloc(N2 * sizeof *D1);
    D2 = malloc(N2 * sizeof *D2);
    if (!D1 || !D2)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }
    #endif

    //solutions for the field and the temporal derivative (we are saving each
    //timestep: 2 * Nx * Nt space)
    field = calloc(N2 * Nt, sizeof *field);
    if (!field)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    // solution for the scale parameter a: Nt space
    a = calloc(Nt, sizeof *a);
    if (!a)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }
}

/*
	make the N fourier spectral gridpoints for the interval [a, b],
	and the spectral operators (N \times N matrices) for first and second order
	derivatives; x, D1 and D2 need to be initialized and memory needs to be
	allocated before calling this function
	see mathematica notebook spectral_operators for further info
*/
void mk_fourier_spectral_operators(size_t N, double a, double b) {
	if (N <= 0 || N%2 == 0)
	{
		fputs("Need an odd number of gridpoints for a fourier grid", stderr);
	}

	// Grid points
	for (size_t i = 0; i < N; ++i)
	{
		x[i] = a + (b - a) * i / N;
	}

	#ifdef SPECTRAL_OPERATOR_DERIVATIVE
	//D1
	for (size_t i = 0; i < N; ++i)
	{
		size_t offset = i * N;
		double tmpD1 = 0.0;
		for (size_t j = 0; j < N; ++j)
		{
			if (i != j)
			{
				int diff = (int)j - (int)i;
				tmpD1 = PI / ( sin( (diff*PI) / (double)N ) * (b - a) );
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
	#endif

	// Console output for debugging
	#ifdef PRINT_SPECTRAL_OPERATORS
		puts("x");
        print_vector(x, N);
        #ifdef SPECTRAL_OPERATOR_DERIVATIVE
        puts("\nD1");
        print_matrix(D1, N);
        puts("\nD2");
        print_matrix(D2, N);
        #endif
        puts("\n");
	#endif
}

/*
setup initial conditions for the field
*/
void mk_initial_conditions(size_t N, double (*f_init)(double),
								double (*df_init)(double)) {
    for (size_t i = 0; i < N; ++i)
    {
        field[i] = f_init(x[i]);
        field[N+i] = df_init(x[i]);
    }

    a[0] = 1.0;

    // Console output for debugging
    #ifdef PRINT_INITIAL_CONDITIONS
    	puts("phi");
        print_vector(field, N);
        puts("\ndphi");
        print_vector(field + N, N);
        puts("\n");
    #endif
}

/*
example test functions for initial conditions
when using fourier grid points, those need to be periodic in the spatial domain
specified by LOW_BND and UP_BND in main.h
*/
double phi_init(double x) {
	// return sin(x);
	return tanh(pow(x, 8));
}

double dphi_init(double x) {
	// return -cos(x);
	return 0.0;
}

/*
free all allocated memory
*/
void free_all_external() {
	free(x);
	#ifdef SPECTRAL_OPERATOR_DERIVATIVE
	free(D1);
	free(D2);
	#endif

	free(field);
	DEBUG(puts("Memory from all external variables freed.\n"));
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