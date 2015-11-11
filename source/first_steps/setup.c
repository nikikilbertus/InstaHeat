#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "setup.h"
#include "main.h"

void allocate_and_initialize_all(parameters_t *pars) {
	initialize_parameters(pars);
	RUNTIME_INFO(puts("Initialized parameters.\n"));
	size_t Nx = pars->Nx;
	size_t Nt = pars->Nt;
    double a = pars->a;
	double b = pars->b;
    allocate_external(Nx, Nt);
    RUNTIME_INFO(puts("Allocated memory for all external variables.\n"));
    mk_fourier_spectral_operators(Nx, a, b);
    RUNTIME_INFO(puts("Constructed fourier gridpoint.\n"));
    mk_initial_conditions(Nx, phi_init, dphi_init);
    RUNTIME_INFO(puts("Initialized the field and its temporal derivative.\n"));
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
    pars->cutoff_fraction = CUTOFF_FRACTION;
}

/*
allocate memory and initialize all external variables
*/
void allocate_external(size_t Nx, size_t Nt) {
    size_t N2 = Nx * Nx;

    //grid points
    x  = malloc(Nx * sizeof *x);
    if (!x)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    //solutions for the field and the temporal derivative (we are saving each
    //timestep: 2 * Nx * Nt space)
    field = calloc(N2 * Nt, sizeof *field);
    if (!field)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    // solution for the scale parameter a: Nt space
    frw_a = calloc(Nt, sizeof *frw_a);
    if (!frw_a)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    // T00 of the scalar field over time: Nt space
    rho = calloc(Nt, sizeof *rho);
    if (!rho)
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    // default array to save fourier coefficients from real to complex transform
    size_t nc = Nx / 2 + 1;
    cfftw_tmp = fftw_malloc(nc * sizeof *cfftw_tmp);
    if (!cfftw_tmp)
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    // general purpose 2*N double memory block for temporary use
    dmisc_tmp = calloc(N2, sizeof *dmisc_tmp);
    if (!dmisc_tmp)
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

    // Console output for debugging
#ifdef DEBUG
		puts("x");
        print_vector(x, N);
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
    frw_a[0] = 1.0;

    // Console output for debugging
#ifdef DEBUG
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
	return sin(x);
	// return tanh(pow(x, 8));
}

double dphi_init(double x) {
	// return -sin(x);
	return 0.0;
}

/*
free all allocated memory
*/
void free_all_external() {
	free(x);
	free(field);
    free(frw_a);
    free(rho);
    free(cfftw_tmp);
	RUNTIME_INFO(puts("Memory from all external variables freed.\n"));
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