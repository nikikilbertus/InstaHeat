#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <time.h>
#include "main.h"
#include "setup.h"
#include "RK4_stepper.h"
#include "filehandling.h"
#include "tests.h"

// simulation parameters
parameters_t pars;

// spatial grid points
double *grid;

// evolution of the field and spatial derivatives (2*Nx * Nt space required)
double *field;

// evolution of the scale parameter a for the FRW equations
double *frw_a;

// evolution of the total energy of the field
double *rho;

// default array to save fourier coefficients from real to complex transform
complex *cfftw_tmp_x;
complex *cfftw_tmp_y;
complex *cfftw_tmp_z;

// general purpose double memory blocks for temporary use
double *dtmp_x;
double *dtmp_y;
double *dtmp_z;

// times all dfts for timing analysis
double fftw_time = 0.0;

/*
--------------------------------main--------------------------------------------
*/
int main(int argc, const char * argv[]) {
    clock_t start = clock();

    int threadnum, threadinit;
    threadinit = fftw_init_threads();
    if (threadinit == 0)
    {
        fputs("Could not initialize fftw threads.", stderr);
        exit(EXIT_FAILURE);
    }
    threadnum = GRIDPOINTS_TOTAL > 5000 ? omp_get_max_threads() : 1;
    fftw_plan_with_nthreads(threadnum);
    RUNTIME_INFO(printf("Initialized fftw with %d thread(s)\n\n", threadnum));

#ifdef RUN_TESTS_ONLY
    pars.dt = 0.1;
    allocate_and_initialize_all(&pars);
    run_all_tests(&pars);
    free_all_external();
    return 0;
#endif

    int count = 0;
    for (double dt = 0.01; dt > 1e-3; dt /= 2., count++)
    {
    	pars.dt = dt;
    	allocate_and_initialize_all(&pars);

        size_t Nx = pars.x.N;
        size_t Ny = pars.y.N;
        size_t Nz = pars.z.N;
        size_t Ntot = Nx * Ny * Nz;
        size_t os = 2 * Ntot * (pars.Nt - 1);

    	run_RK4_stepper(&pars);

        char *prefix_field = "field";
		print_vector_to_file(field + os, Ntot, 1,
                                prefix_field, count);

        char *prefix_frw_a = "a";
        print_vector_to_file(frw_a, pars.Nt, 1, prefix_frw_a, count);

        char *prefix_rho = "energy";
        print_vector_to_file(rho, pars.Nt, 1, prefix_rho, count);

    	free_all_external();
        // break;
    }
    fftw_cleanup_threads();

    clock_t end = clock();
    double secs = (double)(end - start) / CLOCKS_PER_SEC;
    RUNTIME_INFO(printf("main took %f seconds.\n", secs));
    RUNTIME_INFO(printf("fftw took %f seconds (%.2f %%).\n\n",
                                    fftw_time, 100.*(fftw_time/secs)));
    return 0;
}