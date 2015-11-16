#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// current field and temporal derivatives
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
        // size_t Ntot = pars.x.N * pars.y.N * pars.z.N;

        char filename[32];
        sprintf(filename, "field_%03d.txt", count);
        strcpy(pars.field_name, DATAPATH);
        strcat(pars.field_name, filename);
        file_create_empty(pars.field_name);

        run_RK4_stepper(&pars);

        char filename_a[32];
        char filepath_a[64];
        sprintf(filename_a, "a_%03d.txt", count);
        strcpy(filepath_a, DATAPATH);
        strcat(filepath_a, filename_a);
        file_create_empty(filepath_a);
        file_append_1d(frw_a, pars.Nt, 1, filepath_a);

        char filename_rho[32];
        char filepath_rho[64];
        sprintf(filename_rho, "rho_%03d.txt", count);
        strcpy(filepath_rho, DATAPATH);
        strcat(filepath_rho, filename_rho);
        file_create_empty(filepath_rho);
        file_append_1d(rho, pars.Nt, 1, filepath_rho);

    	free_all_external(&pars);
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