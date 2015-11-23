#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <time.h>
#include <sys/time.h>
#include <gperftools/profiler.h>
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

// power spectrum
double *pow_spec;

// default array to save fourier coefficients from real to complex transform
complex *cfftw_tmp;
complex *cfftw_tmp_x;
complex *cfftw_tmp_y;
complex *cfftw_tmp_z;

// general purpose double memory blocks for temporary use
double *dtmp_x;
double *dtmp_y;
double *dtmp_z;
double *dtmp_grad2;
double *dtmp_lap;

// frequently reused fftw plans
fftw_plan p_fw_3d;
fftw_plan p_bw_3d;

// times all dfts for timing analysis
double fftw_time_exe  = 0.0;
double fftw_time_plan = 0.0;

/*
--------------------------------main--------------------------------------------
*/
int main(int argc, const char * argv[]) {

    double start, end;
    start = get_wall_time();

    int threadnum, threadinit;
    threadinit = fftw_init_threads();
    if (threadinit == 0)
    {
        fputs("Could not initialize fftw threads.", stderr);
        exit(EXIT_FAILURE);
    }
    threadnum = 1;
    omp_set_num_threads(threadnum);
    // threadnum = GRIDPOINTS_TOTAL > 5000 ? omp_get_max_threads() : 1;
    fftw_plan_with_nthreads(threadnum);
    RUNTIME_INFO(printf("Initialized openmp with %d thread(s)\n\n", threadnum));

#ifdef RUN_TESTS_ONLY
    pars.t.dt = 0.1;
    allocate_and_initialize_all(&pars);
    run_all_tests(&pars);
    free_and_destroy_all(&pars);
    return 0;
#endif

    int count = 0;
    for (double dt = 0.1; dt > 1e-3; dt /= 2.0, count += 1)
    {
    	pars.t.dt = (double)GRIDPOINTS_X / 1000.0;//dt;
    	allocate_and_initialize_all(&pars);

        file_create_empty_by_name(pars.file.name_field);
        file_create_empty_by_name(pars.file.name_powspec);

        #ifdef ENABLE_PROFILER
        ProfilerStart("testprofile.prof");
        #endif

        run_RK4_stepper(&pars);

        #ifdef ENABLE_PROFILER
        ProfilerStop();
        #endif

        file_single_write_by_name_1d(frw_a, pars.t.Nt, 1, "a_000");
        file_single_write_by_name_1d(rho, pars.t.Nt, 1, "rho_000");

    	free_and_destroy_all(&pars);
        break;
    }
    fftw_cleanup_threads();

    end = get_wall_time();
    double secs = end - start;
    RUNTIME_INFO(printf("main took %f seconds.\n", secs));
    RUNTIME_INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
                                fftw_time_exe, 100.*(fftw_time_exe/secs)));
    RUNTIME_INFO(printf("fftw planning took %f seconds (%.2f %%).\n\n",
                                fftw_time_plan, 100.*(fftw_time_plan/secs)));
    return 0;
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        RUNTIME_INFO(fputs("Could not get wall time.", stderr));
        return 0.0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}