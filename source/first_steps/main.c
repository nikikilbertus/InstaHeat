#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <time.h>
#include <sys/time.h>
#include <gperftools/profiler.h>
#include "main.h"
#include "setup.h"
#include "RK4_stepper.h"
#include "dopri853_stepper.h"
#include "filehandling.h"
#include "tests.h"

// -----------------------------global variables--------------------------------
// simulation parameters
parameters_t pars;

// spatial grid points
double *grid;

// time buffer for write out
double *time_buf;

// scalar field phi and temporal derivative plus buffer for write out
double *field, *field_new;
double *dfield, *dfield_new;
double *field_buf;

// buffer for scaling parameter a
double *f_a_buf;

// total energy density rho of the field plus buffer for write out
double rho;
double *rho_buf;

// power spectrum plus buffer for write out
double *pow_spec;
double *pow_spec_buf;

// general purpose complex double memory blocks for temporary use
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

// time all dfts for timing analysis
double fftw_time_exe  = 0.0;
double fftw_time_plan = 0.0;
double h5_time_write = 0.0;

// -----------------------------main--------------------------------------------
int main(int argc, const char * argv[]) {

    #ifdef SHOW_TIMING_INFO
    double start, end;
    start = get_wall_time();
    #endif

    #ifdef RUN_TESTS_ONLY
    allocate_and_initialize_all();
    run_all_tests();
    free_and_destroy_all();
    return 0;
    #endif

    allocate_and_initialize_all();
    h5_create_empty_by_path(DATAPATH);

    #ifdef ENABLE_PROFILER
    ProfilerStart("testprofile.prof");
    #endif

    // main call to integration routine
    // run_rk4();
    run_dopri853();

    #ifdef ENABLE_PROFILER
    ProfilerStop();
    #endif

    h5_close(pars.file.id);
    free_and_destroy_all();
    fftw_cleanup_threads();

    #ifdef SHOW_TIMING_INFO
    end = get_wall_time();
    double secs = end - start;
    RUNTIME_INFO(printf("main took %f seconds.\n", secs));
    RUNTIME_INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
                                fftw_time_exe, 100.*(fftw_time_exe/secs)));
    RUNTIME_INFO(printf("fftw planning took %f seconds (%.2f %%).\n",
                                fftw_time_plan, 100.*(fftw_time_plan/secs)));
    RUNTIME_INFO(printf("h5 write to disk took %f seconds (%.2f %%).\n",
                                h5_time_write, 100.*(h5_time_write/secs)));
    #endif
    return 0;
}

// -----------------------------timing------------------------------------------
#ifdef SHOW_TIMING_INFO
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time, NULL)){
        RUNTIME_INFO(fputs("Could not get wall time.", stderr));
        return 0.0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
#endif