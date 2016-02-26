#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <time.h>
#include <sys/time.h>
/* #include <gperftools/profiler.h> */
#include "main.h"
#include "setup.h"
#include "RK4_stepper.h"
#include "dopri853_stepper.h"
#include "filehandling.h"
#include "tests.h"

// -----------------------------global variables--------------------------------
// explanation in main.h
parameters_t pars;
double *grid;
double *time_buf;
double *field;
double *field_new;
double *dfield;
double *dfield_new;
double *phi_buf;
double *dphi_buf;
double *psi_buf;
double *dpsi_buf;
double phi_mean;
double phi_var;
double *phi_mean_buf;
double *phi_var_buf;
double dphi_mean;
double dphi_var;
double *dphi_mean_buf;
double *dphi_var_buf;
double psi_mean;
double psi_var;
double *psi_mean_buf;
double *psi_var_buf;
double dpsi_mean;
double dpsi_var;
double *dpsi_mean_buf;
double *dpsi_var_buf;
double *f_a_buf;
double *rho;
double *rho_buf;
double rho_mean;
double rho_var;
double *rho_mean_buf;
double *rho_var_buf;
double *pow_spec;
double *pow_spec_buf;
temporary_t tmp;
fftw_plan p_fw;
fftw_plan p_bw;
double fftw_time_exe  = 0.0;
double fftw_time_plan = 0.0;
double filter_time = 0.0;
double poisson_time = 0.0;
double h5_time_write = 0.0;

// -----------------------------main--------------------------------------------
int main(int argc, const char * argv[]) {

    #ifdef RUN_TESTS_ONLY
    allocate_and_initialize_all();
    run_all_tests();
    free_and_destroy_all();
    return 0;
    #endif

    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif

    allocate_and_initialize_all();
    h5_create_empty_by_path(DATAPATH);

    #ifdef ENABLE_PROFILER
    ProfilerStart("testprofile.prof");
    #endif

    // main call to integration routine
    run_rk4();
    /* run_dopri853(); */

    #ifdef ENABLE_PROFILER
    ProfilerStop();
    #endif

    h5_close(pars.file.id);
    free_and_destroy_all();

    #ifdef SHOW_TIMING_INFO
    double secs = get_wall_time() - start;
    RUNTIME_INFO(printf("main took %f seconds.\n", secs));
    RUNTIME_INFO(puts("as percentage of total, not mutually disjoint:"));
    RUNTIME_INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
                        fftw_time_exe, 100. * (fftw_time_exe / secs)));
    RUNTIME_INFO(printf("fftw planning took %f seconds (%.2f %%).\n",
                        fftw_time_plan, 100. * (fftw_time_plan / secs)));
    RUNTIME_INFO(printf("fft filtering took %f seconds (%.2f %%).\n",
                        filter_time, 100. * (filter_time / secs)));
    RUNTIME_INFO(printf("poisson equation took %f seconds (%.2f %%).\n",
                        poisson_time, 100. * (poisson_time / secs)));
    RUNTIME_INFO(printf("h5 write to disk took %f seconds (%.2f %%).\n",
                        h5_time_write, 100. * (h5_time_write / secs)));
    #endif
    return 0;
}

// -----------------------------timing------------------------------------------
#ifdef SHOW_TIMING_INFO
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        RUNTIME_INFO(puts("Could not get wall time, reurning 0.\n"));
        return 0.0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}
#endif
