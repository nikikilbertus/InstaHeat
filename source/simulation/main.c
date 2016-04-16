#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#ifdef ENABLE_PROFILER
#include <gperftools/profiler.h>
#endif
#include "setup.h"
#include "RK4_stepper.h"
#include "dopri853_stepper.h"
#include "rkf45.h"
#include "filehandling.h"
#include "tests.h"

/**
 * @file main.c
 * @brief Declaration of global variables, the main routine and one function for
 * wall clock time.
 */

struct parameters pars;
struct output t_out;
struct output phi;
struct output dphi;
struct output psi;
struct output dpsi;
struct output rho_out;
struct output phi_smry;
struct output dphi_smry;
struct output psi_smry;
struct output dpsi_smry;
struct output rho_smry;
struct output a_out;
struct output phi_ps;
struct output cstr;
double *field;
double *field_new;
double *dfield;
double *dfield_new;
double *rho;
double rho_mean;
double *pressure;
double pressure_mean;
double *filter;
struct k_grid kvec;
struct temporary tmp;
fftw_plan p_fw;
fftw_plan p_bw;
double fftw_time_exe  = 0.0;
double fftw_time_plan = 0.0;
double filter_time = 0.0;
double poisson_time = 0.0;
double h5_time_write = 0.0;

/**
 * @brief Main routine: Calls setup, integration and cleanup routines.
 */
int main(int argc, const char * argv[])
{
    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif

    allocate_and_initialize_all();

    #ifdef RUN_TESTS_ONLY
    run_all_tests();
    free_and_destroy_all();
    return 0;
    #endif

    #ifdef ENABLE_PROFILER
    ProfilerStart("testprofile.prof");
    #endif

    /* #if INTEGRATION_METHOD == RK4 */
    /* run_rk4(); */
    /* #elif INTEGRATION_METHOD == DOPRI853 */
    /* run_dopri853(); */
    /* #endif */
    run_rkf45();

    #ifdef ENABLE_PROFILER
    ProfilerStop();
    #endif

    #ifdef SHOW_TIMING_INFO
    double secs = get_wall_time() - start;
    INFO(printf("main took %f seconds.\n", secs));
    INFO(puts("as percentage of total, not mutually disjoint:"));
    INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
                        fftw_time_exe, 100. * (fftw_time_exe / secs)));
    INFO(printf("fftw planning took %f seconds (%.2f %%).\n",
                        fftw_time_plan, 100. * (fftw_time_plan / secs)));
    INFO(printf("fft filtering took %f seconds (%.2f %%).\n",
                        filter_time, 100. * (filter_time / secs)));
    INFO(printf("poisson equation took %f seconds (%.2f %%).\n",
                        poisson_time, 100. * (poisson_time / secs)));
    INFO(printf("h5 write to disk took %f seconds (%.2f %%).\n",
                        h5_time_write, 100. * (h5_time_write / secs)));

    INFO(puts("Writing runtimes to disk\n"));
    h5_write_parameter(H5_RUNTIME_TOTAL_NAME, &secs, 1);
    h5_write_parameter(H5_RUNTIME_FFTW_NAME, &fftw_time_exe, 1);
    h5_write_parameter(H5_RUNTIME_FFTWPLAN_NAME, &fftw_time_plan, 1);
    h5_write_parameter(H5_RUNTIME_FILTER_NAME, &filter_time, 1);
    h5_write_parameter(H5_RUNTIME_ELLIPTIC_NAME, &poisson_time, 1);
    h5_write_parameter(H5_RUNTIME_WRITEOUT_NAME, &h5_time_write, 1);
    #endif

    free_and_destroy_all();
    return 0;
}

#ifdef SHOW_TIMING_INFO
/**
 * @brief Get wall clock time for timing analysis.
 *
 * Only differences in wall clock times are used for timing analysis.
 *
 * @return The current wall clock time in seconds since a fixed time.
 */
double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        INFO(puts("Could not get wall time, reurning 0.\n"));
        return 0.0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}
#endif
