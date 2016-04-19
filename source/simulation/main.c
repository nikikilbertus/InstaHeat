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
#include "rkf45.h"
#include "rk4.h"
#include "dopri853.h"
#include "io.h"
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
struct monitor mon;

/**
 * @brief Main routine: Calls setup, integration and cleanup routines.
 */
int main(int argc, const char * argv[])
{
    double secs = 0.0;
    TIME(secs = -get_wall_time());

    allocate_and_initialize_all();

    #ifdef RUN_TESTS_ONLY
    run_all_tests();
    free_and_destroy_all();
    return 0;
    #endif

    #ifdef ENABLE_PROFILER
    ProfilerStart("testprofile.prof");
    #endif

    #if INTEGRATION_METHOD == RK4
    run_rk4();
    #elif INTEGRATION_METHOD == DOPRI853
    run_dopri853();
    #elif INTEGRATION_METHOD == RKF45
    run_rkf45();
    #endif

    #ifdef ENABLE_PROFILER
    ProfilerStop();
    #endif

    #ifdef SHOW_TIMING_INFO
    secs += get_wall_time();
    INFO(printf("main took %f seconds.\n", secs));
    INFO(puts("as percentage of total, not mutually disjoint:"));
    INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
            mon.fftw_time_exe, 100. * (mon.fftw_time_exe / secs)));
    INFO(printf("fftw planning took %f seconds (%.2f %%).\n",
            mon.fftw_time_plan, 100. * (mon.fftw_time_plan / secs)));
    INFO(printf("fft filtering took %f seconds (%.2f %%).\n",
            mon.filter_time, 100. * (mon.filter_time / secs)));
    INFO(printf("poisson equation took %f seconds (%.2f %%).\n",
            mon.poisson_time, 100. * (mon.poisson_time / secs)));
    INFO(printf("computing constraints took %f seconds (%.2f %%).\n",
            mon.cstr_time, 100. * (mon.cstr_time / secs)));
    INFO(printf("computing summaries took %f seconds (%.2f %%).\n",
            mon.smry_time, 100. * (mon.smry_time / secs)));
    INFO(printf("copying buffers took %f seconds (%.2f %%).\n",
            mon.copy_buffer_time, 100. * (mon.copy_buffer_time / secs)));
    INFO(printf("h5 write to disk took %f seconds (%.2f %%).\n",
            mon.h5_time_write, 100. * (mon.h5_time_write / secs)));

    INFO(puts("Writing runtimes to disk\n"));
    h5_write_parameter(H5_RUNTIME_TOTAL_NAME, &secs, 1);
    h5_write_parameter(H5_RUNTIME_FFTW_NAME, &mon.fftw_time_exe, 1);
    h5_write_parameter(H5_RUNTIME_FFTWPLAN_NAME, &mon.fftw_time_plan, 1);
    h5_write_parameter(H5_RUNTIME_FILTER_NAME, &mon.filter_time, 1);
    h5_write_parameter(H5_RUNTIME_ELLIPTIC_NAME, &mon.poisson_time, 1);
    h5_write_parameter(H5_RUNTIME_COPY_BUFFER_NAME, &mon.copy_buffer_time, 1);
    h5_write_parameter(H5_RUNTIME_CSTR_NAME, &mon.cstr_time, 1);
    h5_write_parameter(H5_RUNTIME_SMRY_NAME, &mon.smry_time, 1);
    h5_write_parameter(H5_RUNTIME_WRITEOUT_NAME, &mon.h5_time_write, 1);
    #endif

    INFO(printf("rhs was called %zu times.\n\n", mon.calls_rhs));
    double tmp = (double) mon.calls_rhs;
    h5_write_parameter(H5_COUNTER_RHS, &tmp, 1);

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
