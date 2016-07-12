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
#include "rk4.h"
#include "dopri853.h"
#include "io.h"

/**
 * @file main.c
 * @brief Declaration of global variables, the main routine and a function for
 * timing purposes.
 */

#ifdef ENABLE_TIMING
static void write_monitoring();
#endif

#define LOGO ("********************************************************\n" \
      "      _____           _                         _       \n" \
      "      \\_   \\_ __  ___| |_ __ _  /\\  /\\___  __ _| |_     \n" \
      "       / /\\/ '_ \\/ __| __/ _` |/ /_/ / _ \\/ _` | __|    \n" \
      "    /\\/ /_ | | | \\__ \\ || (_| / __  /  __/ (_| | |_     \n" \
      "    \\____/ |_| |_|___/\\__\\__,_\\/ /_/ \\___|\\__,_|\\__|    \n" \
      "                                                        \n" \
      "                    Heating up...                       \n" \
      "********************************************************\n\n")

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
struct output p_smry;
struct output h1_smry;
struct output h2_smry;
struct output a_out;
struct output phi_ps;
struct output psi_ps;
struct output rho_ps;
struct output cstr;
struct output gw;
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
    mon.total = -get_wall_time();
    INFO(printf("%s", LOGO));
    allocate_and_init_all();

    #ifdef ENABLE_PROFILER
    ProfilerStart("profile.prof");
    #endif

    #if INTEGRATION_METHOD == RK4
    run_rk4();
    #elif INTEGRATION_METHOD == DOPRI853
    run_dopri853();
    #endif

    #ifdef ENABLE_PROFILER
    ProfilerStop();
    #endif

    mon.total += get_wall_time();
    INFO(printf("InstaHeat took %f seconds.\n\n", mon.total));
    h5_write_simple(H5_RUNTIME_TOTAL_NAME, &mon.total, 1, H5D_COMPACT);
    TIME(write_monitoring());
    INFO(printf("The rhs was evaluated %zu times.\n\n", mon.calls_rhs));
    double tmp = (double) mon.calls_rhs;
    INFO(puts("Final write to disk.\n"));
    h5_write_simple(H5_COUNTER_RHS, &tmp, 1, H5D_COMPACT);
    free_and_destroy_all();
    INFO(puts("Heat off...\n"));
    return 0;
}

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

#ifdef ENABLE_TIMING
/**
 * @brief Write monitoring info to h5 file and print on screen.
 */
static void write_monitoring()
{
    INFO(puts("as percentage of total, not mutually disjoint:"));
    INFO(printf("fftw execution took %f seconds (%.2f %%).\n",
            mon.fftw_exe, 100. * (mon.fftw_exe / mon.total)));
    INFO(printf("fftw planning took %f seconds (%.2f %%).\n",
            mon.fftw_plan, 100. * (mon.fftw_plan / mon.total)));
    INFO(printf("fft filtering took %f seconds (%.2f %%).\n",
            mon.filter, 100. * (mon.filter / mon.total)));
    INFO(printf("poisson equation took %f seconds (%.2f %%).\n",
            mon.elliptic, 100. * (mon.elliptic / mon.total)));
    INFO(printf("S_{ij}^{TT} took %f seconds (%.2f %%).\n",
            mon.gw_sources, 100. * (mon.gw_sources / mon.total)));
    INFO(printf("computing constraints took %f seconds (%.2f %%).\n",
            mon.cstr, 100. * (mon.cstr / mon.total)));
    INFO(printf("computing summaries took %f seconds (%.2f %%).\n",
            mon.smry, 100. * (mon.smry / mon.total)));
    INFO(printf("copying buffers took %f seconds (%.2f %%).\n",
            mon.cpy_buffers, 100. * (mon.cpy_buffers / mon.total)));
    INFO(printf("h5 write to disk took %f seconds (%.2f %%).\n",
            mon.h5_write, 100. * (mon.h5_write / mon.total)));
    INFO(printf("stiffness check took %f seconds (%.2f %%).\n",
            mon.stiffcheck, 100. * (mon.stiffcheck / mon.total)));

    INFO(puts("Writing runtimes to disk.\n"));
    h5_write_simple(H5_RUNTIME_FFTW_NAME, &mon.fftw_exe, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_FFTWPLAN_NAME, &mon.fftw_plan, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_FILTER_NAME, &mon.filter, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_ELLIPTIC_NAME, &mon.elliptic, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_STT_NAME, &mon.gw_sources, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_COPY_BUFFER_NAME, &mon.cpy_buffers, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_CSTR_NAME, &mon.cstr, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_SMRY_NAME, &mon.smry, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_STIFFCHECK_NAME, &mon.stiffcheck, 1, H5D_COMPACT);
    h5_write_simple(H5_RUNTIME_WRITEOUT_NAME, &mon.h5_write, 1, H5D_COMPACT);
}
#endif
