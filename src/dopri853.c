#include "main.h"
#if INTEGRATION_METHOD == DOPRI853

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <fftw3.h>
#include "dopri853.h"
#include "toolbox.h"
#include "io.h"

/**
 * @file dopri853.c
 *
 * @brief An 8th order explicit Dormand Prince integrator with 5th and 3rd order
 * error estimates for adaptive stepsizes.
 *
 * Only the function run_dopri853() is called from outside this file. It is the
 * only one that needs to be visible.
 *
 * @see <a href="http://numerical.recipes">Numerical Recipes</a>
 */

static void initialize_dopri853();
static void write_start_info();
static void wrap_up_dopri853();
static int perform_step(const double dt_try);
static void try_step(const double dt);
static double error(const double dt);
static int success(const double err, double *dt);
#ifdef ENABLE_STIFFNESSCHECK
static void check_for_stiffness(const double dt);
#endif
static void allocate_dopri853_values();
static void allocate_and_initialize_tolerances();
static void free_dopri853();

/**
 * @brief Holds intermediate evaluations of the right hand side and errors for
 * the Dormand Prince integrator.
 *
 * Holds pointers to memory blocks for the intermediate evaluations of
 * the right hand side of the partial differential equation (as determined in
 * mk_rhs(const double t, double *f, double *result) in toolbox.c as
 * well as memory blocks for the error estimates (5th and 3rd order).
 */
struct dopri853_values
{
    double *k2, *k3, *k4, *k5, *k6, *k7, *k8, *k9, *k10, *k_tmp;
    double *yerr, *yerr2;
};

/**
 * @brief Holds parameters for the Dormand Prince integrator.
 *
 * Most of these parameters should be self explanatory by their names. For more
 * information see <a href="http://numerical.recipes">Numerical Recipes</a>.
 */
struct dopri853_control
{
    double t; ///< The current time
    double t_old; ///< The previous time (on last time slice)
    double ti; ///< The initial time
    double tf; ///< The final time
    double dt; ///< The time step size
    double dt_did; ///< The previously used time step size
    double dt_next; ///< The proposed next time step size
    double dt_min; ///< The minimal permissible time step size
    size_t max_steps; ///< The maximal number of steps
    int n_stp; ///< The number of performed steps
    int n_ok; ///< The number of successful steps
    int n_bad; ///< The number of unsuccessful steps
    double beta; ///< Internal parameter for the error estimates
    double alpha; ///< Internal parameter for the error estimates
    double safe; ///< Internal parameter for the error estimates
    double minscale; ///< Minimal permissible rescaling of the time step size
    double maxscale; ///< Maximal permissible rescaling of the time step size
    double *a_tol; ///< Absolute tolerances
    double *r_tol; ///< Relative tolerances
    double err_old; ///< The previous error (on the last time slice)
    int reject; ///< Flag whether time step is rejected or accepted
    double eps; ///< Epsilon value for comparisons
    double dtlamb; ///< The multiple of dt used for stiffness detection
    int n_stiff; ///< Skips for stiffness detection
    int stiff; ///< Counter for stiffness detection (set)
    int nonstiff; ///< Counter for stiffness detection (reset)
};

/**
 * @brief The parameters for the Dormand Prince intergrator.
 * @see dopri853.h
 */
struct dopri853_control dp;

/**
 * @brief Holds intermediate fields as well as the error estimates.
 * @see dopri853.h
 */
struct dopri853_values dpv;

/**
 * @brief Run the Dormand Prince integrator to evolve the fields.
 *
 * This is the main routine of dopri853.c. It initializes the
 * integration parameters, allocates necessary memory, integrates the field
 * (assuming initial values are already present), writes the evolution to disk
 * as specified and frees the memory allocated earlier. All steps are reported
 * via the standard output.
 */
void run_dopri853()
{
    initialize_dopri853();
    write_start_info();
    prepare_and_save_timeslice();
    TIME(mon.integration = -get_wall_time());
    for (dp.n_stp = 0; ; ++dp.n_stp) {
        if (dp.t + dp.dt * 1.0001 > dp.tf) {
            dp.dt = dp.tf - dp.t;
            INFO(printf("\nOvershoot, new dt = %f\n", dp.dt));
        }
        if (perform_step(dp.dt)) {
            break;
        }
        if (dp.dt_did == dp.dt) {
            ++dp.n_ok;
        } else {
            ++dp.n_bad;
        }
        #ifdef DEBUG
        INFO(printf("did step: %d with dt: %f\n", dp.n_stp, dp.dt_did));
        #endif
        if ((dp.n_stp + 1) % pars.file.skip == 0) {
            save();
        }
        if (dp.t >= dp.tf) {
            puts("\nReached specified final time.\n");
            break;
        }
        if (dp.n_stp >= dp.max_steps) {
            puts("\nReached maximal number of steps.\n");
            break;
        }
        double time = get_wall_time() + mon.total;
        if (time > pars.max_runtime) {
            puts("\nReached specified maximal runtime.\n");
            break;
        }
        if (fabs(dp.dt_next) <= dp.dt_min) {
            puts("\n!!! Stepsize underflow.\n");
            break;
        }
        dp.dt = dp.dt_next;
    }
    TIME(mon.integration += get_wall_time());
    wrap_up_dopri853();
}

/**
 * @brief Initializes parameters in the struct dopri853_control dp of the
 * Dormand Prince integrator, and calls further allocation and initialization
 * functions.
 *
 * All fields of the struct dopri853_control dp are set according either as
 * fixed initial values or according to parameters entered by the user. All
 * parameters of the integration routine are saved in this struct.
 *
 * @note Changes in here are not recommended. All parameters that can/should be
 * chosen by the user are determined somewhere else and only copied here.
 */
static void initialize_dopri853()
{
    dp.t = pars.t.ti;
    dp.t_old = pars.t.ti;
    dp.ti = pars.t.ti;
    dp.tf = pars.t.tf;
    dp.dt = pars.t.dt;
    dp.dt_did = 0.0;
    dp.dt_next = pars.t.dt;
    dp.dt_min = MINIMAL_DELTA_T;
    dp.max_steps = MAX_STEPS;
    dp.n_stp = 0;
    dp.n_ok = 0;
    dp.n_bad = 0;
    dp.beta = BETA;
    dp.alpha = 1.0/8.0 - dp.beta * 0.2;
    dp.safe = SAFE;
    dp.minscale = SMALLEST_SCALING;
    dp.maxscale = LARGEST_SCALING;
    dp.err_old = 1.0e-4;
    dp.reject = 0;
    dp.eps = DBL_EPSILON;
    dp.dtlamb = 0.0;
    dp.n_stiff = 10;
    dp.stiff = 0;
    dp.nonstiff = 0;
}

/**
 * @brief Print information about the integration to stdout.
 */
static void write_start_info()
{
    INFO(puts("Initialized dopri853 parameters.\n"));
    allocate_dopri853_values();
    allocate_and_initialize_tolerances();
    INFO(puts("Starting dopri853 integration with:"));
    INFO(printf("  Initial time: %.17g\n", dp.ti));
    INFO(printf("  Final time: %.17g\n", dp.tf));
    INFO(printf("  Initial time step dt: %.17g\n", dp.dt));
    INFO(printf("  Minimal time step dt: %.17g\n", dp.dt_min));
    INFO(printf("  Max number of steps: %zu\n", dp.max_steps));
    // TODO: change this once I have more definitions
    INFO(printf("  Relative tolerance: %.17g\n", RELATIVE_TOLERANCE));
    INFO(printf("  Absolute tolerance: %.17g\n\n", ABSOLUTE_TOLERANCE));
}

/**
 * @brief Wraps up the integration by writing meta data to disk, print
 * information and free memory.
 */
static void wrap_up_dopri853()
{
    size_t index = pars.file.index;
    if (index != 0 && t_out.buf[index - 1] < dp.t) {
        prepare_and_save_timeslice();
    }
    #ifdef ENABLE_FOLLOWUP
    h5_write_followup();
    #endif
    free_dopri853();

    INFO(puts("Writing simulation meta data to disk.\n"));
    double val[1];
    val[0] = (double)dp.n_stp;
    h5_write_simple(H5_STEPS_TOTAL_NAME, val, 1, H5D_COMPACT);
    val[0] = (double)dp.n_ok;
    h5_write_simple(H5_STEPS_OK_NAME, val, 1, H5D_COMPACT);
    val[0] = (double)dp.n_bad;
    h5_write_simple(H5_STEPS_BAD_NAME, val, 1, H5D_COMPACT);

    INFO(puts("Finished dopri853:"));
    #ifdef ENABLE_TIMING
    INFO(printf("  Time: %f seconds\n", mon.integration));
    h5_write_simple(H5_RUNTIME_STEPPER_NAME, &mon.integration, 1, H5D_COMPACT);
    #endif
    INFO(printf("  Total steps: %d\n", dp.n_stp + 1));
    INFO(printf("  Good steps:  %d\n", dp.n_ok));
    INFO(printf("  Bad steps:   %d\n\n", dp.n_bad));
}

/**
 * @brief Evolves field forward in time by one step.
 *
 * @param[in] dt_try The initially proposed stepsize for this step.
 * @return Returns 0 in case of success and 1 in case of failure, in this case a
 * underflow of the specified minimal stepsize.
 *
 * Starting with the input @p dt_try this function repeatedly calls
 * try_step(const double dt), computes the normalized overall error by calling
 * error(const double dt), plugs that into success(const double err, double
 * *dt) to determine, whether the step is accepted or rejected. In the first
 * case, the fiels are updated and filtered and/or written to disk if required.
 * In case of failure, the step is retried with the updated stepsize.
 *
 * @note Various counters and global values (like the current time) are updated.
 * If the stepsize is decreased up to the minimal specified stepsize without
 * succeeding, a error message is reported and the exit code 1 is returned.
 */
static int perform_step(const double dt_try)
{
    double dt = dt_try;
    for ( ; ; ) {
        try_step(dt);
        double err = error(dt);
        if (success(err, &dt)) {
            break;
        }
        if (fabs(dt) <= fabs(dp.t) * dp.eps) {
            fputs("\n!!! Stepsize underflow.\n", stderr);
            return 1;
        }
    }
    #ifdef ENABLE_FFT_FILTER
    evo_flags.filter = 1;
    #endif
    if ((dp.n_stp + 1) % pars.file.skip == 0) {
        evo_flags.output = 1;
    }
    mk_rhs(dp.t + dt, field_new, dfield_new);
    evo_flags.output = 0;
    evo_flags.filter = 0;

    #ifdef ENABLE_STIFFNESSCHECK
    if ((dp.n_ok % dp.n_stiff) == 0 || dp.stiff > 0) {
        check_for_stiffness(dt_try);
    }
    #endif

    #pragma omp parallel for
    for (size_t i = 0; i < pars.Ntot; ++i) {
        field[i] = field_new[i];
        dfield[i] = dfield_new[i];
    }
    dp.t_old = dp.t;
    dp.t += (dp.dt_did = dt);
    pars.t.t = dp.t;
    t_out.tmp[0] = dp.t;
    return 0;
}

/**
 * @brief Tries one timestep with the specified stepsize.
 *
 * @param[in] dt The stepsize to use in the integration step.
 *
 * Computes all necessary (12) intermediate evaluations of the right hand side
 * of the differential equation (specified by mk_rhs(const double t, double *f,
 * double *result) in toolbox.c) to perform one step of the Dormand
 * Prince integrator, evolves the field forward in time by the given stepsize
 * and computes the error estimates w.r.t. the previous timeslice (still as
 * vectors, i.e. not a single value). Intermediate evaluations and errors are
 * stored in members of the struct dopri853_values dp.
 *
 * @note Contrary to perform_step(const double dt_try) this function is not
 * garuanteed to actually result in a step forward. It only proposes a candidate
 * for a step, which is evaluated by the errors it produced. It might get
 * discarded and recomputed with a different dt.
 */
static void try_step(const double dt)
{
    const size_t Ntot = pars.Ntot;
    const double t = dp.t;
    size_t i;
    // ------------ 1 ------------
    // is already done in perform_step

    // ------------ 2 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt * dpc.a21 * dfield[i];
    }
    mk_rhs(t + dpc.c2 * dt, dpv.k_tmp, dpv.k2);

    // ------------ 3 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a31 * dfield[i] + dpc.a32 * dpv.k2[i]);
    }
    mk_rhs(t + dpc.c3 * dt, dpv.k_tmp, dpv.k3);

    // ------------ 4 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a41 * dfield[i] + dpc.a43 * dpv.k3[i]);
    }
    mk_rhs(t + dpc.c4 * dt, dpv.k_tmp, dpv.k4);

    // ------------ 5 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a51 * dfield[i] + dpc.a53 * dpv.k3[i] + dpc.a54 * dpv.k4[i]);
    }
    mk_rhs(t + dpc.c5 * dt, dpv.k_tmp, dpv.k5);

    // ------------ 6 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a61 * dfield[i] + dpc.a64 * dpv.k4[i] + dpc.a65 * dpv.k5[i]);
    }
    mk_rhs(t + dpc.c6 * dt, dpv.k_tmp, dpv.k6);

    // ------------ 7 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a71 * dfield[i] + dpc.a74 * dpv.k4[i] + dpc.a75 * dpv.k5[i] +
             dpc.a76 * dpv.k6[i]);
    }
    mk_rhs(t + dpc.c7 * dt, dpv.k_tmp, dpv.k7);

    // ------------ 8 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a81 * dfield[i] + dpc.a84 * dpv.k4[i] + dpc.a85 * dpv.k5[i] +
             dpc.a86 * dpv.k6[i] + dpc.a87 * dpv.k7[i]);
    }
    mk_rhs(t + dpc.c8 * dt, dpv.k_tmp, dpv.k8);

    // ------------ 9 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a91 * dfield[i] + dpc.a94 * dpv.k4[i] + dpc.a95 * dpv.k5[i] +
             dpc.a96 * dpv.k6[i] + dpc.a97 * dpv.k7[i] + dpc.a98 * dpv.k8[i]);
    }
    mk_rhs(t + dpc.c9 * dt, dpv.k_tmp, dpv.k9);

    // ------------ 10 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a101 * dfield[i] + dpc.a104 * dpv.k4[i] + dpc.a105 * dpv.k5[i] +
             dpc.a106 * dpv.k6[i] + dpc.a107 * dpv.k7[i] + dpc.a108 * dpv.k8[i] +
             dpc.a109 * dpv.k9[i]);
    }
    mk_rhs(t + dpc.c10 * dt, dpv.k_tmp, dpv.k10);

    // ------------ 11 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a111 * dfield[i] + dpc.a114 * dpv.k4[i] + dpc.a115 * dpv.k5[i] +
             dpc.a116 * dpv.k6[i] + dpc.a117 * dpv.k7[i] + dpc.a118 * dpv.k8[i] +
             dpc.a119 * dpv.k9[i] + dpc.a1110 * dpv.k10[i]);
    }
    mk_rhs(t + dpc.c11 * dt, dpv.k_tmp, dpv.k2);

    // ------------ new dt ------------
    const double tpdt = t + dt;

    // ------------ 12 ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k_tmp[i] = field[i] + dt *
            (dpc.a121 * dfield[i] + dpc.a124 * dpv.k4[i] + dpc.a125 * dpv.k5[i] +
             dpc.a126 * dpv.k6[i] + dpc.a127 * dpv.k7[i] + dpc.a128 * dpv.k8[i] +
             dpc.a129 * dpv.k9[i] + dpc.a1210 * dpv.k10[i] + dpc.a1211 * dpv.k2[i]);
    }
    mk_rhs(tpdt, dpv.k_tmp, dpv.k3);

    // ------------ step ahead ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.k4[i] = dpc.b1 * dfield[i] + dpc.b6 * dpv.k6[i] + dpc.b7 * dpv.k7[i] +
            dpc.b8 * dpv.k8[i] + dpc.b9 * dpv.k9[i] + dpc.b10 * dpv.k10[i] +
            dpc.b11 * dpv.k2[i] + dpc.b12 * dpv.k3[i];

        field_new[i] = field[i] + dt * dpv.k4[i];
    }

    // ------------ error estimates ------------
    #pragma omp parallel for
    for (i = 0; i < Ntot; ++i) {
        dpv.yerr[i] = dpv.k4[i] - dpc.bhh1 * dfield[i] - dpc.bhh2 * dpv.k9[i] -
                      dpc.bhh3 * dpv.k3[i];
        dpv.yerr2[i] = dpc.er1 * dfield[i] + dpc.er6 * dpv.k6[i] +
                       dpc.er7 * dpv.k7[i] + dpc.er8 * dpv.k8[i] +
                       dpc.er9 * dpv.k9[i] + dpc.er10 * dpv.k10[i] +
                       dpc.er11 * dpv.k2[i] + dpc.er12 * dpv.k3[i];
    }
}

/**
 * @brief Computes the error of the previously tried step.
 *
 * @param[in] dt The previously tried stepsize.
 * @return The error of the previously tried stepsize.
 *
 * Computes the collective error of all the fields in the integration routine
 * normalized such that the threshold value is 1.
 */
static double error(const double dt)
{
    double err = 0.0, err2 = 0.0;
    #pragma omp parallel for reduction(+: err, err2)
    for (size_t i = 0; i < pars.Ntot; ++i) {
        double sk = dp.a_tol[i] +
            dp.r_tol[i] * MAX(fabs(field[i]), fabs(field_new[i]));
        double tmp = dpv.yerr[i] / sk;
        err2 += tmp * tmp;
        tmp = dpv.yerr2[i] / sk;
        err += tmp * tmp;
    }
    double deno = err + 0.01 * err2;
    if (deno <= 0.0) {
        deno = 1.0;
    }
    return dt * err * sqrt(1.0 / (pars.Ntot * deno));
}

/**
 * @brief Reports whether a step was successful and adjusts the stepsize
 * accordingly.
 *
 * @param[in] err The error estimate of the previously tried stepsize as
 * computed by error().
 * @param[inout] dt The stepsize that has previously been tried. Serves as input
 * to compute next stepsize in case of success and is changed for a new trial of
 * the current step in case of failure.
 * @return Returns 1 in case of success and 0 in case of failure.
 *
 * When the error @p err is smaller than 1, the current step is accepted. In
 * this case the new stepsize is computed and set accordingly. Otherwise the
 * current stepsize @p dt is adjusted and the current step is retried.
 * Various rules enter the scaling of the stepsize and various counters and
 * flags are set.
 */
static int success(const double err, double *dt)
{
    double scale;
    if (err <= 1.0) {
        if (err == 0.0) {
            scale = dp.maxscale;
        } else {
            scale = dp.safe * pow(err, - dp.alpha) * pow(dp.err_old, dp.beta);
            if (scale < dp.minscale) {
                scale = dp.minscale;
            }
            if (scale > dp.maxscale) {
                scale = dp.maxscale;
            }
        }
        if (dp.reject) {
            dp.dt_next = (*dt) * MIN(scale, 1.0);
        } else {
            dp.dt_next = (*dt) * scale;
        }
        #ifdef MAX_DT_HUBBLE_FRACTION
        const double minstep = MAX_DT_HUBBLE_FRACTION * sqrt(3.0 / rho_mean);
        if (dp.dt_next > minstep) {
            dp.dt_next = minstep;
        }
        #endif
        dp.err_old = MAX(err, 1.0e-4);
        dp.reject = 0;
        return 1;
    } else {
        scale = MAX(dp.safe * pow(err, - dp.alpha), dp.minscale);
        (*dt) *= scale;
        dp.reject = 1;
        return 0;
    }
}

#ifdef ENABLE_STIFFNESSCHECK
/**
 * @brief Performs a check whether the evolution becomes stiff and aborts if
 * necessary.
 *
 * If stiffness requirements are fulfilled on 15 subsequent time steps, the
 * integration is aborted.
 */
static void check_for_stiffness(const double dt)
{
    TIME(mon.stiffcheck -= get_wall_time());
    const size_t Ntot = pars.Ntot;
    double num = 0.0;
    double den = 0.0;
    #pragma omp parallel for private(tmp) reduction(+: num, den)
    for (size_t i = 0; i < Ntot; ++i) {
        double tmp = dfield_new[i] - dpv.k3[i];
        num += tmp * tmp;
        tmp = field_new[i] - dpv.k_tmp[i];
        den += tmp * tmp;
    }
    if (den > 0.0) {
        dp.dtlamb = dt * sqrt(num / den);
    }
    if (dp.dtlamb > 6.1) {
        dp.nonstiff = 0;
        dp.stiff += 1;
        if (dp.stiff == 15) {
            fprintf(stderr, "Potential stiffness detected at t = %f\n", dp.t);
            exit(EXIT_FAILURE);
        }
    } else {
        dp.nonstiff += 1;
        if (dp.nonstiff == 6) {
            dp.stiff = 0;
        }
    }
    TIME(mon.stiffcheck += get_wall_time());
}
#endif

/**
 * @brief Allocates memory for the Dormand Prince integrator.
 *
 * Allocates memory for the intermediate evaluations of the Dormand Prince
 * integration routine as well as temporary memory for the erros.
 */
static void allocate_dopri853_values()
{
    const size_t Ntot = pars.Ntot;
    dpv.k2 = fftw_malloc(Ntot * sizeof *dpv.k2);
    dpv.k3 = fftw_malloc(Ntot * sizeof *dpv.k3);
    dpv.k4 = fftw_malloc(Ntot * sizeof *dpv.k4);
    dpv.k5 = fftw_malloc(Ntot * sizeof *dpv.k5);
    dpv.k6 = fftw_malloc(Ntot * sizeof *dpv.k6);
    dpv.k7 = fftw_malloc(Ntot * sizeof *dpv.k7);
    dpv.k8 = fftw_malloc(Ntot * sizeof *dpv.k8);
    dpv.k9 = fftw_malloc(Ntot * sizeof *dpv.k9);
    dpv.k10 = fftw_malloc(Ntot * sizeof *dpv.k10);
    dpv.k_tmp = fftw_malloc(Ntot * sizeof *dpv.k_tmp);
    dpv.yerr = fftw_malloc(Ntot * sizeof *dpv.yerr);
    dpv.yerr2 = fftw_malloc(Ntot * sizeof *dpv.yerr2);
    if (!(dpv.k2 && dpv.k3 && dpv.k4 && dpv.k5 && dpv.k6 && dpv.k7 && dpv.k8 &&
          dpv.k9 && dpv.k10 && dpv.k_tmp && dpv.yerr && dpv.yerr2)) {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < Ntot; ++i) {
        dpv.k2[i] = 0.0;
        dpv.k3[i] = 0.0;
        dpv.k4[i] = 0.0;
        dpv.k5[i] = 0.0;
        dpv.k6[i] = 0.0;
        dpv.k7[i] = 0.0;
        dpv.k8[i] = 0.0;
        dpv.k9[i] = 0.0;
        dpv.k10[i] = 0.0;
        dpv.yerr[i] = 0.0;
        dpv.yerr2[i] = 0.0;
        dpv.k_tmp[i] = 0.0;
    }
    INFO(puts("Allocated memory for dopri853 variables.\n"));
}

/**
 * @brief Initializes the vectorized absolute and relative tolerances.
 *
 * From the definitions in the parameter file initializes the absolute and
 * relative tolerances used to compute the error. The error determines whether
 * a step size is accepted or rejected and how to choose the next adaptive step
 * size.
 */
static void allocate_and_initialize_tolerances()
{
    const size_t Ntot = pars.Ntot;
    /* const size_t N = pars.N; */
    dp.a_tol = fftw_malloc(Ntot * sizeof *dp.a_tol);
    dp.r_tol = fftw_malloc(Ntot * sizeof *dp.r_tol);
    // TODO: initialize values
    #pragma omp parallel for
    for (size_t i = 0; i < Ntot; ++i) {
        dp.a_tol[i] = ABSOLUTE_TOLERANCE;
        dp.r_tol[i] = RELATIVE_TOLERANCE;
    }
    /* double rtol = RELATIVE_TOLERANCE; */
    /* #pragma omp parallel for */
    /* for (size_t i = 0; i < 2 * N; ++i) { */
    /*     dp.r_tol[i] = rtol; */
    /* } */
    /* /1* rtol = 1.0e4 * RELATIVE_TOLERANCE; *1/ */
    /* #pragma omp parallel for */
    /* for (size_t i = 2 * N; i < 4 * N; ++i) { */
    /*     dp.r_tol[i] = rtol; */
    /* } */
    /* #pragma omp parallel for */
    /* for (size_t i = 4 * N; i < Ntot; ++i) { */
    /*     dp.r_tol[i] = RELATIVE_TOLERANCE; */
    /* } */
}

/**
 * @brief Frees memory for the Dormand Prince integrator.
 *
 * Frees all memory that was allocated for the intermediate evaluations of the
 * Dormand Prince integration routine as well as temporary memory for the
 * errors.
 */
static void free_dopri853()
{
    fftw_free(dp.a_tol);
    fftw_free(dp.r_tol);
    fftw_free(dpv.k2);
    fftw_free(dpv.k3);
    fftw_free(dpv.k4);
    fftw_free(dpv.k5);
    fftw_free(dpv.k6);
    fftw_free(dpv.k7);
    fftw_free(dpv.k8);
    fftw_free(dpv.k9);
    fftw_free(dpv.k10);
    fftw_free(dpv.k_tmp);
    fftw_free(dpv.yerr);
    fftw_free(dpv.yerr2);
    INFO(puts("Freed memory of dopri853 variables.\n"));
}

#endif
