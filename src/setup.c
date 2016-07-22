#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "setup.h"
#include "toolbox.h"
#include "io.h"
#include "main.h"

/**
 * @file setup.c
 * @brief One time call to setup/initialization destroy/cleanup before and after
 * the simulation respectively.
 *
 * Only the functions allocate_and_init_all() and free_and_destroy_all()
 * are called from outside this file. Each of them is called exactly once per
 * simulation. Therefore, performance is not a priority in this file.
 */

static void init_rng();
static void init_threading();
static void init_parameters();
static void init_grid_pars();
static void init_time_pars();
static void init_file_pars();
static void init_monitoring();
static void allocate_external();
static void init_output(struct output *out, const size_t dim, const int mode);
static fftw_plan mk_fftw_plan(double *r, complex *c,
        const size_t Nx, const size_t Ny, const size_t Nz, const int dir);
static void check_simd_alignment();
static int get_simd_alignment_of(double *f);
static void mk_k_grid();
#ifdef ENABLE_FFT_FILTER
static void mk_filter_mask();
static double filter_window(const double x);
#endif
static void mk_initial_conditions();
#ifdef IC_FROM_DAT_FILE
static void init_from_dat();
#endif
#if INITIAL_CONDITIONS == IC_FROM_BUNCH_DAVIES
static void init_from_bunch_davies();
static void mk_bunch_davies(double *f, const double H, const double homo,
        const double gamma);
static void mk_kernel(double *ker, double *rr, const size_t N, gsl_function *f);
static double kernel_integrand(double k, void *params);
static double kernel_integrand_origin(double k, void *params);
static complex box_muller();
static void embed_grid(const complex *s, complex *d,
        const size_t nx, const size_t ny, const size_t nz,
        const size_t mx, const size_t my, const size_t mz);
#endif
#if INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
static void init_from_internal_function();
static void mk_x_grid(double *grid, const size_t Nx, const size_t Ny,
        const size_t Nz);
static double phi_init(const double x, const double y, const double z,
        const double *ph);
static double dphi_init(const double x, const double y, const double z,
        const double *ph);
static double wrapped_gaussian(const double x, const double y, const double z);
#endif
#ifdef SHOW_RUNTIME_INFO
static void print_flag_status();
#endif
static void mk_initial_psi();
static void destroy_and_cleanup_fftw();
static void free_external();

static gsl_rng *rng;

/**
 * @brief Successively calls the subroutines in this file necessary to setup
 * everything for the simulation.
 *
 * A single call to this function sets up everything for the simulation. After
 * this call, on of the available integration routines can be called. This
 * should be the first function called, see main.c.
 *
 * @see main.c
 */
void allocate_and_init_all()
{
    init_rng();
    init_threading();
    init_parameters();
    allocate_external();
    p_fw = mk_fftw_plan(field, tmp.phic, pars.x.N, pars.y.N, pars.z.N,
            FFTW_FORWARD);
    p_bw = mk_fftw_plan(field, tmp.phic, pars.x.N, pars.y.N, pars.z.N,
            FFTW_BACKWARD);
    check_simd_alignment();
    mk_k_grid();
    #ifdef ENABLE_FFT_FILTER
    mk_filter_mask();
    #endif
    mk_initial_conditions();
    h5_create_empty_by_path();
    gsl_rng_free(rng);
    INFO(print_flag_status());
}

/**
 * @brief Allocate memory for and set the seed of the GSL pseudo random number
 * generator.
 *
 * We use the _Mersenne Twister_, i.e. the MT19937 generator of Makoto
 * Matsumoto and Takuji Nishimura as implemented by GNU GSL.
 */
static void init_rng()
{
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (SEED < 1) {
        fputs("\n\nNeed positive seed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    gsl_rng_set(rng, SEED);
    INFO(puts("Initialized random number generator.\n"));
}

/**
 * @brief Initialize FFTW3 with the specified numbers of threads.
 *
 * If the parameter `THREAD_NUMBER` is `0`, FFTW3 is initialized with the return
 * value of `omp_set_num_threads()` threads. Note that this might not be the
 * optimal number. It is worthwhile to empirically find the optimal number of
 * threads.
 */
static void init_threading()
{
    int threadinit = fftw_init_threads();
    if (threadinit == 0) {
        fputs("\n\nCould not initialize FFTW3 threads.\n", stderr);
        exit(EXIT_FAILURE);
    }
    int threadnum = THREAD_NUMBER <= 0 ? omp_get_max_threads() : THREAD_NUMBER;
    omp_set_num_threads(threadnum);
    fftw_plan_with_nthreads(threadnum);
    INFO(printf("Running omp & FFTW3 with %d thread(s).\n\n", threadnum));
}

/**
 * @brief Initialize values in the struct parameters pars defined in
 * `main_template.h`.
 *
 * Most of the values come from preprocessor defines, which in turn are filled
 * from the external `parameters.sh` file before compilation. However, some
 * parameters are computed from others in non trivial ways. The struct pars
 * provides a flexible way to access all parameters in a global scope.
 */
static void init_parameters()
{
    init_grid_pars();
    init_time_pars();
    init_file_pars();
    init_monitoring();
    if (BUNCH_DAVIES_CUTOFF < 0) {
        fputs("\n\nNeed positive cutoff for bunch davies vacuum.\n", stderr);
        exit(EXIT_FAILURE);
    }
    pars.bunch_davies_cutoff = BUNCH_DAVIES_CUTOFF;
    if (MAX_RUNTIME < 0) {
        fputs("\n\nNeed positive maximal runtime.\n", stderr);
        exit(EXIT_FAILURE);
    }
    pars.max_runtime = MAX_RUNTIME;
}

/**
 * @brief Initialize grid related parameters.
 *
 * The simulation volume is specified by a number of parameters for each
 * dimension. One can adjust the box size in each dimension as well as the
 * number of gridpoints and the strides for the output. From those values some
 * more parameters are derived like the wave number (squared). For each
 * dimension we bunlde the values into a grid_dimension struct. Moreover we
 * check the effictive dimension of the problem and several output parameters.
 *
 * @see `parameters` and `grid_dimension` in `main_template.h`
 */
static void init_grid_pars()
{
    pars.x.N = GRIDPOINTS_X;
    pars.x.a = SPATIAL_LOWER_BOUND_X;
    pars.x.b = SPATIAL_UPPER_BOUND_X;
    pars.x.k = TWOPI / (pars.x.b - pars.x.a);
    pars.x.k2 = TWOPI * TWOPI / ((pars.x.b - pars.x.a) * (pars.x.b - pars.x.a));
    pars.x.stride = STRIDE_X;

    pars.y.N = GRIDPOINTS_Y;
    pars.y.a = SPATIAL_LOWER_BOUND_Y;
    pars.y.b = SPATIAL_UPPER_BOUND_Y;
    pars.y.k = TWOPI / (pars.y.b - pars.y.a);
    pars.y.k2 = TWOPI * TWOPI / ((pars.y.b - pars.y.a) * (pars.y.b - pars.y.a));
    pars.y.stride = STRIDE_Y;

    pars.z.N = GRIDPOINTS_Z;
    pars.z.a = SPATIAL_LOWER_BOUND_Z;
    pars.z.b = SPATIAL_UPPER_BOUND_Z;
    pars.z.k = TWOPI / (pars.z.b - pars.z.a);
    pars.z.k2 = TWOPI * TWOPI / ((pars.z.b - pars.z.a) * (pars.z.b - pars.z.a));
    pars.z.stride = STRIDE_Z;

    if (pars.x.N < 1 || pars.y.N < 1 || pars.z.N < 1) {
        fputs("\n\nNeed positive number of gridpoints.\n", stderr);
        exit(EXIT_FAILURE);
    }
    if (pars.x.N < pars.y.N || pars.y.N < pars.z.N) {
        fputs("\n\nOrder number of gridpoints Nx >= Ny >= Nz.\n", stderr);
        exit(EXIT_FAILURE);
    }
    if (pars.x.a >= pars.x.b || pars.y.a > pars.y.b || pars.z.a > pars.z.b) {
        fputs("\n\nInvalid spatial bounds, need a <= b.\n", stderr);
        exit(EXIT_FAILURE);
    }

    pars.N = pars.x.N * pars.y.N * pars.z.N;

    pars.x.outN = (pars.x.N + pars.x.stride - 1) / pars.x.stride;
    pars.y.outN = (pars.y.N + pars.y.stride - 1) / pars.y.stride;
    pars.z.outN = (pars.z.N + pars.z.stride - 1) / pars.z.stride;
    pars.outN = pars.x.outN * pars.y.outN * pars.z.outN;

    pars.dim = 3;
    if (pars.z.N == 1) {
        pars.dim = 2;
        pars.z.a = 0.0; pars.z.b = 0.0; pars.z.k = 0.0; pars.z.k2 = 0.0;
        if (pars.y.N == 1) {
            pars.dim = 1;
            pars.y.a = 0.0; pars.y.b = 0.0; pars.y.k = 0.0; pars.y.k2 = 0.0;
            if (pars.x.N == 1) {
                fputs("\n\nNeed at least two grid points.\n", stderr);
                exit(EXIT_FAILURE);
            }
        }
    }

    // due to the memory layout of fftw, we need different upper bounds in for
    // loops depending on the dimension, (the N gridpoints from the last
    // dimension are transformed to floor(N/2)+1 points in fourier space)
    switch (pars.dim) {
        case 1:
            pars.x.M = pars.x.N / 2 + 1;
            pars.y.M = 1;
            pars.z.M = 1;
            break;
        case 2:
            pars.x.M = pars.x.N;
            pars.y.M = pars.y.N / 2 + 1;
            pars.z.M = 1;
            break;
        case 3:
            pars.x.M = pars.x.N;
            pars.y.M = pars.y.N;
            pars.z.M = pars.z.N / 2 + 1;
            break;
    }
    pars.M = pars.x.M * pars.y.M * pars.z.M;
    pars.Next = 2 * pars.M;
    pars.Ntot = 4 * pars.N + 4 * pars.Next + 1;

    INFO(printf("Initialized grid in %zu dimension(s).\n", pars.dim));
    INFO(printf("  Gridpoints: X: %zu, Y: %zu, Z: %zu.\n",
                pars.x.N, pars.y.N, pars.z.N));
    INFO(printf("  N: %zu, Next: %zu, Ntot: %zu\n\n",
                pars.N, pars.Next, pars.Ntot));
}

/**
 * @brief Initialize time evolution related parameters.
 *
 * @see `timing` in `main_template.h`
 */
static void init_time_pars()
{
    pars.t.dt = DELTA_T;
    pars.t.t = INITIAL_TIME;
    pars.t.ti = INITIAL_TIME;
    pars.t.tf = FINAL_TIME;
    pars.t.Nt = ceil((pars.t.tf - pars.t.ti) / pars.t.dt) + 1;
    if (pars.t.dt <= 0.0 || pars.t.ti > pars.t.tf) {
        fputs("\n\nNeed positive delta t and initial < final time.\n", stderr);
        exit(EXIT_FAILURE);
    }
    #if INTEGRATION_METHOD == RK4
    if (pars.t.Nt > MAX_STEPS) {
        puts("For current settings MAX_STEPS would be exceeded!\n");
    }
    #endif
    #ifdef GSL_STEPPER
    if (GSL_OUTPUT_NUMBER > MAX_STEPS) {
        puts("For current settings MAX_STEPS would be exceeded!\n");
    }
    #endif
    INFO(puts("Initialized time parameters.\n"));
}

/**
 * @brief Initialize output related parameters.
 *
 * @see `file_parameters` in `main_template.h`
 */
static void init_file_pars()
{
    pars.file.index = 0;
    if (WRITE_OUT_BUFFER_NUMBER < 1) {
        fputs("\n\nNeed positive number for buffer size.\n", stderr);
        exit(EXIT_FAILURE);
    }
    pars.file.buf_size = WRITE_OUT_BUFFER_NUMBER;
    if (TIME_STEP_SKIPS < 1) {
        fputs("\n\nNeed positive number for time step skips.\n", stderr);
        exit(EXIT_FAILURE);
    }
    pars.file.skip = TIME_STEP_SKIPS;
}

/**
 * @brief Initialize monitoring and timing variables.
 *
 * @see `monitor` in `main_template.h`
 */
static void init_monitoring()
{
    mon.calls_rhs = 0;
    mon.fftw_exe = 0.0;
    mon.fftw_plan = 0.0;
    mon.filter = 0.0;
    mon.elliptic = 0.0;
    mon.integration = 0.0;
    mon.gw_sources = 0.0;
    mon.h5_write = 0.0;
    mon.cpy_buffers = 0.0;
    mon.cstr = 0.0;
    mon.smry = 0.0;
    mon.stiffcheck = 0.0;
    INFO(puts("Initialized monitoring variables.\n"));
}

/**
 * @brief Allocate memory for all external (i.e. global) variables.
 */
static void allocate_external()
{
    const size_t N = pars.N;
    const size_t Ntot = pars.Ntot;
    const size_t M = pars.M;
    #ifdef LARGE_OUTPUT
    const size_t outN = pars.outN;
    #endif

    // ---------------------------time, a---------------------------------------
    init_output(&t_out, 1, 1);
    init_output(&a_out, 1, 1);

    // ---------------------------full fields: phi, dphi, psi, dpsi, rho--------
    field = fftw_malloc(Ntot * sizeof *field);
    field_new = fftw_malloc(Ntot * sizeof *field_new);
    dfield = fftw_malloc(Ntot * sizeof *dfield);
    dfield_new = fftw_malloc(Ntot * sizeof *dfield_new);
    rho = fftw_malloc(N * sizeof *rho);
    pressure = fftw_malloc(N * sizeof *pressure);

    #ifdef OUTPUT_PHI
    init_output(&phi, outN, 0);
    #endif
    #ifdef OUTPUT_DPHI
    init_output(&dphi, outN, 0);
    #endif
    #ifdef OUTPUT_PSI
    init_output(&psi, outN, 0);
    #endif
    #ifdef OUTPUT_DPSI
    init_output(&dpsi, outN, 0);
    #endif
    #ifdef OUTPUT_RHO
    init_output(&rho_out, outN, 0);
    #endif

    // ---------------------------summaries-------------------------------------
    #ifdef OUTPUT_PHI_SMRY
    init_output(&phi_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    init_output(&dphi_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    init_output(&psi_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    init_output(&dpsi_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    init_output(&rho_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_PRESSURE_SMRY
    init_output(&p_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_H1_SMRY
    init_output(&h1_smry, SUMMARY_VALUES, 1);
    #endif
    #ifdef OUTPUT_H2_SMRY
    init_output(&h2_smry, SUMMARY_VALUES, 1);
    #endif

    // ---------------------------power spectra---------------------------------
    #ifdef OUTPUT_PHI_PS
    init_output(&phi_ps, POWER_SPECTRUM_BINS, 1);
    #endif
    #ifdef OUTPUT_PSI_PS
    init_output(&psi_ps, POWER_SPECTRUM_BINS, 1);
    #endif
    #ifdef OUTPUT_RHO_PS
    init_output(&rho_ps, POWER_SPECTRUM_BINS, 1);
    #endif
    #ifdef ENABLE_GW
    init_output(&gw, POWER_SPECTRUM_BINS, 1);
    #endif

    // ---------------------------constraints-----------------------------------
    #ifdef OUTPUT_CONSTRAINTS
    init_output(&cstr, NUMBER_CONSTRAINTS, 1);
    #endif

    // ---------------------------k grids---------------------------------------
    kvec.sq = fftw_malloc(M * sizeof *kvec.sq);
    kvec.x = fftw_malloc(M * sizeof *kvec.x);
    kvec.y = fftw_malloc(M * sizeof *kvec.y);
    kvec.z = fftw_malloc(M * sizeof *kvec.z);
    kvec.xf = fftw_malloc(M * sizeof *kvec.xf);
    kvec.yf = fftw_malloc(M * sizeof *kvec.yf);
    kvec.zf = fftw_malloc(M * sizeof *kvec.zf);
    #ifdef ENABLE_FFT_FILTER
    filter = fftw_malloc(M * sizeof *filter);
    #endif

    // ---------------------------default arrays for fourier transforms--------
    tmp.phic = fftw_malloc(M * sizeof *tmp.phic);
    tmp.xphic = fftw_malloc(M * sizeof *tmp.xphic);
    tmp.yphic = fftw_malloc(M * sizeof *tmp.yphic);
    tmp.zphic = fftw_malloc(M * sizeof *tmp.zphic);
    tmp.psic = fftw_malloc(M * sizeof *tmp.psic);
    tmp.dpsic = fftw_malloc(M * sizeof *tmp.dpsic);
    tmp.fc = fftw_malloc(M * sizeof *tmp.fc);
    tmp.deltarhoc = fftw_malloc(M * sizeof *tmp.deltarhoc);

    // ------------general purpose memory blocks for temporary use-------------
    tmp.xphi = fftw_malloc(N * sizeof *tmp.xphi);
    tmp.yphi = fftw_malloc(N * sizeof *tmp.yphi);
    tmp.zphi = fftw_malloc(N * sizeof *tmp.zphi);
    tmp.grad = fftw_malloc(N * sizeof *tmp.grad);
    tmp.lap = fftw_malloc(N * sizeof *tmp.lap);
    tmp.f = fftw_malloc(N * sizeof *tmp.f);
    tmp.deltarho = fftw_malloc(N * sizeof *tmp.deltarho);

    if (!(field && field_new && dfield && dfield_new && rho && pressure &&
        tmp.phic && tmp.xphic && tmp.yphic && tmp.zphic &&
        tmp.xphi && tmp.yphi && tmp.zphi && tmp.grad && tmp.lap && tmp.psic &&
        tmp.fc && tmp.deltarhoc && tmp.dpsic && tmp.f && tmp.deltarho &&
        kvec.sq && kvec.x && kvec.y && kvec.z && kvec.xf && kvec.yf &&
        kvec.zf)) {
        fputs("\n\nAllocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("Allocated memory for external variables.\n"));
}

/**
 * @brief Allocate and initialize a `output` struct.
 *
 * @param[out] out The `output` struct for which to set the dimension and
 * allocate memory.
 * @param[in] dim The dimension of the `output` structure
 * @param[in] mode If @p mode == 0, the field `out.tmp` of @p out is not
 * allocated. For @p mode != 0, memory for `out.tmp` will be allocated.
 *
 * @see `output` struct in `main_template.h`
 */
static void init_output(struct output *out, const size_t dim, const int mode)
{
    const size_t Nbuf = pars.file.buf_size;
    out->dim = dim;
    if (mode != 0) {
        out->tmp = calloc(out->dim, sizeof *out->tmp);
    }
    out->buf = calloc(Nbuf * out->dim, sizeof *out->buf);
    if (!out->buf || (mode != 0 && !out->tmp)) {
        fputs("\n\nAllocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Setup an FFTW3 plan for discrete Fourier transforms.
 *
 * @param[in] r A pointer to the real array used to setup the plan.
 * @param[in] c A pointer to the complex array used to setup the plan.
 * @param[in] Nx The number of gridpoints in the x direction.
 * @param[in] Ny The number of gridpoints in the y direction.
 * @param[in] Nz The number of gridpoints in the z direction.
 * @param[in] dir The direction of the Fourier transform, i.e. `FFTW_FORWARD` or
 * `FFTW_BACKWARD`.
 * @return The created FFTW3 plan.
 *
 * The plans have to be created __before__ the arrays involved in the planning
 * procedure are initialized. In some planning modes, the values of the arrays
 * are destroyed during planning. We create the plans for fixed global arrays
 * and reuse them for various different arrays. One has to be careful that later
 * arrays fulfil the memory alignment.
 *
 * @see `check_simd_alignment()`
 */
static fftw_plan mk_fftw_plan(double *r, complex *c,
        const size_t Nx, const size_t Ny, const size_t Nz, const int dir)
{
    TIME(mon.fftw_plan -= get_wall_time());
    fftw_plan plan;
    if (pars.dim == 1) {
        if (dir == FFTW_FORWARD) {
            plan = fftw_plan_dft_r2c_1d(Nx, r, c, FFTW_DEFAULT_FLAG);
        } else {
            plan = fftw_plan_dft_c2r_1d(Nx, c, r, FFTW_DEFAULT_FLAG);
        }
    } else if (pars.dim == 2) {
        if (dir == FFTW_FORWARD) {
            plan = fftw_plan_dft_r2c_2d(Nx, Ny, r, c, FFTW_DEFAULT_FLAG);
        } else {
            plan = fftw_plan_dft_c2r_2d(Nx, Ny, c, r, FFTW_DEFAULT_FLAG);
        }
    } else if (pars.dim == 3) {
        if (dir == FFTW_FORWARD) {
            plan = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, r, c, FFTW_DEFAULT_FLAG);
        } else {
            plan = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, c, r, FFTW_DEFAULT_FLAG);
        }
    } else {
        fputs("\n\nBUG: Dimension must be 1, 2 or 3.\n", stderr);
        exit(EXIT_FAILURE);
    }
    TIME(mon.fftw_plan += get_wall_time());
    INFO(puts("Constructed FFTW3 plan.\n"));
    return plan;
}

/**
 * @brief Checks that the different fields that are bundled in the larger field
 * array are correctly aligned for reusing FFTW3 plans on different fields.
 */
static void check_simd_alignment()
{
    int ref = get_simd_alignment_of(field);
    int a1 = get_simd_alignment_of(dfield);
    int a2 = get_simd_alignment_of(field_new);
    int a3 = get_simd_alignment_of(dfield_new);
    if (ref != a1 || ref != a2 || ref != a3) {
        fputs("\n\nAlignment error!\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("All field arrays are correctly aligned.\n"));
}

/**
 * @brief Gets the alignment of a single field.
 *
 * @return An integer encoding the memory alignment.
 */
static int get_simd_alignment_of(double *f)
{
    const int ref = fftw_alignment_of(f);
    int fail = 0;
    for (size_t i = 1; i < 4; ++i) {
        int test = fftw_alignment_of(f + i * pars.N);
        if (test != ref) {
            fail = 1;
            break;
        }
    }
    if (fail == 1) {
        fputs("\n\nAlignment error!\n", stderr);
        exit(EXIT_FAILURE);
    }
    return ref;
}

/**
 * @brief Construct grids for the k vectors and its square.
 *
 * Fills the members of the global `k_grid` struct `kvec` defined in
 * `main_template.h`.
 *
 * @note The Nx/2, Ny/2, Nz/2 entries of kvec.x, kvec.y, kvec.z are set to zero
 * for differentiation via discrete Fourier transforms. The structure members
 * kvec.xf, kvec.yf, kvec.zf contain those entries. They are also used normally
 * for k^2.
 *
 * @see `k_grid` in `main_template.h`
 */
static void mk_k_grid()
{
    const size_t Nx = pars.x.N, Ny = pars.y.N, Nz = pars.z.N;
    const size_t Mx = pars.x.M, My = pars.y.M, Mz = pars.z.M;
    #pragma omp parallel for
    for (size_t i = 0; i < Mx; ++i) {
        size_t osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j) {
            size_t osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k) {
                size_t id = osy + k;
                double k2 = pars.z.k2 * k * k;

                if (i > Nx / 2) {
                    kvec.x[id] = pars.x.k * ((int)i - (int)Nx);
                    kvec.xf[id] = kvec.x[id];
                    k2 += pars.x.k2 * (Nx - i) * (Nx - i);
                } else if (2 * i == Nx) {
                    kvec.x[id] = 0.0;
                    kvec.xf[id] = pars.x.k * i;
                    k2 += pars.x.k2 * i * i;
                } else {
                    kvec.x[id] = pars.x.k * i;
                    kvec.xf[id] = kvec.x[id];
                    k2 += pars.x.k2 * i * i;
                }

                if (j > Ny / 2) {
                    kvec.y[id] = pars.y.k * ((int)j - (int)Ny);
                    kvec.yf[id] = kvec.y[id];
                    k2 += pars.y.k2 * (Ny - j) * (Ny - j);
                } else if (2 * j == Ny) {
                    kvec.y[id] = 0.0;
                    kvec.yf[id] = pars.y.k * j;
                    k2 += pars.y.k2 * j * j;
                } else {
                    kvec.y[id] = pars.y.k * j;
                    kvec.yf[id] = kvec.y[id];
                    k2 += pars.y.k2 * j * j;
                }

                if (2 * k == Nz) {
                    kvec.z[id] = 0.0;
                    kvec.zf[id] = pars.z.k * k;
                } else {
                    kvec.z[id] = pars.z.k * k;
                    kvec.zf[id] = kvec.z[id];
                }
                kvec.sq[id] = k2;
            }
        }
    }
    kvec.k2_max = pars.x.k2 * (pars.x.N/2) * (pars.x.N/2) +
                  pars.y.k2 * (pars.y.N/2) * (pars.y.N/2) +
                  pars.z.k2 * (pars.z.N/2) * (pars.z.N/2);
    kvec.k_max = sqrt(kvec.k2_max);
    INFO(puts("Constructed grids for wave vectors.\n"));
}

/**
 * @brief Constructs the spatial grid.
 *
 * @param[out] grid A double array of size `Nx + Ny + Nz` that is filled up
 * with the grid values in each direction.
 * @param[in] Nx The number of gridpoints in the x direction.
 * @param[in] Ny The number of gridpoints in the y direction.
 * @param[in] Nz The number of gridpoints in the z direction.
 *
 * Since the grid is rectangular with uniform spacing in each direction, only
 * the discrete x, y and z values are computed.
 */
static void mk_x_grid(double *grid, const size_t Nx, const size_t Ny,
        const size_t Nz)
{
    const double ax = pars.x.a, ay = pars.y.a, az = pars.z.a;
    const double bx = pars.x.b, by = pars.y.b, bz = pars.z.b;
    #pragma omp parallel for
    for (size_t i = 0; i < Nx; ++i) {
        grid[i] = ax + (bx - ax) * i / Nx;
    }
    #pragma omp parallel for
    for (size_t j = Nx; j < Nx + Ny; ++j) {
        grid[j] = ay + (by - ay) * (j - Nx) / Ny;
    }
    #pragma omp parallel for
    for (size_t k = Nx + Ny; k < Nx + Ny + Nz; ++k) {
        grid[k] = az + (bz - az) * (k - Nx - Ny) / Nz;
    }
}

#ifdef ENABLE_FFT_FILTER
/**
 * @brief Construct an array for filtering out high wave number modes in
 * Fourier space.
 *
 * Fills the global array `filter` with multiplicative factors that can be
 * applied pointwise to a field in Fourier space to cut off high frequency
 * modes.
 *
 * @see `filter_window(const double x)`
 */
static void mk_filter_mask()
{
    const double kxmax = (pars.x.N / 2) * pars.x.k;
    const double kymax = (pars.y.N / 2) * pars.y.k;
    const double kzmax = (pars.z.N / 2) * pars.z.k;
    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        double frac = fabs(kvec.xf[i]) / kxmax;
        if (pars.dim > 1) {
            frac = MAX(frac, fabs(kvec.yf[i]) / kymax);
            if (pars.dim > 2) {
                frac = MAX(frac, fabs(kvec.zf[i]) / kzmax);
            }
        }
        filter[i] = filter_window(frac);
    }
    INFO(puts("Constructed filter mask.\n"));
}

/**
 * @brief The specific shape of the cutoff for high frequency modes.
 *
 * @param[in] x A ratio \f$k/k_{\max}\f$ for a given mode \f$k\f$.
 * @return The multiplier for the given mode @p x in the filtering process.
 *
 * We have found the exponential cutoff function described in the references to
 * work well for our purposes. It keeps more modes than the common two thirds
 * rule.
 *
 * @see `mk_filter_mask()`
 * @see [Computing Nearly Singular Solutions Using Pseudo-Spectral
 * Methods](http://arxiv.org/abs/math/0701337)
 * @see [Numerical Study of Nearly Singular Solutions of the 3-D Incompressible
 * Euler Equations](http://arxiv.org/abs/physics/0608126)
 * @see [On the stability of the unsmoothed Fourier method for hyperbolic
 * equations](http://link.springer.com/article/10.1007%2Fs002110050019)
 */
static double filter_window(const double x)
{
    // exponential cutoff smoothing
    return exp(-36.0 * pow(x, 36));

    // two thirds rule
    /* return x < 2.0 / 3.0 ? 1.0 : 0.0; */
}
#endif

/**
 * @brief Initialize the fields based on preprocessor defines given in the
 * parameter file.
 *
 * Depending on the chosen initial conditions from `parameter.sh` the
 * corresponding functions are called. Once the function returns, \f$\phi\f$,
 * \f$\dot{\phi}\f$, \f$\psi\f$, \f$\dot{\psi}\f$ and \f$a\f$ (bundled in
 * `field`) have been assigned their initial values.
 *
 * @see `h5_read_followup()`
 * @see `init_from_dat()`
 * @see `init_from_bunch_davies()`
 * @see `init_from_internal_function()`
 */
static void mk_initial_conditions()
{
    #pragma omp parallel for
    for (size_t i = 0; i < pars.Ntot; ++i) {
        field[i] = 0.0;
        dfield[i] = 0.0;
        field_new[i] = 0.0;
        dfield_new[i] = 0.0;
    }
    #if INITIAL_CONDITIONS == IC_FROM_H5_FILE
    h5_read_followup();
    #elif defined(IC_FROM_DAT_FILE)
    init_from_dat();
    #elif INITIAL_CONDITIONS == IC_FROM_BUNCH_DAVIES
    init_from_bunch_davies();
    #elif INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
    init_from_internal_function();
    #endif
    t_out.tmp[0] = pars.t.ti;
    a_out.tmp[0] = field[pars.Ntot - 1];
    INFO(puts("Initialized all fields on first time slice.\n"));
}

#ifdef IC_FROM_DAT_FILE
/**
 * @brief Read initial conditions from a .dat file.
 *
 * Currently for internal use only. TODO[fill]
 *
 * @see `read_initial_data()` in `io.c`
 * @see `mk_initial_psi()`
 */
static void init_from_dat()
{
    read_initial_data();
    INFO(puts("Read initial data from file.\n"));
    /* center(field + 2 * pars.N, pars.N); */
    /* center(field + 3 * pars.N, pars.N); */
    field[pars.Ntot - 1] = A_INITIAL;
    #if INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITHOUT_PSI
    mk_initial_psi();
    #endif
}
#endif

/**
 * @brief Given that the initial \f$\phi\f$, \f$\dot{\phi}\f$ and \f$a\f$ are
 * already provided in `field`, construct the corresponding \f$\psi\f$ and
 * \f$\dot{\psi}\f$.
 *
 * We use an elliptic equation from the Hamiltonian constraint combined with
 * the momentum contraint to compute \f$\psi\f$ and \f$\dot{\psi}\f$ from given
 * \f$\phi\f$, \f$\dot{\phi}\f$ and \f$a\f$. TODO[link]
 */
static void mk_initial_psi()
{
    TIME(mon.elliptic -= get_wall_time());
    const size_t N = pars.N;
    #pragma omp parallel for
    for (size_t i = 2 * N; i < 4 * N; ++i) {
        field[i] = 0.0;
    }
    evo_flags.output = 0;
    evo_flags.filter = 0;

    //TODO: new routine
    mk_rhs(pars.t.ti, field, dfield);
    const double a2 = field[pars.Ntot - 1] * field[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
    const double phi_mean = mean(field, N);
    const double dphi_mean = mean(field + N, N);
    const double ddphi_mean = mean(dfield + N, N);

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        tmp.deltarho[i] = field[N + i] - dphi_mean;
        tmp.f[i] = field[i] - phi_mean;
    }
    fft(tmp.deltarho, tmp.deltarhoc);
    fft(tmp.f, tmp.fc);

    #pragma omp parallel for
    for (size_t i = 1; i < pars.M; ++i) {
        tmp.phic[i] = (-dphi_mean * tmp.deltarhoc[i] + ddphi_mean * tmp.fc[i]) /
                ((-dphi_mean * dphi_mean + 2.0 * kvec.sq[i] / a2) * N);
    }
    tmp.phic[0] = 0.0;

    ifft(tmp.phic, field + 2 * N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        field[3 * N + i] = 0.5 * dphi_mean * tmp.f[i] - hubble * field[2 * N + i];
    }

    //TODO: old routine (gives same results up to tiny deviations)
    /* mk_gradient_squared_and_laplacian(field); */
    /* mk_rho_and_p(field); */
    /* const double a2 = field[pars.Ntot - 1] * field[pars.Ntot - 1]; */
    /* const double hubble = sqrt(rho_mean / 3.0); */
    /* const double phi_mean = mean(field, N); */
    /* const double dphi_mean = mean(field + N, N); */
    /* double extra1 = 0.0, extra2 = 0.0; */

    /* #pragma omp parallel for reduction(+: extra1, extra2) */
    /* for (size_t i = 0; i < N; ++i) { */
    /*     tmp.deltarho[i] = rho[i] - rho_mean; */
    /*     tmp.f[i] = dphi_mean * (field[i] - phi_mean); */
    /*     extra1 += field[N + i] * field[N + i]; */
    /*     extra2 += tmp.grad[i]; */
    /* } */
    /* const double extra = 0.5 * (extra1 - extra2 / a2) / N; */

    /* printf("a2=%g hubble=%g phimean=%g dphimean=%g extra1=%g, extra2=%g" */
    /*         "extra=%g\n", a2, hubble, phi_mean, dphi_mean, extra1/N, extra2/N, */
    /*         extra); */

    /* fft(tmp.deltarho, tmp.deltarhoc); */
    /* fft(tmp.f, tmp.fc); */
    /* #pragma omp parallel for */
    /* for (size_t i = 1; i < pars.M; ++i) { */
    /*     tmp.phic[i] = 0.5 * (tmp.deltarhoc[i] + */
    /*         3.0 * hubble * tmp.fc[i]) / ((- kvec.sq[i] / a2 + extra) * N); */
    /* } */
    /* tmp.phic[0] = 0.0; */

    /* ifft(tmp.phic, field + 2 * N); */
    /* #pragma omp parallel for */
    /* for (size_t i = 0; i < N; ++i) { */
    /*     field[3 * N + i] = 0.5 * tmp.f[i] - hubble * field[2 * N + i]; */
    /* } */
    TIME(mon.elliptic += get_wall_time());
    INFO(puts("Constructed psi and dot psi from existing phi and dot phi.\n"));
}

#if INITIAL_CONDITIONS == IC_FROM_BUNCH_DAVIES
/**
 * @brief Construct a Bunch Davies vacuum as initial conditions and then call
 * `mk_initial_psi()` for the corresponding \f$\psi\f$, \f$\dot{\psi}\f$.
 *
 * We build upon the values in the DEFROST paper (in the references) for the end
 * of inflation. DEFROST uses `INFLATON_MASS=5e-6` and `MASS=1`. We can freely
 * adjust the `INFLATON_MASS` to change the amplitude of the initial
 * fluctuations. By scaling `MASS`, we can shift the initial Hubble length.
 * (Initial values for H and \f$\dot{\phi}\f$ are scaled automatically according
 * to the value of `MASS`.) We recommend keeping the bos length fixed at 10 and
 * only rescaling `MASS`.
 *
 * @see [DEFROST: A New Code for Simulating Preheating after
 * Inflation](http://arxiv.org/abs/0809.4904)
 * @see `mk_bunch_davies(double *f, const double H, const double homo, const
 * double gamma)`
 * @see `mk_initial_psi()`
 */
static void init_from_bunch_davies()
{
    // directly from DEFROST(v1.0), factor in dphi0 and H0 adjusts modes
    const double phi0 = 1.0093430384226378929425913902459;
    const double dphi0 = - MASS * 0.7137133070120812430962278466136;
    const double hubble = MASS * 0.5046715192113189464712956951230;
    const double meff2 = MASS * MASS - 2.25 * hubble * hubble;
    if (meff2 <= 0.0) {
        fputs("\n\nThe effective mass squared is negative.\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("Initializing phi."));
    mk_bunch_davies(field, meff2, phi0, -0.25);
    INFO(puts("Initializing dot phi."));
    mk_bunch_davies(field + pars.N, meff2, dphi0, 0.25);
    field[pars.Ntot - 1] = A_INITIAL;
    mk_initial_psi();
}

/**
 * @brief Construct a Bunch Davies vacuum for \f$\phi\f$ and \f$\dot{\phi}\f$ as
 * initial conditions following the description in DEFROST.
 *
 * @param[out] f The field to initialize with a Bunch Davies spectrum.
 * @param[in] meff2 The effective mass.
 * @param[in] homo The homogeneous part of the initial field (background).
 * @param[in] gamma The exponent in equation TODO[link] of the initial spectrum.
 *
 * While `mk_kernel(double *ker, double *rr, const size_t N, gsl_function *f)`
 * computes the kernel on some support points, we use natural cubic splines to
 * interpolate the kernel to specific grid points. The interpolation is
 * implemented by GNU GSL.
 *
 * @see [DEFROST: A New Code for Simulating Preheating after
 * Inflation](http://arxiv.org/abs/0809.4904)
 * @see `mk_kernel(double *ker, double *rr, const size_t N, gsl_function *f)`
 * @see `box_muller()`
 * @see `embed_grid(const complex *s, complex *d, const size_t Nx, const size_t
 * Ny, const size_t Nz, const size_t Mx, const size_t My, const size_t Mz)`
 */
static void mk_bunch_davies(double *f, const double meff2, const double homo,
        const double gamma)
{
    size_t Nx, Ny, Nz, Mx, My, Mz;
    size_t cutoff = pars.bunch_davies_cutoff;
    double dk, dkx, dky, dkz, kcut;
    if (cutoff < 1) {
        Nx = pars.x.N, Ny = pars.y.N, Nz = pars.z.N;
    } else {
        Nx = 2 * cutoff;
        Ny = pars.y.N > 1 ? Nx : 1;
        Nz = pars.z.N > 1 ? Nx : 1;
        if (pars.x.N < Nx || pars.y.N < Ny || pars.z.N < Nz) {
            fputs("\n\nCutoff too large for the grid size.\n", stderr);
            exit(EXIT_FAILURE);
        }
    }
    Mx = (pars.y.N == 1 && pars.z.N == 1) ? Nx / 2 + 1 : Nx;
    My = pars.z.N == 1 ? Ny / 2 + 1 : Ny;
    Mz = Nz / 2 + 1;

    dkx = TWOPI / (pars.x.b - pars.x.a);
    dk = dkx;
    kcut = (Nx / 2 + 1) * dkx;
    if (pars.dim > 1) {
        dky = TWOPI / (pars.y.b - pars.y.a);
        dk *= dky;
        kcut = MIN(kcut, (Ny / 2 + 1) * dky);
        if (pars.dim > 2) {
            dkz = TWOPI / (pars.z.b - pars.z.a);
            dk *= dkz;
            kcut = MIN(kcut, (Nz / 2 + 1) * dkz);
        }
    }

    double params[] = {meff2, gamma, kcut};
    gsl_function func;
    func.params = params;
    const size_t Nker = 1e3;
    double *rr = malloc(Nker * sizeof *rr);
    double *ker = malloc(Nker * sizeof *ker);
    mk_kernel(ker, rr, Nker, &func);

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *fker = gsl_spline_alloc(gsl_interp_cspline, Nker);
    int s = gsl_spline_init(fker, rr, ker, Nker);
    if (s != GSL_SUCCESS) {
        fputs("\n\nInterpolation for Bunch Davies failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("  Initialized spline interpolation for Bunch Davies kernel."));

    const double fac = INFLATON_MASS / (Nx * Ny * Nz * sqrt(2.0 * dk));
    double *grid = malloc((Nx + Ny + Nz) * sizeof *grid);
    mk_x_grid(grid, Nx, Ny, Nz);

    double *ftmp = fftw_malloc(Nx * Ny * Nz * sizeof *ftmp);
    for (size_t i = 0; i < Nx; ++i) {
        size_t osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j) {
            size_t osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k) {
                const double r = sqrt(grid[i] * grid[i] +
                                      grid[Nx + j] * grid[Nx + j] +
                                      grid[Nx + Ny + k] * grid[Nx + Ny + k]);
                ftmp[osy + k] = fac * gsl_spline_eval(fker, r, acc);
            }
        }
    }
    free(grid);
    gsl_spline_free(fker);
    gsl_interp_accel_free(acc);
    free(rr);
    free(ker);
    complex *ftmpc = fftw_malloc(Mx * My * Mz * sizeof *ftmpc);
    fftw_plan p = mk_fftw_plan(ftmp, ftmpc, Nx, Ny, Nz, FFTW_FORWARD);
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(ftmp);

    for (size_t i = 0; i < Mx * My * Mz; ++i) {
        ftmpc[i] *= box_muller();
    }
    ftmpc[0] = homo;

    embed_grid(ftmpc, tmp.phic, Nx, Ny, Nz, pars.x.N, pars.y.N, pars.z.N);
    fftw_free(ftmpc);
    ifft(tmp.phic, f);
    INFO(puts("  Done with this field.\n"));
}

/**
 * @brief Construct the kernel for the Bunch Davies spectrum.
 *
 * @param[out] ker The kernel at @p N support points.
 * @param[out] rr The @p N support points for the kernel @p ker.
 * @param[in] N The Number of support points on which to compute the kernel.
 * @param[in] f The kernel integrand as a gsl_function structure.
 *
 * The kernel is given in TODO[link] and contains itself an integration over
 * wave vectors k. This integration is performed We use an adaptive integration
 * routine for oscillatory integrands from GNU GSL for the integration.
 *
 * @see `kernel_integrand(double k, void *params)`
 * @see `kernel_integrand_origin(double k, void *params)`
 */
static void mk_kernel(double *ker, double *rr, const size_t N, gsl_function *f)
{
    double a = 0.0;
    if (pars.dim == 1) {
        a = 1.0e-8;
    }
    const double abs = 1.0e-8;
    const double rel = 1.0e-6;
    const size_t limit = 3e3;
    double L = 0.0;
    size_t trig_levels = 3e3;
    const double rx = MAX(fabs(pars.x.a), fabs(pars.x.b));
    const double ry = MAX(fabs(pars.y.a), fabs(pars.y.b));
    const double rz = MAX(fabs(pars.z.a), fabs(pars.z.b));
    const double rmax = sqrt(rx * rx + ry * ry + rz * rz);
    const double dr = rmax / (N - 2);
    double r = dr;

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(limit);
    gsl_integration_workspace *c_ws = gsl_integration_workspace_alloc(limit);
    gsl_integration_qawo_table *wf =
        gsl_integration_qawo_table_alloc(1.0, L, GSL_INTEG_SINE, trig_levels);

    int s;
    double res, err;
    f->function = &kernel_integrand_origin;
    s = gsl_integration_qagiu(f, a, abs, rel, limit, ws, &res, &err);
    if (s != GSL_SUCCESS) {
        fputs("\n\nIntegration for Bunch Davies failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    rr[0] = 0.0;
    ker[0] = res;

    f->function = &kernel_integrand;
    for (size_t i = 1; i < N; ++i, r += dr) {
        rr[i] = r;
        gsl_integration_qawo_table_set(wf, r, L, GSL_INTEG_SINE);
        s = gsl_integration_qawf(f, a, abs, limit, ws, c_ws, wf, &res, &err);
        if (s != GSL_SUCCESS) {
            fputs("\n\nIntegration for Bunch Davies failed.\n", stderr);
            exit(EXIT_FAILURE);
        }
        ker[i] = res / r;
    }
    gsl_integration_qawo_table_free(wf);
    gsl_integration_workspace_free(c_ws);
    gsl_integration_workspace_free(ws);
    INFO(puts("  Integrated Bunch Davies kernel on support points."));
}

/**
 * @brief The integrand of the kernel for the Bunch Davies spectrum.
 *
 * @param[in] k The wave number k at which to compute the kernel.
 * @param[in] params The optional parameters in a gsl_function structure.
 * @return The kernel value at @p k.
 *
 * The parameter strucutre @p params contains the effective mass squared, the
 * exponent gamma and the cutoff wave number (in this order).
 *
 * @see TODO[link]
 * @see `mk_kernel()`
 */
static double kernel_integrand(double k, void *params)
{
    const double meff2 = *(double *) params;
    const double gamma = *(((double *) params) + 1);
    const double kcut = *(((double *) params) + 2);
    const double tmp = pow(meff2 + k * k, gamma) * exp(- pow(k / kcut, 36));
    if (pars.dim == 3) {
        return k * tmp / sqrt(PI);
    } else if (pars.dim == 2) {
        return tmp / sqrt(2.0);
    } else {
        return tmp / (k * sqrt(PI));
    }
}

/**
 * @brief The integrand of the kernel for the Bunch Davies spectrum at the
 * origin.
 *
 * @param[in] k The wave number k at which to compute the kernel.
 * @param[in] params The optional parameters in a gsl_function structure.
 * @return The kernel value at @p k.
 *
 * The parameter strucutre @p params contains the effective mass squared, the
 * exponent gamma and the cutoff wave number (in this order).
 *
 * @see TODO[link]
 * @see `kernel_integrand(double k, void *params)`
 */
static double kernel_integrand_origin(double k, void *params)
{
    return k * kernel_integrand(k, params);
}

/**
 * @brief Compute a sample from two independent standard normal distributed
 * random variables via the Box-Muller-Transform.
 *
 * @return A complex number whose real and imaginary part are the samples from
 * two independent standard normal random variables.
 *
 * @note The two independent uniformly distributed random values are computed
 * via the Mersenne Twister MT19937 generator of Makoto Matsumoto and Takuji
 * Nishimura.
 */
static complex box_muller()
{
    const double u1 = gsl_rng_uniform(rng), u2 = gsl_rng_uniform(rng);
    return sqrt(-2 * log(u1)) * cexp(TWOPI * u2 * I);
}

/**
 * @brief Embed a smaller grid in the Fourier domain into a larger grid in the
 * Fourier domain.
 *
 * @param[in] s The smaller complex source grid in the Fourier domain.
 * @param[in] d The larger complex destination grid in the Fourier domain.
 * @param[in] Nx Number of gridpoints in the x direction of the smaller grid.
 * @param[in] Ny Number of gridpoints in the y direction of the smaller grid.
 * @param[in] Nz Number of gridpoints in the z direction of the smaller grid.
 * @param[in] Mx Number of gridpoints in the x direction of the larger grid.
 * @param[in] My Number of gridpoints in the y direction of the larger grid.
 * @param[in] Mz Number of gridpoints in the z direction of the larger grid.
 */
static void embed_grid(const complex *s, complex *d,
        const size_t Nx, const size_t Ny, const size_t Nz,
        const size_t Mx, const size_t My, const size_t Mz)
{
    if (Nx > Mx || Ny > My || Nz > Mz) {
        fputs("\n\nCan't embed larger grid into smaller one.\n", stderr);
        exit(EXIT_FAILURE);
    }
    if (Nx % 2 != 0 || (Ny % 2 != 0 && Ny > 1) || (Nz % 2 != 0 && Nz > 1) ||
        Mx % 2 != 0 || (My % 2 != 0 && My > 1) || (Mz % 2 != 0 && Mz > 1)) {
        fputs("\n\nGrid embedding works only for even grids.\n", stderr);
        exit(EXIT_FAILURE);
    }
    if (((Ny == 1 || My == 1) && Ny != My) ||
        ((Nz == 1 || Mz == 1) && Nz != Mz)) {
        fputs("\n\nGrid embedding needs grids with equal dimension.\n", stderr);
        exit(EXIT_FAILURE);
    }
    const size_t nx = (Ny == 1 && Nz == 1) ? Nx / 2 + 1 : Nx;
    const size_t ny = Nz == 1 ? Ny / 2 + 1 : Ny;
    const size_t nz = Nz / 2 + 1;
    const size_t mx = (My == 1 && Mz == 1) ? Mx / 2 + 1 : Mx;
    const size_t my = Mz == 1 ? My / 2 + 1 : My;
    const size_t mz = Mz / 2 + 1;
    #pragma omp parallel for
    for (size_t i = 0; i < mx * my * mz; ++i) {
        d[i] = 0.0;
    }
    // TODO: check that this works for all dimensions
    #pragma omp parallel for
    for (size_t i = 0; i < nx; ++i) {
        size_t x1 = i * ny * nz;
        size_t x2a = (2 * i <= Nx ? i : mx - nx + i) * my * mz;
        size_t x2b = (2 * i == Nx && my != 1 ? mx - nx + i : 0) * my * mz;
        for (size_t j = 0; j < ny; ++j) {
            size_t y1 = j * nz;
            size_t y2a = (2 * j <= Ny ? j : my - ny + j) * mz;
            size_t y2b = (2 * j == Ny && mz != 1 ? my - ny + j : 0) * mz;
            for (size_t k = 0; k < nz; ++k) {
                complex val = s[x1 + y1 + k];
                if (x2b) {
                    d[x2a + y2a + k] = 0.5 * val;
                    d[x2b + y2a + k] = 0.5 * val;
                }
                if (y2b) {
                    d[x2a + y2a + k] = 0.5 * val;
                    d[x2a + y2b + k] = 0.5 * val;
                }
                if (x2b && y2b) {
                    d[x2a + y2a + k] = 0.25 * val;
                    d[x2a + y2b + k] = 0.25 * val;
                    d[x2b + y2a + k] = 0.25 * val;
                    d[x2b + y2b + k] = 0.25 * val;
                }
                if (!x2b && !y2b) {
                    d[x2a + y2a + k] = val;
                }
            }
        }
    }
}
#endif

#if INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
/**
 * @brief Construct initial conditions from internally defined functions.
 *
 * First the actual spatial grid values are computed by calling
 * `mk_x_grid(double *grid, const size_t Nx, const size_t Ny, const size_t Nz)`
 * and then \f$\phi\f$ and \f$\dot{\phi}\f$ are computed by `phi_init(const
 * double x, const double y, const double z, const double *ph)`,
 * `dphi_init(const double x, const double y, const double z, const double *ph)`
 * potentially using some random phases `ph`. The initial scale factor \f$a\f$
 * comes from the parameter file and \f$\psi\f$, \f$\dot{\psi}\f$ are then
 * computed from \f$\phi\f$ and \f$\dot{\phi}\f$.
 */
static void init_from_internal_function()
{
    const size_t Nx = pars.x.N, Ny = pars.y.N, Nz = pars.z.N;
    double *grid = malloc((Nx + Ny + Nz) * sizeof *grid);
    mk_x_grid(grid, Nx, Ny, Nz);
    const size_t Nmodes = 16;
    // random phases
    double *theta = calloc(Nmodes, sizeof *theta);
    for (size_t i = 0; i < Nmodes; ++i) {
        theta[i] = TWOPI * gsl_rng_uniform(rng);
    }
    for (size_t i = 0; i < Nx; ++i) {
        double x = grid[i];
        size_t osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j) {
            double y = grid[Nx + j];
            size_t osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k) {
                double z = grid[Nx + Ny + k];
                field[osy + k] = phi_init(x, y, z, theta);
                field[pars.N + osy + k] = dphi_init(x, y, z, theta);
            }
        }
    }
    free(grid);
    free(theta);
    field[pars.Ntot - 1] = A_INITIAL;
    mk_initial_psi();
}

/**
 * @brief The initial value of \f$phi\f$.
 *
 * @param[in] x The x coordinate where we want to evaluate \f$\phi\f$.
 * @param[in] y The y coordinate where we want to evaluate \f$\phi\f$.
 * @param[in] z The z coordinate where we want to evaluate \f$\phi\f$.
 * @param[in] ph Random phases for various modes.
 *
 * @return The value of phi at the specified coordinates @p x, @p y, @p z.
 *
 * @note This function was mostly used for getting started and debugging and
 * is subject to constant change. If one has an analytic/algorithmic expression
 * for physically interesting initial conditions, one can implement them here,
 * or call another function (like `wrapped_gaussian(const double x, const double
 * y, const double z)`
 */
static double phi_init(const double x, const double y, const double z,
        const double *ph)
{
    // localized for higgs metastability potential
    /* double phi0 = 0.04; */
    /* double lambda = 20.0; */
    /* if (pars.dim == 1) */
    /*     return phi0 * 0.5 * (1.0 + cos(x)) * exp(-lambda * x * x); */
    /* else if (pars.dim == 2) */
    /*     return phi0 * 0.25 * (1.0 + cos(x)) * (1.0 + cos(y)) * */
    /*         exp(-lambda * (x * x + y * y)); */
    /* else */
    /*     return phi0 * 0.125 * (1.0 + cos(x)) * (1.0 + cos(y)) * (1.0 + cos(z)) * */
    /*        exp(-lambda * (x * x + y * y + z * z)); */

    // some simple waves for notch or step potential simulations
    /* double frac = 0.4; // vary the ratio between \phi_0 and \delta \phi */
    /* double phi0 = 0.73 * frac; // only vary if you know exactly why */
    /* double deltaphi = phi0 / frac; */
    /* if (pars.dim == 1) */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1])); */
    /* else if (pars.dim == 2) */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1]) + */
    /*                  cos(1.0 * y + ph[2]) + cos(-1.0 * y + ph[3])); */
    /* else */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1]) + */
    /*                  cos(1.0 * y + ph[2]) + cos(-1.0 * y + ph[3]) + */
    /*                  cos(1.0 * z + ph[4]) + cos(-1.0 * z + ph[5])); */

    // very simple one for testing with phi squared potential
    /* double mean = 14.1421356; // for 50 e-fold inflation */
    /* double mean = 6.319569842; // somewhere at the end of 50 e-fold inflation */

    // compare_2, pos= 5500
    /* double mean = 0.0202977; */
    /* double amplitude = -2.26961e-06; */

    // compare_2, pos= 6000
    const double scale = 1.0e0;
    const double mean = 0.0510864;
    const double amplitude = -3.743790000000000e-07 * scale;

    if (pars.dim == 1) {
        return mean + amplitude * cos(x);
        /* return mean - amplitude * wrapped_gaussian(x, y, z); */
    } else if (pars.dim == 2) {
        /* return mean + amplitude * cos(x + y + ph[0]); */
        return mean - amplitude * wrapped_gaussian(x, y, z);
    } else {
        /* return mean + amplitude * */
        /*     (cos(x + y + z + ph[0]) + cos(-x + y + z + ph[1]) + */
        /*      cos(x - y + z + ph[2]) + cos(x + y - z + ph[3]) + */
        /*      cos(2.0 * x + y + z + ph[4]) + cos(x + 2.0 * y + z + ph[5]) + */
        /*      cos(x + y + 2.0 * z + ph[6]) + cos(-2.0 * x + y + z + ph[7]) + */
        /*      cos(x - 2.0 * y + z + ph[8]) + cos(x + y - 2.0 * z + ph[9]) + */
        /*      cos(3.0 * x + y + z + ph[10]) + cos(x + 3.0 * y + z + ph[11]) + */
        /*      cos(x + y + 3.0 * z + ph[12]) + cos(-3.0 * x + y + z + ph[13]) + */
        /*      cos(x - 3.0 * y + z + ph[14]) + cos(x + y - 3.0 * z + ph[15]) */
        /*      ); */
        return mean - amplitude * wrapped_gaussian(x, y, z);
    }
}

/**
 * @brief The initial value of \f$\dot{\phi}\f$.
 *
 * @param[in] x The x coordinate where we want to evaluate \f$\dot{\phi}\f$.
 * @param[in] y The y coordinate where we want to evaluate \f$\dot{\phi}\f$.
 * @param[in] z The z coordinate where we want to evaluate \f$\dot{\phi}\f$.
 * @param[in] ph Random phases for various modes.
 *
 * @return The value of dphi at the specified coordinates @p x, @p y, @p z.
 *
 * @note This function was mostly used for getting started and debugging and
 * is subject to constant change. If one has an analytic/algorithmic expression
 * for physically interesting initial conditions, one can implement them here,
 * or call another function (like `wrapped_gaussian(const double x, const double
 * y, const double z)`
 */
static double dphi_init(const double x, const double y, const double z,
        const double *ph)
{
    // compare_2, pos= 5500
    /* double mean = -0.00475989; */
    /* double amplitude = -2.91473e-09; */

    // compare_2, pos= 6000
    const double scale = 1.0e0;
    const double mean = 3.255190000000000e-04;
    const double amplitude = 1.742130000000000e-08 * scale;

    if (pars.dim == 1) {
        return (mean + amplitude * cos(x)) * MASS / MASS_KARSTEN;
    } else if (pars.dim == 2) {
        return (mean + amplitude * cos(x + y + ph[0])) *
            MASS / MASS_KARSTEN;
    } else {
        return (mean + amplitude *
            (cos(x + y + z + ph[0]) + cos(-x + y + z + ph[1]) +
             cos(x - y + z + ph[2]) + cos(x + y - z + ph[3]) +
             cos(2.0 * x + y + z + ph[4]) + cos(x + 2.0 * y + z + ph[5]))) *
             MASS / MASS_KARSTEN;
        /* return (mean + amplitude * cos(x + y + z + ph[0])) * MASS / 1.0e-2; */
    }

    /* return -0.089318193; // somewhere at end of 50 e-fold inflation */
}

/**
 * @brief A periodic version of a Gaussian for localized initial conditions.
 */
static double wrapped_gaussian(const double x, const double y, const double z)
{
    const double s = 0.5;
    double res = 0.0;
    size_t max;
    if (pars.dim == 1) {
        max = 32;
        for (size_t i = 1; i <= max; ++i) {
            res += exp(-0.5 * i * i * s * s) * (cos(i * x) + pow(-1.0, i + 1));
        }
    }
    if (pars.dim == 2) {
        max = 16;
        for (size_t i = 1; i <= max; ++i) {
            for (size_t j = 1; j <= max; ++j) {
                res += exp(-0.5 * (i * i + j * j) * s * s) *
                    (cos(i * x) + pow(-1.0, i + 1)) *
                    (cos(j * y) + pow(-1.0, j + 1));
            }
        }
    }
    if (pars.dim == 3) {
        max = 16;
        for (size_t i = 1; i <= max; ++i) {
            for (size_t j = 1; j <= max; ++j) {
                for (size_t k = 1; k <= max; ++k) {
                    res += exp(-0.5 * (i * i + j * j + k * k) * s * s) *
                        (cos(i * x) + pow(-1.0, i + 1)) *
                        (cos(j * y) + pow(-1.0, j + 1)) *
                        (cos(k * z) + pow(-1.0, k + 1));
                }
            }
        }
    }
    return res / TWOPI;
}
#endif

#ifdef SHOW_RUNTIME_INFO
/**
 * @brief Print status of simulation parameters defined via compiler macros to
 * stdout.
 */
static void print_flag_status()
{
    #ifdef ENABLE_FFT_FILTER
    INFO(puts("Frequency cutoff filtering enabled.\n"));
    #else
    INFO(puts("Filtering disabled.\n"));
    #endif
    #ifdef ENABLE_GW
    INFO(puts("Gravitational wave extraction enabled.\n"));
    #else
    INFO(puts("Gravitational waves disabled.\n"));
    #endif
    #ifdef IC_FROM_BUNCH_DAVIES
    INFO(printf("Cutting off Bunch Davies spectrum at N = %zu.\n\n",
                pars.bunch_davies_cutoff));
    #endif
    #ifdef ENABLE_FOLLOWUP
    INFO(puts("Output for followup simulation enabled.\n"));
    #else
    INFO(puts("Output for followup simulation disabled.\n"));
    #endif
}
#endif

/**
 * @brief Successively calls the subroutines in this file necessary to cleanup
 * everything after the simulation is done.
 *
 * A single call to this function cleans up everything after the simulation.
 * After this function returns, the program can exit. It should be the last
 * function called.
 *
 * @see `main.c`.
 */
void free_and_destroy_all()
{
    h5_close();
    destroy_and_cleanup_fftw();
    free_external();
}

/**
 * @brief Destroy FFTW3 plans and clean up FFTW3 threads.
 */
static void destroy_and_cleanup_fftw()
{
    fftw_destroy_plan(p_fw);
    fftw_destroy_plan(p_bw);
    fftw_cleanup_threads();
    INFO(puts("Destroyed FFTW3 plans.\n"));
}

/**
 * @brief Free memory of all external (i.e. global) variables.
 *
 * @note Everything allocated in `allocate_external()` must be freed here.
 */
static void free_external()
{
    fftw_free(field);
    fftw_free(field_new);
    fftw_free(dfield);
    fftw_free(dfield_new);
    fftw_free(rho);
    fftw_free(pressure);
    free(t_out.tmp);
    free(t_out.buf);
    free(a_out.tmp);
    free(a_out.buf);
    #ifdef OUTPUT_PHI
    free(phi.buf);
    #endif
    #ifdef OUTPUT_DPHI
    free(dphi.buf);
    #endif
    #ifdef OUTPUT_PSI
    free(psi.buf);
    #endif
    #ifdef OUTPUT_DPSI
    free(dpsi.buf);
    #endif
    #ifdef OUTPUT_RHO
    free(rho_out.buf);
    #endif

    #ifdef OUTPUT_PHI_SMRY
    free(phi_smry.tmp);
    free(phi_smry.buf);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    free(dphi_smry.tmp);
    free(dphi_smry.buf);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    free(psi_smry.tmp);
    free(psi_smry.buf);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    free(dpsi_smry.tmp);
    free(dpsi_smry.buf);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    free(rho_smry.tmp);
    free(rho_smry.buf);
    #endif
    #ifdef OUTPUT_PRESSURE_SMRY
    free(p_smry.tmp);
    free(p_smry.buf);
    #endif
    #ifdef OUTPUT_H1_SMRY
    free(h1_smry.tmp);
    free(h1_smry.buf);
    #endif
    #ifdef OUTPUT_H2_SMRY
    free(h2_smry.tmp);
    free(h2_smry.buf);
    #endif

    #ifdef OUTPUT_PHI_PS
    free(phi_ps.tmp);
    free(phi_ps.buf);
    #endif
    #ifdef OUTPUT_PSI_PS
    free(psi_ps.tmp);
    free(psi_ps.buf);
    #endif
    #ifdef OUTPUT_RHO_PS
    free(rho_ps.tmp);
    free(rho_ps.buf);
    #endif
    #ifdef ENABLE_GW
    free(gw.tmp);
    free(gw.buf);
    #endif

    #ifdef OUTPUT_CONSTRAINTS
    free(cstr.tmp);
    free(cstr.buf);
    #endif

    #ifdef ENABLE_FFT_FILTER
    fftw_free(filter);
    #endif
    fftw_free(kvec.sq);
    fftw_free(kvec.x);
    fftw_free(kvec.y);
    fftw_free(kvec.z);
    fftw_free(kvec.xf);
    fftw_free(kvec.yf);
    fftw_free(kvec.zf);
    fftw_free(tmp.phic);
    fftw_free(tmp.xphic);
    fftw_free(tmp.yphic);
    fftw_free(tmp.zphic);
    fftw_free(tmp.xphi);
    fftw_free(tmp.yphi);
    fftw_free(tmp.zphi);
    fftw_free(tmp.grad);
    fftw_free(tmp.lap);
    fftw_free(tmp.psic);
    fftw_free(tmp.dpsic);
    fftw_free(tmp.fc);
    fftw_free(tmp.deltarhoc);
    fftw_free(tmp.f);
    fftw_free(tmp.deltarho);
    INFO(puts("Freed memory of external variables.\n"));
}
