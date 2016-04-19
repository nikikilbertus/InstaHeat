#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "setup.h"
#include "toolbox.h"
#include "io.h"
#include "main.h"

/**
 * @file setup.c
 * @brief One time call to setup/initialization destroy/cleanup before and after
 * the simulation respectively.
 *
 * Only the functions allocate_and_initialize_all() and free_and_destroy_all()
 * are called from outside this file. Each of them is called exactly once per
 * simulation. Therefore, performance is not an issue in this file.
 */

/**
 * @brief Successively calls the subroutines in this file necessary to setup
 * everything for the simulation.
 *
 * A single call to this function sets up everything for the simulation. After
 * this call, on of the available integration routines can be started. This
 * should be the first function called, see main.c.
 */
void allocate_and_initialize_all()
{
    initialize_threading();
    initialize_parameters();
    allocate_external();
    mk_fftw_plans();
    mk_k_grid();
    #ifdef ENABLE_FFT_FILTER
    mk_filter_mask();
    #endif
    mk_initial_conditions();
    h5_create_empty_by_path();
    #ifdef ENABLE_FFT_FILTER
    INFO(puts("Frequency cutoff filtering enabled.\n"));
    #else
    INFO(puts("Filtering disabled.\n"));
    #endif
    #if PSI_METHOD == PSI_ELLIPTIC
    INFO(puts("Solving elliptic equation for psi at each timesetp.\n"));
    #elif PSI_METHOD == PSI_PARABOLIC
    INFO(puts("Integrating psi using the parabolic constraint.\n"));
    #elif PSI_METHOD == PSI_HYPERBOLIC
    INFO(puts("Integrating psi using the hyperbolic constraint.\n"));
    #endif
}

/**
 * @brief Initialize fftw3 with the specified numbers of threads.
 *
 * If the parameter THREAD_NUMBER is 0, we initialize fftw3 with
 * omp_set_num_threads() threads.
 */
void initialize_threading()
{
    int threadnum, threadinit;
    threadinit = fftw_init_threads();
    if (threadinit == 0) {
        fputs("\n\nCould not initialize fftw threads.\n", stderr);
        exit(EXIT_FAILURE);
    }
    threadnum = THREAD_NUMBER <= 0 ? omp_get_max_threads() : THREAD_NUMBER;
    omp_set_num_threads(threadnum);
    fftw_plan_with_nthreads(threadnum);
    INFO(printf("\n\nRunning omp & fftw with %d thread(s)\n\n", threadnum));
}

/**
 * @brief Initialize the values in the struct parameters pars.
 *
 * Most of the values come from preprocessor defines, which in turn are filled
 * from the external `parameters.sh` file before compilation. However, some
 * parameters are computed from others in non trivial ways. The struct pars
 * provides a flexible way to access all parameters in a global scope.
 */
void initialize_parameters()
{
    pars.x.N = GRIDPOINTS_X;
    pars.x.a = SPATIAL_LOWER_BOUND_X;
    pars.x.b = SPATIAL_UPPER_BOUND_X;
    pars.x.k = TWOPI / (pars.x.b - pars.x.a);
    pars.x.k2 = -TWOPI * TWOPI / ((pars.x.b - pars.x.a) * (pars.x.b - pars.x.a));
    pars.x.stride = STRIDE_X;

    pars.y.N = GRIDPOINTS_Y;
    pars.y.a = SPATIAL_LOWER_BOUND_Y;
    pars.y.b = SPATIAL_UPPER_BOUND_Y;
    pars.y.k = TWOPI / (pars.y.b - pars.y.a);
    pars.y.k2 = -TWOPI * TWOPI / ((pars.y.b - pars.y.a) * (pars.y.b - pars.y.a));
    pars.y.stride = STRIDE_Y;

    pars.z.N = GRIDPOINTS_Z;
    pars.z.a = SPATIAL_LOWER_BOUND_Z;
    pars.z.b = SPATIAL_UPPER_BOUND_Z;
    pars.z.k = TWOPI / (pars.z.b - pars.z.a);
    pars.z.k2 = -TWOPI * TWOPI / ((pars.z.b - pars.z.a) * (pars.z.b - pars.z.a));
    pars.z.stride = STRIDE_Z;

    pars.N = pars.x.N * pars.y.N * pars.z.N;
    pars.Nall = 4 * pars.N + 2;

    pars.x.outN = (pars.x.N + pars.x.stride - 1) / pars.x.stride;
    pars.y.outN = (pars.y.N + pars.y.stride - 1) / pars.y.stride;
    pars.z.outN = (pars.z.N + pars.z.stride - 1) / pars.z.stride;
    pars.outN = pars.x.outN * pars.y.outN * pars.z.outN;

    pars.dim = 3;
    if (pars.z.N == 1) {
        pars.dim = 2;
        if (pars.y.N == 1) {
            pars.dim = 1;
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

    // the evolution scheme for psi dictates which fields are evolved
    #if PSI_METHOD == PSI_ELLIPTIC
    pars.Ntot = 2 * pars.N + 1;
    #elif PSI_METHOD == PSI_PARABOLIC
    pars.Ntot = 3 * pars.N + 2;
    #elif PSI_METHOD == PSI_HYPERBOLIC
    pars.Ntot = 4 * pars.N + 2;
    #endif

    pars.t.dt = DELTA_T;
    pars.t.t = INITIAL_TIME;
    pars.t.ti = INITIAL_TIME;
    pars.t.tf = FINAL_TIME;
    pars.t.Nt = ceil((pars.t.tf - pars.t.ti) / pars.t.dt) + 1;
    if (pars.t.Nt > MAX_STEPS) {
        fputs("Exeeding MAX_STEPS, decrease DELTA_T.\n", stderr);
        exit(EXIT_FAILURE);
    }

    pars.file.index = 0;
    pars.file.buf_size = WRITE_OUT_BUFFER_NUMBER;
    pars.file.skip = TIME_STEP_SKIPS;

    mon.calls_rhs = 0;
    mon.fftw_time_exe = 0.0;
    mon.fftw_time_plan = 0.0;
    mon.filter_time = 0.0;
    mon.poisson_time = 0.0;
    mon.h5_time_write = 0.0;
    mon.copy_buffer_time = 0.0;
    mon.cstr_time = 0.0;
    mon.smry_time = 0.0;

    INFO(printf("Initialized parameters for %zu dimension(s).\n\n", pars.dim));
}

/**
 * @brief Allocate memory for all external (i.e. global) variables.
 */
void allocate_external()
{
    const size_t N = pars.N;
    const size_t Nall = pars.Nall;
    const size_t M = pars.M;
    const size_t bins = POWER_SPECTRUM_BINS;
    const size_t Nbuf = pars.file.buf_size;

    #ifdef LARGE_OUTPUT
    const size_t outN = pars.outN;
    #endif

    // ---------------------------time, a---------------------------------------
    t_out.dim = 1;
    t_out.tmp = calloc(t_out.dim, sizeof *t_out.tmp);
    t_out.buf = calloc(Nbuf * t_out.dim, sizeof *t_out.buf);
    a_out.dim = 1;
    a_out.tmp = calloc(a_out.dim, sizeof *a_out.tmp);
    a_out.buf = calloc(Nbuf * a_out.dim, sizeof *a_out.buf);

    // ---------------------------full fields: phi, dphi, psi, dpsi, rho--------
    field = fftw_malloc(Nall * sizeof *field);
    field_new = fftw_malloc(Nall * sizeof *field_new);
    dfield = fftw_malloc(Nall * sizeof *dfield);
    dfield_new = fftw_malloc(Nall * sizeof *dfield_new);
    rho = fftw_malloc(N * sizeof *rho);
    #if PSI_METHOD == PSI_HYPERBOLIC
    pressure = fftw_malloc(N * sizeof *pressure);
    #endif

    #ifdef OUTPUT_PHI
    phi.dim = outN;
    phi.buf = calloc(Nbuf * phi.dim, sizeof *phi.buf);
    #endif
    #ifdef OUTPUT_DPHI
    dphi.dim = outN;
    dphi.buf = calloc(Nbuf * dphi.dim, sizeof *dphi.buf);
    #endif
    #ifdef OUTPUT_PSI
    psi.dim = outN;
    psi.buf = calloc(Nbuf * psi.dim, sizeof *psi.buf);
    #endif
    #ifdef OUTPUT_DPSI
    dpsi.dim = outN;
    dpsi.buf = calloc(Nbuf * dpsi.dim, sizeof *dpsi.buf);
    #endif
    #ifdef OUTPUT_RHO
    rho_out.dim = outN;
    rho_out.buf = calloc(Nbuf * rho_out.dim, sizeof *rho_out.buf);
    #endif

    // ---------------------------summaries-------------------------------------
    #ifdef OUTPUT_PHI_SMRY
    phi_smry.dim = SUMMARY_VALUES;
    phi_smry.tmp = calloc(phi_smry.dim, sizeof *phi_smry.tmp);
    phi_smry.buf = calloc(Nbuf * phi_smry.dim, sizeof *phi_smry.buf);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    dphi_smry.dim = SUMMARY_VALUES;
    dphi_smry.tmp = calloc(dphi_smry.dim, sizeof *dphi_smry.tmp);
    dphi_smry.buf = calloc(Nbuf * dphi_smry.dim, sizeof *dphi_smry.buf);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    psi_smry.dim = SUMMARY_VALUES;
    psi_smry.tmp = calloc(psi_smry.dim, sizeof *psi_smry.tmp);
    psi_smry.buf = calloc(Nbuf * psi_smry.dim, sizeof *psi_smry.buf);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    dpsi_smry.dim = SUMMARY_VALUES;
    dpsi_smry.tmp = calloc(dpsi_smry.dim, sizeof *dpsi_smry.tmp);
    dpsi_smry.buf = calloc(Nbuf * dpsi_smry.dim, sizeof *dpsi_smry.buf);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    rho_smry.dim = SUMMARY_VALUES;
    rho_smry.tmp = calloc(rho_smry.dim, sizeof *rho_smry.tmp);
    rho_smry.buf = calloc(Nbuf * rho_smry.dim, sizeof *rho_smry.buf);
    #endif

    // ---------------------------power spectra---------------------------------
    #ifdef OUTPUT_PHI_PS
    phi_ps.dim = bins;
    phi_ps.tmp = calloc(phi_ps.dim, sizeof *phi_ps.tmp);
    phi_ps.buf = calloc(Nbuf * phi_ps.dim, sizeof *phi_ps.buf);
    #endif

    // ---------------------------constraints-----------------------------------
    #ifdef OUTPUT_CONSTRAINTS
    cstr.dim = NUMBER_CONSTRAINTS;
    cstr.tmp = calloc(cstr.dim, sizeof *cstr.tmp);
    cstr.buf = calloc(Nbuf * cstr.dim, sizeof *cstr.buf);
    #endif

    // ---------------------------k grids---------------------------------------
    kvec.sq = fftw_malloc(M * sizeof *kvec.sq);
    kvec.x = fftw_malloc(M * sizeof *kvec.x);
    kvec.y = fftw_malloc(M * sizeof *kvec.y);
    kvec.z = fftw_malloc(M * sizeof *kvec.z);
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

    if (!(field && field_new && dfield && dfield_new &&
        tmp.phic && tmp.xphic && tmp.yphic && tmp.zphic &&
        tmp.xphi && tmp.yphi && tmp.zphi && tmp.grad && tmp.lap && tmp.psic &&
        tmp.fc && tmp.deltarhoc && tmp.dpsic && tmp.f && tmp.deltarho)) {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("Allocated memory for external variables.\n"));
}

/**
 * @brief Setup the fftw plans for discrete fourier transforms.
 *
 * The plans have to be created __before__ the arrays involved in the planning
 * procedure are initialized. In some planning modes, the values of the arrays
 * are destroyed during planning. We create the plans for fixed global arrays
 * and reuse them for various different arrays. One has to be careful that later
 * arrays fulfil the memory alignment.
 *
 * @see FFTW3 documentation for more information on memory alignment (for SIMD
 * operations).
 */
void mk_fftw_plans()
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;

    TIME(mon.fftw_time_plan -= get_wall_time());
    switch (pars.dim) {
        case 1:
            p_fw = fftw_plan_dft_r2c_1d(Nx, field, tmp.phic,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_1d(Nx, tmp.phic, field,
                    FFTW_DEFAULT_FLAG);
            break;
        case 2:
            p_fw = fftw_plan_dft_r2c_2d(Nx, Ny, field, tmp.phic,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_2d(Nx, Ny, tmp.phic, field,
                    FFTW_DEFAULT_FLAG);
            break;
        case 3:
            p_fw = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, tmp.phic,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, tmp.phic, field,
                    FFTW_DEFAULT_FLAG);
            break;
    }
    TIME(mon.fftw_time_plan += get_wall_time());
    INFO(puts("Created fftw plans.\n"));
}

/**
 * @brief Construct grids for the k vector and its square.
 *
 * Fills the global arrays kvec.x, kvec.y, kvec.z and kvec.sq with the
 * corresponding components of the k vector and its square respectively.
 *
 * @note The Nx/2, Ny/2, Nz/2 entries of k_x, k_y, k_z are set to zero for
 * differentiation via discrete fourier transforms. However, those entries are
 * used normally for k^2.
 */
void mk_k_grid()
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t Mx = pars.x.M;
    const size_t My = pars.y.M;
    const size_t Mz = pars.z.M;

    double k2;
    size_t osx, osy, id;
    #pragma omp parallel for private(osx, osy, id, k2)
    for (size_t i = 0; i < Mx; ++i) {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j) {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k) {
                id = osy + k;
                k2 = pars.z.k2 * k * k;

                if (i > Nx / 2) {
                    kvec.x[id] = pars.x.k * ((int)i - (int)Nx);
                    k2 += pars.x.k2 * (Nx - i) * (Nx - i);
                } else if (2 * i == Nx) {
                    kvec.x[id] = 0.0;
                    k2 += pars.x.k2 * i * i;
                } else {
                    kvec.x[id] = pars.x.k * i;
                    k2 += pars.x.k2 * i * i;
                }

                if (j > Ny / 2) {
                    kvec.y[id] = pars.y.k * ((int)j - (int)Ny);
                    k2 += pars.y.k2 * (Ny - j) * (Ny - j);
                } else if (2 * j == Ny) {
                    kvec.y[id] = 0.0;
                    k2 += pars.y.k2 * j * j;
                } else {
                    kvec.y[id] = pars.y.k * j;
                    k2 += pars.y.k2 * j * j;
                }

                if (2 * k == Nz) {
                    kvec.z[id] = 0.0;
                } else {
                    kvec.z[id] = pars.z.k * k;
                }
                kvec.sq[id] = k2;
            }
        }
    }
    INFO(puts("Constructed grids for wave vectors.\n"));
}

/**
 * @brief Construct an arrray for filtering out high modes in Fourier space.
 *
 * Fills the global array `filter` with multiplicative factors that can be applied
 * pointwise to a field in Fourier space to cut off high frequency modes.
 */
void mk_filter_mask()
{
    const size_t N = pars.N;
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t Mx = pars.x.M;
    const size_t My = pars.y.M;
    const size_t Mz = pars.z.M;

    double tmp;
    size_t osx, osy;
    #pragma omp parallel for private(osx, osy, tmp)
    for (size_t i = 0; i < Mx; ++i) {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j) {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k) {
                tmp = filter_window(2.0 *
                    (i > Nx / 2 ? (int)Nx - (int)i : i) / (double) Nx);
                tmp *= filter_window(2.0 *
                    (j > Ny / 2 ? (int)Ny - (int)j : j) / (double) Ny);
                tmp *= filter_window(2.0 * k / (double) Nz);
                filter[osy + k] = tmp / (double) N;
            }
        }
    }
    INFO(puts("Constructed filter mask.\n"));
}

/**
 * @brief The specific shape of the cutoff for high frequency modes.
 *
 * We have found the exponential cutoff function described in the references to
 * work well for our purposes. It keeps more modes than the common two thirds
 * rule.
 *
 * @see `mk_filter_mask()`
 * @see [Computing Nearly Singular Solutions Using Pseudo-Spectral Methods](http://arxiv.org/abs/math/0701337)
 * @see [Numerical Study of Nearly Singular Solutions of the 3-D Incompressible Euler Equations](http://arxiv.org/abs/physics/0608126)
 * @see [On the stability of the unsmoothed Fourier method for hyperbolic equations](http://link.springer.com/article/10.1007%2Fs002110050019)
 */
inline double filter_window(const double x)
{
    return exp(-36.0 * pow(x, 36));

    // two thirds rule
    // return x < 2.0/3.0 ? x : 0.0;
}

/**
 * @brief Initializes the fields based on preprocessor defines given in the
 * parameter file.
 *
 * As a first step all fields are set to zero. Then the desired initialization
 * rountine is called. Once the function returns, $$\phi$$, $$\dot{\phi}$$,
 * $$\psi$$, $$\dot{\psi}$$, $$t$$ and $$a$$ have their initial values.
 */
void mk_initial_conditions()
{
    const size_t Nall = pars.Nall;
    #pragma omp parallel for
    for (size_t i = 0; i < Nall; ++i) {
        field[i] = 0.0;
        dfield[i] = 0.0;
        field_new[i] = 0.0;
        dfield_new[i] = 0.0;
    }

    #if INITIAL_CONDITIONS == IC_FROM_H5_FILE
    h5_read_timeslice();
    #elif INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITH_PSI || \
          INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITHOUT_PSI
    initialize_from_dat();
    #elif INITIAL_CONDITIONS == IC_FROM_BUNCH_DAVIES
    initialize_from_bunch_davies();
    #elif INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
    initialize_from_internal_function();
    #endif

    #ifdef EVOLVE_WITHOUT_PSI
    const size_t N2p = 2 * pars.N + 2;
    #pragma omp parallel for
    for (size_t i = N2p; i < Nall; ++i) {
        field[i] = 0.0;
    }
    #endif
    t_out.tmp[0] = pars.t.ti;
    a_out.tmp[0] = field[2 * pars.N];
    INFO(puts("Initialized fields on first time slice.\n"));
}

/**
 * @brief Read initial conditions from a .dat file.
 *
 * This is currently for internal use only. TODO[fill]
 *
 * @see `read_initial_data()` in `io.c`
 */
void initialize_from_dat()
{
    read_initial_data();
    /* center(field + 2 * pars.N + 2, pars.N); */
    /* center(field + 3 * pars.N + 2, pars.N); */
    field[2 * pars.N] = A_INITIAL;
    #if PSI_METHOD != PSI_ELLIPTIC && \
        INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITHOUT_PSI
    mk_initial_psi();
    #endif
}

/**
 * @brief Given that the initial $$\phi$$, $$\dot{\phi}$$ and $$a$$ are already
 * provided in `field`, construct the corresponding $$\psi$$ and
 * $$\dot{\psi}$$.
 */
void mk_initial_psi()
{
    const size_t N = pars.N;
    const size_t N2p = 2 * N + 2;
    const size_t Nall = pars.Nall;
    #pragma omp parallel for
    for (size_t i = N2p; i < Nall; ++i) {
        field[i] = 0.0;
    }
    mk_gradient_squared_and_laplacian(field);
    mk_rho(field);
    mk_psi(field);
}

/**
 * @brief Construct a Bunch Davies vacuum as initial conditions if the
 * parameters satisfy the conditions and then construct corresponding $$\psi$$,
 * $$\dot{psi}$$.
 *
 * We use the exact same values as in the DEFROST paper (in the references), hence we
 * need 3 dimensions, `MASS=1`, a box length of 10 in each direction. DEFROST uses
 * `INFLATON_MASS=5e-6`.
 *
 * @see [DEFROST: A New Code for Simulating Preheating after
 * Inflation](http://arxiv.org/abs/0809.4904)
 * @see `mk_bunch_davies(double *f, const double H, const double homo, const
 * double gamma)`
 */
void initialize_from_bunch_davies()
{
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    if (pars.dim != 3) {
        fputs("Bunch Davies vacuum works only in three dimensions.\n", stderr);
        exit(EXIT_FAILURE);
    }
    if (Nx != Ny || Nx != Nz || Ny != Nz) {
        fputs("Bunch Davies vacuum works only for Nx = Ny = Nz.\n", stderr);
        exit(EXIT_FAILURE);
    }
    double lx = fabs(pars.x.b - pars.x.a - 10.0);
    double ly = fabs(pars.y.b - pars.y.a - 10.0);
    double lz = fabs(pars.z.b - pars.z.a - 10.0);
    if (lx > DBL_EPSILON || ly > DBL_EPSILON || lz > DBL_EPSILON) {
        fputs("Bunch Davies vacuum works only for box size 10.0.\n", stderr);
        exit(EXIT_FAILURE);
    }
    const double phi0 = 1.0093430384226378929425913902459;
    const double dphi0 = -0.713706915863227;
    const double hubble = sqrt((dphi0 * dphi0 + MASS * MASS * phi0 * phi0) / 6.0);
    mk_bunch_davies(field, hubble, phi0, -0.25);
    mk_bunch_davies(field + pars.N, hubble, dphi0, 0.25);
    field[2 * pars.N] = A_INITIAL;
    #if PSI_METHOD != PSI_ELLIPTIC
    mk_initial_psi();
    #endif
}

/**
 * @brief Construct a Bunch Davies vacuum for $$\phi$$ and $$\dot{\phi}$$ as
 * initial conditions following the description in DEFROST.
 *
 * @see [DEFROST: A New Code for Simulating Preheating after Inflation](http://arxiv.org/abs/0809.4904)
 */
void mk_bunch_davies(double *f, const double H, const double homo,
        const double gamma)
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t N = pars.N;
    const size_t nn = Nx / 2 + 1;
    const size_t os = 16;
    const size_t nos = Nx * os * os;
    const double dx = (pars.x.b - pars.x.a) / Nx;
    const double dxos = dx / os;
    const double dk = TWOPI / (pars.x.b - pars.x.a);
    const double dkos = 0.5 * dk / os;
    //TODO: pspectre uses kcutpspectre = 2 * kcutdefrost
    const double kcut2 = 0.25 * nn * nn * dk * dk;
    const double meff2 = MASS * MASS - 2.25 * H * H;
    const double norm = 0.5 * INFLATON_MASS / (N * sqrt(TWOPI * pow(dk, 3))) *
        (dkos / dxos);

    if (meff2 <= 0.0) {
        fputs("The effective mass turned out to be negative.\n", stderr);
        exit(EXIT_FAILURE);
    }
    double *ker = fftw_malloc(nos * sizeof *ker);

    double kk;
    #pragma omp parallel for private(kk)
    for (size_t i = 0; i < nos; ++i) {
        kk = (i + 0.5) * dkos;
        ker[i] = kk * pow(kk * kk + meff2, gamma) *
            exp(-kk * kk / kcut2);
    }

    fftw_plan pl = fftw_plan_r2r_1d(nos, ker, ker, FFTW_RODFT10, FFTW_ESTIMATE);
    fftw_execute(pl);
    fftw_destroy_plan(pl);

    #pragma omp parallel for
    for (size_t i = 0; i < nos; ++i) {
        ker[i] *= norm / (i + 1);
    }

    size_t osx, osy, l;
    #pragma omp parallel for private(osx, osy, kk, l)
    for (int i = 0; i < Nx; ++i) {
        osx = i * Ny * Nz;
        for (int j = 0; j < Ny; ++j) {
            osy = osx + j * Nz;
            for (int k = 0; k < Nz; ++k) {
                kk = sqrt((double)( (i + 1 - nn) * (i + 1 - nn) +
                                    (j + 1 - nn) * (j + 1 - nn) +
                                    (k + 1 - nn) * (k + 1 - nn))) * os;
                l = (size_t) floor(kk);

                if (l > 0) {
                    f[osy + k] = ker[l - 1] + (kk - l) * (ker[l] - ker[l - 1]);
                } else {
                    f[osy + k] = (4.0 * ker[0] - ker[1]) / 3.0;
                }
            }
        }
    }

    fftw_free(ker);
    fftw_execute_dft_r2c(p_fw, f, tmp.phic);

    #pragma omp parallel for private(osx, osy)
    for (size_t i = 0; i < Nx; ++i) {
        osx = i * Ny * nn;
        for (size_t j = 0; j < Ny; ++j) {
            osy = osx + j * nn;
            for (size_t k = 0; k < nn; ++k) {
                tmp.phic[osy + k] *= box_muller();
            }
        }
    }

    tmp.phic[0] = homo;
    fftw_execute_dft_c2r(p_bw, tmp.phic, f);
}

/**
 * @brief Compute a sample from two independent standard normal distributed
 * random variables via the Box-Muller-Transform.
 *
 * @return A complex number whose real and imaginary part are the samples from
 * two independent standard normal random variables.
 *
 * @note The two independent uniformly distributed random values are computed
 * via the standard library rand() function.
 */
inline complex box_muller()
{
    const double u1 = (double)rand() / (double)RAND_MAX;
    const double u2 = (double)rand() / (double)RAND_MAX;
    return sqrt(-2 * log(u1)) * cexp(TWOPI * u2 * I);
}

/**
 * @brief Construct initial conditions for phi from internally defined functions
 *
 * First the actual spatial grid values are computed by calling
 * `mk_x_grid(double *grid)` and then $$\phi$$ and $$\dot{\phi}$$ are computed
 * by `phi_init(const double x, const double y, const double z, const double
 * *ph)`, `dphi_init(const double x, const double y, const double z, const
 * double *ph)` potentially using some random phases `ph`. The initial scale
 * factor $$a$$ comes from the parameter file and $$\psi$$, $$\dot{\psi}$$ are
 * then computed from $$\phi$$ and $$\dot{\phi}$$.
 */
void initialize_from_internal_function()
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t N = pars.N;
    size_t osx, osy;
    double x, y, z;
    double *grid = malloc((Nx + Ny + Nz) * sizeof *grid);
    mk_x_grid(grid);
    const size_t Nmodes = 16;
    double *theta = calloc(Nmodes, sizeof *theta);

    // random phases
    srand(SEED);
    for (size_t i = 0; i < Nmodes; ++i) {
        theta[i] = TWOPI * (double)rand() / (double)RAND_MAX;
    }

    for (size_t i = 0; i < Nx; ++i) {
        x = grid[i];
        osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j) {
            y = grid[Nx + j];
            osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k) {
                z = grid[Nx + Ny + k];
                field[osy + k] = phi_init(x, y, z, theta);
                field[N + osy + k] = dphi_init(x, y, z, theta);
            }
        }
    }
    free(grid);
    free(theta);
    field[2 * N] = A_INITIAL;
    #if PSI_METHOD != PSI_ELLIPTIC
    mk_initial_psi();
    #endif
}

/**
 * @brief Construct the spatial grid.
 *
 * @param[out] grid A double array of size `Nx + Ny + Nz` that is filled up
 * with the grid values in each direction.
 *
 * Since the grid is rectangular with uniform spacing in each direction, only
 * the x, y and z values are computed.
 */
void mk_x_grid(double *grid)
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const double ax = pars.x.a;
    const double bx = pars.x.b;
    const double ay = pars.y.a;
    const double by = pars.y.b;
    const double az = pars.z.a;
    const double bz = pars.z.b;

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
    INFO(puts("Constructed spatial grid.\n"));
}

/**
 * @brief The initial value of phi.
 *
 * @param[in] x The x coordinate where we want to evaluate $$\phi$$.
 * @param[in] y The y coordinate where we want to evaluate $$\phi$$.
 * @param[in] z The z coordinate where we want to evaluate $$\phi$$.
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
double phi_init(const double x, const double y, const double z,
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

    // karstens first set original
    /* double mean = 1.0; */
    /* double amplitude = 1.0e-5; */

    // karstens first set rescaling field instead of mass
    /* double mean = 6.0e3; */
    /* double amplitude = 6.0e-2; */

    // 25063
    /* double mean = 0.001543576559919; */
    /* double amplitude = 2.194936266463790e-6; */

    //10138
    /* double mean = 0.003801532800616; */
    /* double amplitude = 2.200533462675648e-7; */

    //1833
    /* double mean = 0.021019511647747; */
    /* double amplitude = 1.890808635066822e-6; */

    //30
    /* double mean = 0.999993224388493; */
    /* double amplitude = 1.000511520852772e-05; */

    //499
    /* double mean = 0.077381458703679; */
    /* double amplitude = -2.306228596956429e-07; */

    // compare_2, pos= 7671
    /* double mean = 0.0150052; */
    /* double amplitude = 4.048590000000000e-07; */

    // compare_2, pos= 5500
    /* double mean = 0.0202977; */
    /* double amplitude = -2.26961e-06; */

    // compare_2, pos= 6000
    const double scale = 1.0e0;
    const double mean = 0.0510864;
    const double amplitude = -3.743790000000000e-07 * scale;

    // compare_2, pos= 1
    /* double mean = 5.0; */
    /* double amplitude = 0.01; */

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
 * @brief The initial value of $$\dot{\phi}$$.
 *
 * @param[in] x The x coordinate where we want to evaluate $$\dot{\phi}$$.
 * @param[in] y The y coordinate where we want to evaluate $$\dot{\phi}$$.
 * @param[in] z The z coordinate where we want to evaluate $$\dot{\phi}$$.
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
double dphi_init(const double x, const double y, const double z,
        const double *ph)
{
    /* return 0.0; */

    /* double mean = -1.447595527218249e-8; */
    /* double amplitude = 1.794195724493731e-7; */

    /* double mean = 9.996023030535600e-9; */
    /* double amplitude = 1.794182821708652e-7; */

    /* double mean = -7.814852944111800e-8; */
    /* double amplitude = 1.791222773169365e-7; */

    /* double mean = -5.351102009251914e-06; */
    /* double amplitude = 1.865069892229237e-09; */

    /* double mean = 9.369552970351966e-07; */
    /* double amplitude = 1.768837536606555e-07; */

    /* double mean = -4.397960000000000e-06; */
    /* double amplitude = 1.816140000000000e-08; */

    /* double mean = -0.00475989; */
    /* double amplitude = -2.91473e-09; */

    const double scale = 1.0e0;
    const double mean = 3.255190000000000e-04;
    const double amplitude = 1.742130000000000e-08 * scale;

    /* double mean = -0.00806088; */
    /* double amplitude = -1.134420000000000e-20; */

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
double wrapped_gaussian(const double x, const double y, const double z)
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

/**
 * @brief Successively calls the subroutines in this file necessary to cleanup
 * everything after the simulation is done.
 *
 * A single call to this function cleans up everything after the simulation.
 * After this call, the program can exit. It should be the last function called,
 * see `main.c`.
 */
void free_and_destroy_all()
{
    h5_close();
    destroy_and_cleanup_fftw();
    free_external();
}

/**
 * @brief Destroy fftw plans and clean up threads.
 */
void destroy_and_cleanup_fftw()
{
    fftw_destroy_plan(p_fw);
    fftw_destroy_plan(p_bw);
    fftw_cleanup_threads();
    INFO(puts("Destroyed fftw plans.\n"));
}

/**
 * @brief Free memory from all external (i.e. global) variables.
 *
 * @note Everything allocated in `allocate_external()` must be freed here.
 */
void free_external()
{
    fftw_free(field);
    fftw_free(field_new);
    fftw_free(dfield);
    fftw_free(dfield_new);
    fftw_free(rho);
    #if PSI_METHOD == PSI_HYPERBOLIC
    fftw_free(pressure);
    #endif
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

    #ifdef OUTPUT_PHI_PS
    free(phi_ps.tmp);
    free(phi_ps.buf);
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
    INFO(puts("Freed external variables.\n"));
}

/**
 * @brief Print a vector to standard output.
 *
 * @param[in] vector The vector to print.
 * @param[in] N The length of the vector, i.e. how many etries to print.
 *
 * This function was mostly used for debugging, hence is obsolete/for internal
 * use only.
 */
void print_vector(const double *vector, const size_t N)
{
    for (size_t i = 0; i < N; i++) {
        printf("%f\n", vector[i]);
    }
}
