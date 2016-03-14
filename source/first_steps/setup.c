#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "setup.h"
#include "evolution_toolkit.h"
#include "filehandling.h"
#include "main.h"

void allocate_and_initialize_all()
{
    initialize_threading();
    initialize_parameters();
    allocate_external();
    mk_grid();
    mk_fftw_plans();
    mk_initial_conditions();
    mk_filter_mask();
    h5_create_empty_by_path(DATAPATH);
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
 *  initialize the values in the paramters_t pars variable, mostly from defines
 *  in main.h; using the struct gives more flexibility than using the defines
 *  throughout the code
 */
void initialize_parameters()
{
    pars.x.N  = GRIDPOINTS_X;
    pars.x.a  = SPATIAL_LOWER_BOUND_X;
    pars.x.b  = SPATIAL_UPPER_BOUND_X;
    pars.x.k  = TWOPI * I / (pars.x.b - pars.x.a);
    pars.x.k2 = -TWOPI * TWOPI / ((pars.x.b - pars.x.a) * (pars.x.b - pars.x.a));
    pars.x.stride = STRIDE_X;
    pars.y.N  = GRIDPOINTS_Y;
    pars.y.a  = SPATIAL_LOWER_BOUND_Y;
    pars.y.b  = SPATIAL_UPPER_BOUND_Y;
    pars.y.k  = TWOPI * I / (pars.y.b - pars.y.a);
    pars.y.k2 = -TWOPI * TWOPI / ((pars.y.b - pars.y.a) * (pars.y.b - pars.y.a));
    pars.y.stride = STRIDE_Y;
    pars.z.N  = GRIDPOINTS_Z;
    pars.z.a  = SPATIAL_LOWER_BOUND_Z;
    pars.z.b  = SPATIAL_UPPER_BOUND_Z;
    pars.z.k  = TWOPI * I / (pars.z.b - pars.z.a);
    pars.z.k2 = -TWOPI * TWOPI / ((pars.z.b - pars.z.a) * (pars.z.b - pars.z.a));
    pars.z.stride = STRIDE_Z;

    pars.x.outN = (pars.x.N + pars.x.stride - 1) / pars.x.stride;
    pars.y.outN = (pars.y.N + pars.y.stride - 1) / pars.y.stride;
    pars.z.outN = (pars.z.N + pars.z.stride - 1) / pars.z.stride;
    pars.outN = pars.x.outN * pars.y.outN * pars.z.outN;

    // set the number of dimensions according to gridpoints in each direction
    pars.dim = 3;
    if (pars.z.N == 1) {
        pars.dim = 2;
        if (pars.y.N == 1) {
            pars.dim = 1;
        }
    }

    // due to the memory usage of fftw, we need different upper bounds in for
    // loops depending on  the dimension, (the N gridpoints from the last
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

    // the total number of scalars evolved in the integration routine depends on
    // if and how we evolve psi
    pars.N = pars.x.N * pars.y.N * pars.z.N;
    pars.Nall = 4 * pars.N + 2;
    #if PSI_METHOD == PSI_ELLIPTIC
    pars.Ntot = 2 * pars.N + 1;
    #elif PSI_METHOD == PSI_PARABOLIC
    pars.Ntot = 3 * pars.N + 2;
    #elif PSI_METHOD == PSI_HYPERBOLIC
    pars.Ntot = 4 * pars.N + 2;
    #endif

    pars.t.dt = DELTA_T;
    pars.t.t  = INITIAL_TIME;
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
    pars.file.bins_powspec = POWER_SPECTRUM_BINS;
    INFO(printf("Initialized parameters using %zu dimension(s).\n\n",
            pars.dim));
}

// allocate memory for all external variables
void allocate_external()
{
    const size_t Nx   = pars.x.N;
    const size_t Ny   = pars.y.N;
    const size_t Nz   = pars.z.N;
    const size_t N    = pars.N;
    const size_t Nall = pars.Nall;
    const size_t outN = pars.outN;
    const size_t bins = pars.file.bins_powspec;
    const size_t buf_size = pars.file.buf_size;

    grid       = malloc((Nx + Ny + Nz) * sizeof *grid);
    field      = fftw_malloc(Nall * sizeof *field);
    field_new  = fftw_malloc(Nall * sizeof *field_new);
    dfield     = fftw_malloc(Nall * sizeof *dfield);
    dfield_new = fftw_malloc(Nall * sizeof *dfield_new);
    time_buf   = calloc(buf_size, sizeof *time_buf);
    a_buf      = calloc(buf_size, sizeof *a_buf);
    rho        = fftw_malloc(N * sizeof *rho);
    #if PSI_METHOD == PSI_HYPERBOLIC
    pressure   = fftw_malloc(N * sizeof *pressure);
    #endif
    pow_spec   = calloc(bins, sizeof *pow_spec);

    #ifdef OUTPUT_PHI
    phi_buf       = calloc(buf_size * outN, sizeof *phi_buf);
    #endif
    #ifdef OUTPUT_DPHI
    dphi_buf      = calloc(buf_size * outN, sizeof *dphi_buf);
    #endif
    #ifdef OUTPUT_PHI_MEAN
    phi_mean_buf  = calloc(buf_size, sizeof *phi_mean_buf);
    #endif
    #ifdef OUTPUT_PHI_VARIANCE
    phi_var_buf   = calloc(buf_size, sizeof *phi_var_buf);
    #endif
    #ifdef OUTPUT_DPHI_MEAN
    dphi_mean_buf = calloc(buf_size, sizeof *dphi_mean_buf);
    #endif
    #ifdef OUTPUT_DPHI_VARIANCE
    dphi_var_buf  = calloc(buf_size, sizeof *dphi_var_buf);
    #endif
    #ifdef OUTPUT_PSI
    psi_buf       = fftw_malloc(buf_size * outN * sizeof *psi_buf);
    #endif
    #ifdef OUTPUT_DPSI
    dpsi_buf      = fftw_malloc(buf_size * outN * sizeof *dpsi_buf);
    #endif
    #ifdef OUTPUT_PSI_MEAN
    psi_mean_buf  = calloc(buf_size, sizeof *psi_mean_buf);
    #endif
    #ifdef OUTPUT_PSI_VARIANCE
    psi_var_buf   = calloc(buf_size, sizeof *psi_var_buf);
    #endif
    #ifdef OUTPUT_DPSI_MEAN
    dpsi_mean_buf = calloc(buf_size, sizeof *dpsi_mean_buf);
    #endif
    #ifdef OUTPUT_DPSI_VARIANCE
    dpsi_var_buf  = calloc(buf_size, sizeof *dpsi_var_buf);
    #endif
    #ifdef OUTPUT_RHO
    rho_buf       = fftw_malloc(buf_size * outN * sizeof *rho_buf);
    #endif
    #ifdef OUTPUT_RHO_MEAN
    rho_mean_buf  = calloc(buf_size, sizeof *rho_mean_buf);
    #endif
    #ifdef OUTPUT_RHO_VARIANCE
    rho_var_buf   = calloc(buf_size, sizeof *rho_var_buf);
    #endif
    #ifdef OUTPUT_POWER_SPECTRUM
    pow_spec_buf  = calloc(buf_size * bins, sizeof *pow_spec_buf);
    #endif

    const size_t M = pars.x.M * pars.y.M * pars.z.M;
    #ifdef ENABLE_FFT_FILTER
    filter = fftw_malloc(M * sizeof *filter);
    #endif

    // default arrays to save coefficients of real to complex transforms
    // see fftw3 documentation and Mxyz for this
    tmp.phic  = fftw_malloc(M * sizeof *tmp.phic);
    tmp.xphic = fftw_malloc(M * sizeof *tmp.xphic);
    tmp.yphic = fftw_malloc(M * sizeof *tmp.yphic);
    tmp.zphic = fftw_malloc(M * sizeof *tmp.zphic);
    tmp.psic  = fftw_malloc(M * sizeof *tmp.psic);
    tmp.dpsic = fftw_malloc(M * sizeof *tmp.dpsic);
    tmp.fc    = fftw_malloc(M * sizeof *tmp.fc);
    tmp.deltarhoc  = fftw_malloc(M * sizeof *tmp.deltarhoc);

    // general purpose double memory blocks for temporary use
    tmp.xphi = fftw_malloc(N * sizeof *tmp.xphi);
    tmp.yphi = fftw_malloc(N * sizeof *tmp.yphi);
    tmp.zphi = fftw_malloc(N * sizeof *tmp.zphi);
    tmp.grad = fftw_malloc(N * sizeof *tmp.grad);
    tmp.lap  = fftw_malloc(N * sizeof *tmp.lap);
    tmp.f    = fftw_malloc(N * sizeof *tmp.f);
    tmp.deltarho = fftw_malloc(N * sizeof *tmp.deltarho);

    if (!(grid && field && field_new && dfield && dfield_new &&
        rho && pow_spec && tmp.phic  && tmp.xphic && tmp.yphic && tmp.zphic &&
        tmp.xphi && tmp.yphi && tmp.zphi && tmp.grad && tmp.lap && tmp.psic  &&
        tmp.fc && tmp.deltarhoc && tmp.dpsic && tmp.f && tmp.deltarho)) {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    INFO(puts("Allocated memory for external variables.\n"));
}

// make the N fourier spectral gridpoints for the computational domain
void mk_grid()
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

    if (Nx < 1 || Ny < 1 || Nz < 1) {
        fputs("Need positive number of gridpoints\n", stderr);
        exit(EXIT_FAILURE);
    }

    // set up the grid points
    #pragma omp parallel for
    for (size_t i = 0; i < Nx; ++i) {
        grid[i] = ax + (bx - ax) * i / Nx;
    }
    #pragma omp parallel for
    for (size_t j = Nx; j < Nx+Ny; ++j) {
        grid[j] = ay + (by - ay) * (j - Nx) / Ny;
    }
    #pragma omp parallel for
    for (size_t k = Nx + Ny; k < Nx + Ny+ Nz; ++k) {
        grid[k] = az + (bz - az) * (k - Nx - Ny) / Nz;
    }

    INFO(puts("Constructed gridpoints.\n"));
}

// create the fftw plans, IMPORTANT: create BEFORE initializing arrays, because
// setting up the plans destroys the arrays!
void mk_fftw_plans()
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;

    #ifdef SHOW_TIMING_INFO
    fftw_time_plan -= get_wall_time();
    #endif
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
    #ifdef SHOW_TIMING_INFO
    fftw_time_plan += get_wall_time();
    #endif
    INFO(puts("Created fftw plans.\n"));
}

// setup initial conditions for the field
void mk_initial_conditions()
{
    #if INITIAL_CONDITIONS == IC_FROM_H5_FILE
    h5_read_timeslice();
    #elif INITIAL_CONDITIONS == IC_FROM_DAT_FILE
    read_initial_data();
    field[2 * pars.N] = A_INITIAL;
    field[2 * pars.N + 1] = 0.0;
        #if PSI_METHOD != PSI_ELLIPTIC
        mk_initial_psi();
        #endif
    #elif INITIAL_CONDITIONS == IC_FROM_BUNCH_DAVIES
    //TODO: need correct values at end of inflation here
    /* const double phi0 = -0.008283249078352; */
    /* const double dphi0 = -5.738421024284034e-05; */
    /* const double hubble = sqrt(0.5 * (dphi0 * dphi0 + MASS * MASS * phi0 * phi0) / 3.0); */
    //defrost
    /* double phi0 = 2.339383796213256; */
    /* double dphi0 = -2.736358272992573; */
    /* double hubble = 1.934897490588959; */
    //pspectre defrost style
    const double phi0 = 1.0093430384226378929425913902459/sqrt(8.0 * PI);
    const double dphi0 = 0.0;
    const double hubble = sqrt(4.0 * PI * (dphi0 * dphi0 + MASS * MASS * phi0 * phi0) / 3.0);
    mk_bunch_davies(field, hubble, phi0, -0.25);
    mk_bunch_davies(field + pars.N, hubble, dphi0, 0.25);
    field[2 * pars.N] = A_INITIAL;
    field[2 * pars.N + 1] = 0.0;
        #if PSI_METHOD != PSI_ELLIPTIC
        mk_initial_psi();
        #endif
    #elif INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t N = pars.N;
    size_t osx, osy;
    double x, y, z;

    const size_t Nmodes = 16;
    // random phases
    srand(SEED);
    double *theta = calloc(Nmodes, sizeof *theta);
    for (size_t i = 0; i < Nmodes; ++i) {
        theta[i] = TWOPI * (double)rand() / (double)RAND_MAX;
    }

    // initialize the scalar field and its temporal derivative
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

    free(theta);

    field[2 * N] = A_INITIAL;
    field[2 * N + 1] = 0.0;

    // initialize psi
    #if PSI_METHOD != PSI_ELLIPTIC
    mk_initial_psi();
    #endif
    #endif

    INFO(puts("Initialized fields on first time slice.\n"));
}

// create mask for the fourier filtering
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
                tmp = 1.0;
                if (i != 0) {
                    tmp = filter_window(2.0 *
                        (i > Nx / 2 ? (int)Nx - (int)i : i) / (double) Nx);
                }
                if (pars.dim > 1) {
                    if (j != 0) {
                        tmp *= filter_window(2.0 *
                            (j > Ny / 2 ? (int)Ny - (int)j : j) / (double) Ny);
                    }
                    if (pars.dim > 2 && k != 0) {
                        tmp *= filter_window(2.0 * k / (double) Nz);
                    }
                }
                filter[osy + k] = tmp / (double) N;
            }
        }
    }
}

// the cutoff function for filtering, use either two thirds or fourier smoothing
inline double filter_window(const double x)
{
    // fourier smoothing
    return exp(-36.0 * pow(x, 36));

    // two thirds rule
    // return x < 2.0/3.0 ? x : 0.0;

    // miscellaneous, have been used in testing
    /* return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. ); */
    /* return exp(1. + 1. / ( pow(x, 8) - 1. )); */
    /* return 0.5 * ( 1. + cos( pow(x, 8) * PI ) ); */
    /* return 0.0; */
}

// initial values of the scalar field, make sure its periodic
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
    const double scale= 1.0e4;
    const double mean = 0.0510864;
    const double amplitude = -3.743790000000000e-07 * scale;

    // compare_2, pos= 1
    /* double mean = 5.0; */
    /* double amplitude = 0.01; */

    /* double k = 1.0/6.0e3; */
    const double k = 1.0;
    if (pars.dim == 1) {
        return mean + amplitude * cos(k * x);
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

// initial values of the time deriv. of the scalar field, make sure its periodic
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

    const double scale = 1.0e4;
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

void free_and_destroy_all()
{
    h5_close(pars.file.id);
    destroy_and_cleanup_fftw();
    free_external();
}

// destroy the fftw plans and call cleanup for internal fftw3 cleanup
void destroy_and_cleanup_fftw()
{
    fftw_destroy_plan(p_fw);
    fftw_destroy_plan(p_bw);
    fftw_cleanup_threads();
    INFO(puts("Destroyed fftw plans.\n"));
}

// free memory of all global variables
void free_external()
{
    free(grid);
    fftw_free(field);
    fftw_free(field_new);
    fftw_free(dfield);
    fftw_free(dfield_new);
    free(rho);
    #if PSI_METHOD == PSI_HYPERBOLIC
    free(pressure);
    #endif
    free(pow_spec);
    free(time_buf);
    free(a_buf);
    #ifdef OUTPUT_PHI
    free(phi_buf);
    #endif
    #ifdef OUTPUT_DPHI
    free(dphi_buf);
    #endif
    #ifdef OUTPUT_PHI_MEAN
    free(phi_mean_buf);
    #endif
    #ifdef OUTPUT_PHI_VARIANCE
    free(phi_var_buf);
    #endif
    #ifdef OUTPUT_DPHI_MEAN
    free(dphi_mean_buf);
    #endif
    #ifdef OUTPUT_DPHI_VARIANCE
    free(dphi_var_buf);
    #endif
    #ifdef OUTPUT_PSI
    free(psi_buf);
    #endif
    #ifdef OUTPUT_DPSI
    free(dpsi_buf);
    #endif
    #ifdef OUTPUT_PSI_MEAN
    free(psi_mean_buf);
    #endif
    #ifdef OUTPUT_PSI_VARIANCE
    free(psi_var_buf);
    #endif
    #ifdef OUTPUT_DPSI_MEAN
    free(dpsi_mean_buf);
    #endif
    #ifdef OUTPUT_DPSI_VARIANCE
    free(dpsi_var_buf);
    #endif
    #ifdef OUTPUT_RHO
    free(rho_buf);
    #endif
    #ifdef OUTPUT_RHO_MEAN
    free(rho_mean_buf);
    #endif
    #ifdef OUTPUT_RHO_VARIANCE
    free(rho_var_buf);
    #endif
    #ifdef OUTPUT_POWER_SPECTRUM
    free(pow_spec_buf);
    #endif
    #ifdef ENABLE_FFT_FILTER
    free(filter);
    #endif
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
    fftw_free(tmp.fc);
    fftw_free(tmp.deltarhoc);
    fftw_free(tmp.dpsic);
    fftw_free(tmp.f);
    fftw_free(tmp.deltarho);
    INFO(puts("Freed external variables.\n"));
}

void mk_bunch_davies(double *f, const double H, const double homo,
        const double gamma)
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    if (Nx != Ny || Nx != Nz || Ny != Nz) {
        fputs("Bunch Davies vacuum works only for Nx = Ny = Nz.\n", stderr);
        exit(EXIT_FAILURE);
    }
    const size_t N  = pars.N;
    const size_t nn = Nx / 2 + 1;
    const size_t os = 16;
    const size_t nos = Nx * os * os;
    const double dx = (pars.x.b - pars.x.a) / Nx;
    const double dxos = dx / os;
    const double dk = TWOPI / (pars.x.b - pars.x.a);
    const double dkos = 0.5 * dk / os;
    //TODO: pspectre uses kcutpspectre = 2 * kcutdefrost
    const double kcut2 = 0.25 * nn * nn * dk * dk;
    /* const double kcut2 = 0.01 * nn * nn * dk * dk; */
    const double meff2 = MASS * MASS - 2.25 * H * H;
    /* const double norm = 0.5 / (N * sqrt(TWOPI * pow(dk, 3))) * (dkos / dxos); */
    const double norm = 0.5 / (N * sqrt(TWOPI * pow(dk, 3)) *
            (2.e5/sqrt(8*PI))) * (dkos / dxos);

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
                /* tmp.phic[osy + k] *= box_muller(); */
                tmp.phic[osy + k] *= box_muller() / (8.0 * PI);
            }
        }
    }

    tmp.phic[0] = homo;
    fftw_execute_dft_c2r(p_bw, tmp.phic, f);
}

inline complex box_muller()
{
    const double u1 = (double)rand() / (double)RAND_MAX;
    const double u2 = (double)rand() / (double)RAND_MAX;
    return sqrt(-2 * log(u1)) * cexp(TWOPI * u2 * 1i);
}

// -------------------------printing function-----------------------------------
// for debugging mostly
void print_vector(const double *vector, const size_t N)
{
    for (size_t i = 0; i < N; i++) {
        printf("%f\n", vector[i]);
    }
}
