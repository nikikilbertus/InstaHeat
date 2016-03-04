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

void allocate_and_initialize_all() {
    initialize_threading();
    initialize_parameters();
    allocate_external();
    mk_grid();
    mk_fftw_plans();
    mk_initial_conditions();
    h5_create_empty_by_path(DATAPATH);
}

void initialize_threading() {
    int threadnum, threadinit;
    threadinit = fftw_init_threads();
    if (threadinit == 0)
    {
        fputs("Could not initialize fftw threads.\n", stderr);
        exit(EXIT_FAILURE);
    }
    threadnum = THREAD_NUMBER <= 0 ? omp_get_max_threads() : THREAD_NUMBER;
    omp_set_num_threads(threadnum);
    fftw_plan_with_nthreads(threadnum);
    RUNTIME_INFO(printf("Running omp & fftw with %d thread(s)\n\n", threadnum));
}

/**
 *  initialize the values in the paramters_t pars variable, mostly from defines
 *  in main.h; using the struct gives more flexibility than using the defines
 *  throughout the code
 */
void initialize_parameters() {
    pars.x.N  = GRIDPOINTS_X;
    pars.x.a  = SPATIAL_LOWER_BOUND_X;
    pars.x.b  = SPATIAL_UPPER_BOUND_X;
    pars.x.k  = 2. * PI * I / (pars.x.b - pars.x.a);
    pars.x.k2 = -4. * PI * PI / ((pars.x.b - pars.x.a) * (pars.x.b - pars.x.a));
    pars.x.stride = STRIDE_X;
    pars.y.N  = GRIDPOINTS_Y;
    pars.y.a  = SPATIAL_LOWER_BOUND_Y;
    pars.y.b  = SPATIAL_UPPER_BOUND_Y;
    pars.y.k  = 2. * PI * I / (pars.y.b - pars.y.a);
    pars.y.k2 = -4. * PI * PI / ((pars.y.b - pars.y.a) * (pars.y.b - pars.y.a));
    pars.y.stride = STRIDE_Y;
    pars.z.N  = GRIDPOINTS_Z;
    pars.z.a  = SPATIAL_LOWER_BOUND_Z;
    pars.z.b  = SPATIAL_UPPER_BOUND_Z;
    pars.z.k  = 2. * PI * I / (pars.z.b - pars.z.a);
    pars.z.k2 = -4. * PI * PI / ((pars.z.b - pars.z.a) * (pars.z.b - pars.z.a));
    pars.z.stride = STRIDE_Z;

    pars.x.outN = (pars.x.N + pars.x.stride - 1) / pars.x.stride;
    pars.y.outN = (pars.y.N + pars.y.stride - 1) / pars.y.stride;
    pars.z.outN = (pars.z.N + pars.z.stride - 1) / pars.z.stride;
    pars.outN = pars.x.outN * pars.y.outN * pars.z.outN;

    // set the number of dimensions according to gridpoints in each direction
    pars.dim = 3;
    if (pars.z.N == 1)
    {
        pars.dim = 2;
        if (pars.y.N == 1)
        {
            pars.dim = 1;
        }
    }

    // due to the memory usage of fftw, we need different upper bounds in for
    // loops depending on  the dimension, (the N gridpoints from the last
    // dimension are transformed to floor(N/2)+1 points in fourier space)
    switch (pars.dim)
    {
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
    #if PSI_METHOD == PSI_ELLIPTIC
    pars.Ntot = 2 * N + 1;
    #elif PSI_METHOD == PSI_PARABOLIC
    pars.Ntot = 3 * N + 1;
    #elif PSI_METHOD == PSI_HYPERBOLIC
    pars.Ntot = 4 * N + 1;
    #endif

    pars.t.dt = DELTA_T;
    pars.t.t  = INITIAL_TIME;
    pars.t.ti = INITIAL_TIME;
    pars.t.tf = FINAL_TIME;
    pars.t.Nt = ceil((pars.t.tf - pars.t.ti) / pars.t.dt) + 1;
    if (pars.t.Nt > MAX_STEPS)
    {
        fputs("Exeeding MAX_STEPS, decrease DELTA_T.\n", stderr);
        exit(EXIT_FAILURE);
    }

    pars.file.index = 0;
    pars.file.buf_size = WRITE_OUT_BUFFER_NUMBER;
    pars.file.skip = TIME_STEP_SKIPS;
    pars.file.bins_powspec = POWER_SPECTRUM_BINS;
    RUNTIME_INFO(printf("Initialized parameters using %zu dimension(s).\n\n",
            pars.dim));
}

// allocate memory for all external variables
void allocate_external() {
    size_t Nx   = pars.x.N;
    size_t Ny   = pars.y.N;
    size_t Nz   = pars.z.N;
    size_t N    = pars.N;
    size_t Ntot = pars.Ntot;
    size_t outN = pars.outN;
    size_t bins = pars.file.bins_powspec;
    size_t buf_size = pars.file.buf_size;

    grid       = malloc((Nx + Ny + Nz) * sizeof *grid);
    field      = fftw_malloc(Ntot * sizeof *field);
    field_new  = fftw_malloc(Ntot * sizeof *field_new);
    dfield     = fftw_malloc(Ntot * sizeof *dfield);
    dfield_new = fftw_malloc(Ntot * sizeof *dfield_new);
    time_buf   = calloc(buf_size, sizeof *time_buf);
    a_buf      = calloc(buf_size, sizeof *a_buf);
    rho        = fftw_malloc(N * sizeof *rho);
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

    // default arrays to save coefficients of real to complex transforms
    // see fftw3 documentation and Mxyz for this
    size_t M = pars.x.M * pars.y.M * pars.z.M;
    tmp.phic  = fftw_malloc(M * sizeof *tmp.phic);
    tmp.xphic = fftw_malloc(M * sizeof *tmp.xphic);
    tmp.yphic = fftw_malloc(M * sizeof *tmp.yphic);
    tmp.zphic = fftw_malloc(M * sizeof *tmp.zphic);
    tmp.psic  = fftw_malloc(M * sizeof *tmp.psic);
    tmp.fc    = fftw_malloc(M * sizeof *tmp.fc);
    tmp.deltarhoc  = fftw_malloc(M * sizeof *tmp.deltarhoc);
    tmp.dpsic = fftw_malloc(M * sizeof *tmp.dpsic);

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
        tmp.fc && tmp.deltarhoc && tmp.dpsic && tmp.f && tmp.deltarho))
    {
        fputs("Allocating memory failed.\n", stderr);
        exit(EXIT_FAILURE);
    }
    RUNTIME_INFO(puts("Allocated memory for external variables.\n"));
}

// make the N fourier spectral gridpoints for the computational domain
void mk_grid() {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    double ax = pars.x.a;
    double bx = pars.x.b;
    double ay = pars.y.a;
    double by = pars.y.b;
    double az = pars.z.a;
    double bz = pars.z.b;

    if (Nx < 1 || Ny < 1 || Nz < 1)
    {
        fputs("Need positive number of gridpoints\n", stderr);
        exit(EXIT_FAILURE);
    }

    // set up the grid points
    #pragma omp parallel for
    for (size_t i = 0; i < Nx; ++i)
    {
        grid[i] = ax + (bx - ax) * i / Nx;
    }
    #pragma omp parallel for
    for (size_t j = Nx; j < Nx+Ny; ++j)
    {
        grid[j] = ay + (by - ay) * (j - Nx) / Ny;
    }
    #pragma omp parallel for
    for (size_t k = Nx + Ny; k < Nx + Ny+ Nz; ++k)
    {
        grid[k] = az + (bz - az) * (k - Nx - Ny) / Nz;
    }

    RUNTIME_INFO(puts("Constructed gridpoints.\n"));
}

// create the fftw plans, IMPORTANT: create BEFORE initializing arrays, because
// setting up the plans destroys the arrays!
void mk_fftw_plans() {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;

    #ifdef SHOW_TIMING_INFO
    fftw_time_plan -= get_wall_time();
    #endif
    switch (pars.dim)
    {
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
    RUNTIME_INFO(puts("Created fftw plans.\n"));
}

// setup initial conditions for the field
void mk_initial_conditions() {
    #if INITIAL_CONDITIONS == IC_FROM_H5_FILE
    h5_read_timeslice(pars.t.ti, field);
    #elif INITIAL_CONDITIONS == IC_FROM_DAT_FILE
    read_initial_data();
    mk_initial_psi();
    #elif INITIAL_CONDITIONS == IC_FROM_INTERNAL_FUNCTION
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t osx, osy;
    double x, y, z;

    size_t Nmodes = 16;
    /* double *ks = calloc(pars.dim * Nmodes, sizeof *ks); */
    /* for (size_t i = 0; i < pars.dim * Nmodes; ++i) */
    /* { */
    /*     ks[i] = 1.0; */
    /* } */
    /* for (size_t i = 3; i < Nmodes; i += 3) */
    /* { */
    /*     ks[i] *= -1.0; */
    /* } */

    // random phases
    srand(SEED);
    double *theta = calloc(Nmodes, sizeof *theta);
    for (size_t i = 0; i < Nmodes; ++i)
    {
        theta[i] = 2.0 * PI * (double)rand() / (double)RAND_MAX;
    }

    // initialize the scalar field and its temporal derivative
    for (size_t i = 0; i < Nx; ++i)
    {
        x = grid[i];
        osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j)
        {
            y = grid[Nx + j];
            osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k)
            {
                z = grid[Nx + Ny + k];
                field[osy + k] = phi_init(x, y, z, theta);
                field[N + osy + k] = dphi_init(x, y, z, theta);
            }
        }
    }

    // initialize a
    field[3 * N] = A_INITIAL;
    free(theta);
    #endif

    // initialize psi
    RUNTIME_INFO(puts("Initialized the field and its temporal derivative.\n"));
}

void read_initial_data() {
    size_t N = pars.N;

    FILE *file = fopen(INITIAL_DATAPATH, "r");
    if (!file)
    {
        fputs("Could not read initial data file.\n", stderr);
        exit(EXIT_FAILURE);
    }

    //TODO: adjust to actual file format, this is just a dummy
    for (size_t i = 0; i < N; ++i)
    {
        fscanf(file, "%lf", &field[i]);
    }
    fclose(file);
}

// initial values of the scalar field, make sure its periodic
double phi_init(double x, double y, double z, double *ph) {
    // localized for higgs metastability potential
    /* double phi0 = 0.04; */
    /* double lambda = 20.0; */
    /* if (pars.dim == 1) */
    /* { */
    /*     return phi0 * 0.5 * (1.0 + cos(x)) * exp(-lambda * x * x); */
    /* } */
    /* else if (pars.dim == 2) */
    /* { */
    /*     return phi0 * 0.25 * (1.0 + cos(x)) * (1.0 + cos(y)) * */
    /*         exp(-lambda * (x * x + y * y)); */
    /* } */
    /* else */
    /* { */
    /*     return phi0 * 0.125 * (1.0 + cos(x)) * (1.0 + cos(y)) * (1.0 + cos(z)) * */
    /*        exp(-lambda * (x * x + y * y + z * z)); */
    /* } */

    // some simple waves for notch or step potential simulations
    /* double frac = 0.4; // vary the ratio between \phi_0 and \delta \phi */
    /* double phi0 = 0.73 * frac; // only vary if you know exactly why */
    /* double deltaphi = phi0 / frac; */
    /* if (pars.dim == 1) */
    /* { */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1])); */
    /* } */
    /* else if (pars.dim == 2) */
    /* { */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1]) + */
    /*                  cos(1.0 * y + ph[2]) + cos(-1.0 * y + ph[3])); */
    /* } */
    /* else */
    /* { */
    /*     return phi0 + deltaphi * */
    /*                 (cos(1.0 * x + ph[0]) + cos(-1.0 * x + ph[1]) + */
    /*                  cos(1.0 * y + ph[2]) + cos(-1.0 * y + ph[3]) + */
    /*                  cos(1.0 * z + ph[4]) + cos(-1.0 * z + ph[5])); */
    /* } */

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
    double scale= 1.0e4;
    double mean = 0.0510864;
    double amplitude = -3.743790000000000e-07 * scale;

    // compare_2, pos= 1
    /* double mean = 5.0; */
    /* double amplitude = 0.01; */

    /* double k = 1.0/6.0e3; */
    double k = 1.0;
    if (pars.dim == 1)
    {
        return mean + amplitude * cos(k * x);
        /* return mean - amplitude * wrapped_gaussian(x, y, z); */
    }
    else if (pars.dim == 2)
    {
        /* return mean + amplitude * cos(x + y + ph[0]); */
        return mean - amplitude * wrapped_gaussian(x, y, z);
    }
    else
    {
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
double dphi_init(double x, double y, double z, double *ph) {
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

    double scale = 1.0e4;
    double mean = 3.255190000000000e-04;
    double amplitude = 1.742130000000000e-08 * scale;

    /* double mean = -0.00806088; */
    /* double amplitude = -1.134420000000000e-20; */

    if (pars.dim == 1)
    {
        return (mean + amplitude * cos(x)) * MASS / MASS_KARSTEN;
    }
    else if (pars.dim == 2)
    {
        return (mean + amplitude * cos(x + y + ph[0])) *
            MASS / MASS_KARSTEN;
    }
    else
    {
        return (mean + amplitude *
            (cos(x + y + z + ph[0]) + cos(-x + y + z + ph[1]) +
             cos(x - y + z + ph[2]) + cos(x + y - z + ph[3]) +
             cos(2.0 * x + y + z + ph[4]) + cos(x + 2.0 * y + z + ph[5]))) *
             MASS / MASS_KARSTEN;
        /* return (mean + amplitude * cos(x + y + z + ph[0])) * MASS / 1.0e-2; */
    }

    /* return -0.089318193; // somewhere at end of 50 e-fold inflation */
}

double wrapped_gaussian(double x, double y, double z) {
    double s = 0.5;
    double res = 0.0;
    if (pars.dim == 1)
    {
        size_t max = 32;
        for (size_t i = 1; i <= max; ++i)
        {
            res += exp(-0.5 * i * i * s * s) * (cos(i * x) + pow(-1.0, i + 1)) /
                (2.0 * PI);
        }
    }
    if (pars.dim == 2)
    {
        size_t max = 8;
        for (size_t i = 1; i <= max; ++i)
        {
            for (size_t j = 1; j <= max; ++j)
            {
                res += exp(-0.5 * (i * i + j * j) * s * s) *
                    (cos(i * x) + pow(-1.0, i + 1)) *
                    (cos(j * y) + pow(-1.0, j + 1)) /
                    (2.0 * PI);
            }
        }
    }
    if (pars.dim == 3)
    {
        size_t max = 16;
        for (size_t i = 1; i <= max; ++i)
        {
            for (size_t j = 1; j <= max; ++j)
            {
                for (size_t k = 1; k <= max; ++k)
                {
                    res += exp(-0.5 * (i * i + j * j + k * k) * s * s) *
                        (cos(i * x) + pow(-1.0, i + 1)) *
                        (cos(j * y) + pow(-1.0, j + 1)) *
                        (cos(k * z) + pow(-1.0, k + 1)) /
                        (2.0 * PI);
                }
            }
        }
    }
    return res;
}

void mk_initial_psi() {
    size_t N = pars.N;
    size_t N2 = 2 * N;
    size_t N3 = 3 * N;

    #pragma omp parallel for
    for (size_t i = N2; i < N3; ++i)
    {
        field[i] = 0.0;
    }

    mk_gradient_squared_and_laplacian(field);
    mk_rho(field);
    mk_psi(field);
}

void free_and_destroy_all() {
    h5_close(pars.file.id);
    destroy_and_cleanup_fftw();
    free_external();
}

// destroy the fftw plans and call cleanup for internal fftw3 cleanup
void destroy_and_cleanup_fftw() {
    fftw_destroy_plan(p_fw);
    fftw_destroy_plan(p_bw);
    fftw_cleanup_threads();
    RUNTIME_INFO(puts("Destroyed fftw plans.\n"));
}

// free memory of all global variables
void free_external() {
    free(grid);
    fftw_free(field);
    fftw_free(field_new);
    fftw_free(dfield);
    fftw_free(dfield_new);
    free(rho);
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
    RUNTIME_INFO(puts("Freed external variables.\n"));
}

// -------------------------printing function-----------------------------------
// for debugging mostly
void print_vector(const double *vector, const size_t N) {
    for (size_t i = 0; i < N; i++)
    {
        printf("%f\n", vector[i]);
    }
}
