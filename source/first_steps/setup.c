#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "setup.h"
#include "main.h"

void allocate_and_initialize_all() {
    initialize_threading();
    initialize_parameters();
    allocate_external();
    mk_grid();
    mk_fftw_plans();
    mk_initial_conditions();
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

    pars.N = pars.x.N * pars.y.N * pars.z.N;

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
    size_t N2   = 2 * N;
    size_t Ntot = N2 + 1;
    size_t buf_size = pars.file.buf_size;
    size_t bins = pars.file.bins_powspec;
    size_t outN = pars.outN;

    grid         = malloc((Nx + Ny + Nz) * sizeof *grid);
    // note that the field contains the scalar field, its time deriv. and a
    field        = fftw_malloc(Ntot * sizeof *field);
    field_new    = fftw_malloc(Ntot * sizeof *field_new);
    dfield       = fftw_malloc(Ntot * sizeof *dfield);
    dfield_new   = fftw_malloc(Ntot * sizeof *dfield_new);
    // this buffer only holds the scalar field (not the deriv. or a)
    field_buf    = calloc(buf_size * outN, sizeof *field_buf);
    psi          = fftw_malloc(N * sizeof *psi);
    dpsi         = fftw_malloc(N * sizeof *dpsi);
    psi_buf      = fftw_malloc(buf_size * outN * sizeof *psi_buf);
    time_buf     = calloc(buf_size, sizeof *time_buf);
    f_a_buf      = calloc(buf_size, sizeof *f_a_buf);
    rho          = fftw_malloc(N * sizeof *rho);
    rho_buf      = fftw_malloc(buf_size * outN * sizeof *rho_buf);
    pow_spec     = calloc(bins, sizeof *pow_spec);
    pow_spec_buf = calloc(buf_size * bins, sizeof *pow_spec_buf);

    // default arrays to save coefficients of real to complex transforms
    // see fftw3 documentation and Mxyz for this
    // TODO: maybe wrap this allocation in extra method (slight differences!)
    size_t M = pars.x.M * pars.y.M * pars.z.M;
    tmp_phi.c   = fftw_malloc(M * sizeof *tmp_phi.c);
    tmp_phi.cx  = fftw_malloc(M * sizeof *tmp_phi.cx);
    tmp_phi.cy  = fftw_malloc(M * sizeof *tmp_phi.cy);
    tmp_phi.cz  = fftw_malloc(M * sizeof *tmp_phi.cz);
    tmp_psi.c   = fftw_malloc(M * sizeof *tmp_psi.c);
    tmp_psi.cx  = fftw_malloc(M * sizeof *tmp_psi.cx);
    tmp_psi.cy  = fftw_malloc(M * sizeof *tmp_psi.cy);
    tmp_psi.cz  = fftw_malloc(M * sizeof *tmp_psi.cz);

    // general purpose double memory blocks for temporary use
    tmp_phi.dx   = fftw_malloc(Ntot * sizeof *tmp_phi.dx); // used in dopri853 (dense)
    tmp_phi.dy   = fftw_malloc(N * sizeof *tmp_phi.dy);
    tmp_phi.dz   = fftw_malloc(N * sizeof *tmp_phi.dz);
    tmp_phi.grad = fftw_malloc(N * sizeof *tmp_phi.grad);
    tmp_phi.lap  = fftw_malloc(N * sizeof *tmp_phi.lap);
    tmp_psi.dx   = fftw_malloc(N * sizeof *tmp_psi.dx);
    tmp_psi.dy   = fftw_malloc(N * sizeof *tmp_psi.dy);
    tmp_psi.dz   = fftw_malloc(N * sizeof *tmp_psi.dz);
    tmp_psi.grad = fftw_malloc(N * sizeof *tmp_psi.grad);

    if (!(grid && field && field_new && dfield && dfield_new && field_buf &&
        psi && dpsi && time_buf && rho && rho_buf && pow_spec && pow_spec_buf &&
        tmp_phi.c  && tmp_phi.cx && tmp_phi.cy && tmp_phi.cz &&
        tmp_phi.dx && tmp_phi.dy && tmp_phi.dz && tmp_phi.grad && tmp_phi.lap &&
        tmp_psi.c  && tmp_psi.cx && tmp_psi.cy && tmp_psi.cz &&
        tmp_psi.dx && tmp_psi.dy && tmp_psi.dz && tmp_psi.grad))
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

    #ifdef DEBUG
    puts("x");
    print_vector(grid, Nx);
    puts("\n");
    puts("y");
    print_vector(grid + Nx, Ny);
    puts("\n");
    puts("z");
    print_vector(grid + Nx + Ny, Nz);
    puts("\n");
    #endif
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
            p_fw = fftw_plan_dft_r2c_1d(Nx, field, tmp_phi.c,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_1d(Nx, tmp_phi.c, field,
                    FFTW_DEFAULT_FLAG);
            break;
        case 2:
            p_fw = fftw_plan_dft_r2c_2d(Nx, Ny, field, tmp_phi.c,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_2d(Nx, Ny, tmp_phi.c, field,
                    FFTW_DEFAULT_FLAG);
            break;
        case 3:
            p_fw = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, tmp_phi.c,
                    FFTW_DEFAULT_FLAG);
            p_bw = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, tmp_phi.c, field,
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
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t osx, osy;
    double x, y, z;

    // random phases used for analysis of notch and step potential
    srand(SEED);
    double *theta = calloc(6, sizeof *theta);
    for (size_t i = 0; i < 6; ++i)
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
                field[N + osy + k] = dphi_init(x, y, z);
            }
        }
    }

    // initialize a
    field[2 * N] = 1.0;

    // console output for debugging
    #ifdef DEBUG
    puts("phi");
    print_vector(field, N);
    puts("\ndphi");
    print_vector(field + N, N);
    puts("\na");
    print_vector(field + 2 * N, 1);
    puts("\n");
    #endif
    free(theta);
    RUNTIME_INFO(puts("Initialized the field and its temporal derivative.\n"));
}

// initial values of the scalar field, make sure its periodic
double phi_init(const double x, const double y, const double z,
                                                const double *ph) {
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
    double mean = 1.0;
    double amplitude = 1.0e-5 * mean;
    double k = MASS / 6.0;
    if (pars.dim == 1)
    {
        return mean + amplitude * cos(k * x);
    }
    else if (pars.dim == 2)
    {
        return mean + amplitude * cos(x + ph[0]) * cos(y + ph[1]);
    }
    else
    {
        return mean + amplitude *
            cos(x + ph[0]) * cos(y + ph[1]) * cos(z + ph[2]);
    }
}

// initial values of the time deriv. of the scalar field, make sure its periodic
double dphi_init(const double x, const double y, const double z) {
    return 0.0;

    /* return -0.089318193; // somewhere at end of 50 e-fold inflation */
}

void free_and_destroy_all() {
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
    free(field_buf);
    free(psi);
    free(dpsi);
    free(psi_buf);
    free(time_buf);
    free(f_a_buf);
    free(rho);
    free(rho_buf);
    free(pow_spec);
    free(pow_spec_buf);
    fftw_free(tmp_phi.c);
    fftw_free(tmp_phi.cx);
    fftw_free(tmp_phi.cy);
    fftw_free(tmp_phi.cz);
    fftw_free(tmp_phi.dx);
    fftw_free(tmp_phi.dy);
    fftw_free(tmp_phi.dz);
    fftw_free(tmp_phi.grad);
    fftw_free(tmp_phi.lap);
    fftw_free(tmp_psi.c);
    fftw_free(tmp_psi.cx);
    fftw_free(tmp_psi.cy);
    fftw_free(tmp_psi.cz);
    fftw_free(tmp_psi.dx);
    fftw_free(tmp_psi.dy);
    fftw_free(tmp_psi.dz);
    fftw_free(tmp_psi.grad);
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
