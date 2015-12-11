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
	pars.x.N = GRIDPOINTS_X;
	pars.x.a = SPATIAL_LOWER_BOUND_X;
	pars.x.b = SPATIAL_UPPER_BOUND_X;
    pars.x.L = pars.x.b - pars.x.a;
    pars.y.N = GRIDPOINTS_Y;
    pars.y.a = SPATIAL_LOWER_BOUND_Y;
    pars.y.b = SPATIAL_UPPER_BOUND_Y;
    pars.y.L = pars.y.b - pars.y.a;
    pars.z.N = GRIDPOINTS_Z;
    pars.z.a = SPATIAL_LOWER_BOUND_Z;
    pars.z.b = SPATIAL_UPPER_BOUND_Z;
    pars.z.L = pars.z.b - pars.z.a;

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
    RUNTIME_INFO(puts("Initialized parameters.\n"));
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

    grid = malloc((Nx + Ny + Nz) * sizeof *grid);
    // note that the field contains the scalar field, its time deriv. and a
    field        = fftw_malloc(Ntot * sizeof *field);
    field_new    = fftw_malloc(Ntot * sizeof *field_new);
    dfield       = fftw_malloc(Ntot * sizeof *dfield);
    dfield_new   = fftw_malloc(Ntot * sizeof *dfield_new);
    // this buffer only holds the scalar field (not the deriv. or a) 
    field_buf    = calloc(buf_size * N, sizeof *field_buf);
    time_buf     = calloc(buf_size, sizeof *time_buf);
    f_a_buf      = calloc(buf_size, sizeof *f_a_buf);
    rho_buf      = calloc(buf_size, sizeof *rho_buf);
    pow_spec     = calloc(bins, sizeof *pow_spec);
    pow_spec_buf = calloc(buf_size * bins, sizeof *pow_spec_buf);

    // default arrays to save coefficients of real to complex transforms
    size_t ncz = Nz / 2 + 1; // see fftw3 documentation for this
    cfftw_tmp   = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp);
    cfftw_tmp_x = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_x);
    cfftw_tmp_y = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_y);
    cfftw_tmp_z = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_z);

    // general purpose double memory blocks for temporary use
    dtmp_x = fftw_malloc(Ntot * sizeof *dtmp_x); // used in dopri853 (dense)
    dtmp_y = fftw_malloc(N * sizeof *dtmp_y);
    dtmp_z = fftw_malloc(N * sizeof *dtmp_z);
    dtmp_grad2 = fftw_malloc(N * sizeof *dtmp_grad2);
    dtmp_lap = fftw_malloc(N * sizeof *dtmp_lap);

    if (!(grid && field && field_new && dfield && dfield_new && field_buf &&
        time_buf && rho_buf && pow_spec && pow_spec_buf &&
        cfftw_tmp && cfftw_tmp_x && cfftw_tmp_y && cfftw_tmp_z &&
        dtmp_x && dtmp_y && dtmp_z && dtmp_grad2 && dtmp_lap))
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
	for (size_t i = 0; i < Nx; ++i)
	{
		grid[i] = ax + (bx - ax) * i / Nx;
	}
    for (size_t i = Nx; i < Nx+Ny; ++i)
    {
        grid[i] = ay + (by - ay) * (i-Nx) / Ny;
    }
    for (size_t i = Nx+Ny; i < Nx+Ny+Nz; ++i)
    {
        grid[i] = az + (bz - az) * (i-Nx-Ny) / Nz;
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
    double start =  get_wall_time();
    #endif
    p_fw_3d = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, cfftw_tmp,
                                            FFTW_DEFAULT_FLAG);
    p_bw_3d = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, cfftw_tmp, field,
                                            FFTW_DEFAULT_FLAG);
    #ifdef SHOW_TIMING_INFO
    fftw_time_plan += get_wall_time() - start;
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
    srand(1113);
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
                                                const double *phases) {
    // localized for higgs metastability potential
    double phi0 = 0.04;
    double lambda = 20.0;
    // return phi0 * 0.5 * (1.0 + cos(x)) * exp(-lambda * x * x);
    return phi0 * 0.125 * (1.0 + cos(x)) * (1.0 + cos(y)) * (1.0 + cos(z)) *
       exp(-lambda * (x * x + y * y + z * z));

    // some simple waves for notch or step potential simulations
    // double frac = 0.4;
    // double phi0 = 0.73 * frac;
    // double deltaphi = phi0 / frac;
	// return phi0 + deltaphi *
    //             (cos(1.0 * x + phases[0]) + cos(-1.0 * x + phases[1]) +
    //              cos(1.0 * y + phases[2]) + cos(-1.0 * y + phases[3]) +
    //              cos(1.0 * z + phases[4]) + cos(-1.0 * z + phases[5]));
}

// initial values of the time deriv. of the scalar field, make sure its periodic
double dphi_init(const double x, const double y, const double z) {
	return 0.0;
}

void free_and_destroy_all() {
    destroy_and_cleanup_fftw();
    free_external();
}

// destroy the fftw plans and call cleanup for internal fftw3 cleanup
void destroy_and_cleanup_fftw() {
    fftw_destroy_plan(p_fw_3d);
    fftw_destroy_plan(p_bw_3d);
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
    free(time_buf);
    free(f_a_buf);
    free(rho_buf);
    free(pow_spec);
    free(pow_spec_buf);
    fftw_free(cfftw_tmp);
    fftw_free(cfftw_tmp_x);
    fftw_free(cfftw_tmp_y);
    fftw_free(cfftw_tmp_z);
    fftw_free(dtmp_x);
    fftw_free(dtmp_y);
    fftw_free(dtmp_z);
    fftw_free(dtmp_grad2);
    fftw_free(dtmp_lap);
	RUNTIME_INFO(puts("Freed external variables.\n"));
}

// -------------------------printing functions----------------------------------
// for debugging mostly
void print_vector(const double *vector, const size_t N) {
    for (size_t i = 0; i < N; i++)
    {
        printf("%f\n", vector[i]);
    }
}
