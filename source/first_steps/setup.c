#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "setup.h"
#include "main.h"

void allocate_and_initialize_all(parameters_t *pars) {
    initialize_threading();
	initialize_parameters(pars);
    allocate_external(pars);
    mk_grid(pars);
    mk_fftw_plans(pars);
    mk_initial_conditions(pars);
}

void initialize_threading() {
    int threadnum, threadinit;
    threadinit = fftw_init_threads();
    if (threadinit == 0)
    {
        fputs("Could not initialize fftw threads.", stderr);
        exit(EXIT_FAILURE);
    }
    threadnum = THREAD_NUMBER;
    omp_set_num_threads(threadnum);
    fftw_plan_with_nthreads(threadnum);
    RUNTIME_INFO(printf("Running omp & fftw with %d thread(s)\n\n", threadnum));
}

void initialize_parameters(parameters_t *pars) {
	pars->x.N = GRIDPOINTS_X;
    pars->y.N = GRIDPOINTS_Y;
    pars->z.N = GRIDPOINTS_Z;
	pars->x.a = SPATIAL_LOWER_BOUND_X;
	pars->x.b = SPATIAL_UPPER_BOUND_X;
    pars->y.a = SPATIAL_LOWER_BOUND_Y;
    pars->y.b = SPATIAL_UPPER_BOUND_Y;
    pars->z.a = SPATIAL_LOWER_BOUND_Z;
    pars->z.b = SPATIAL_UPPER_BOUND_Z;

    pars->Ntot = pars->x.N * pars->y.N * pars->z.N;

    pars->t.dt = DELTA_T > 0 ? DELTA_T : pars->t.dt;
    pars->t.t  = INITIAL_TIME;
    pars->t.ti = INITIAL_TIME;
    pars->t.tf = FINAL_TIME;
	pars->t.Nt = ceil((pars->t.tf - pars->t.ti) / pars->t.dt) + 1;

    pars->cutoff_fraction = CUTOFF_FRACTION;

    pars->file.index = 0;
    pars->file.buf_size = WRITE_OUT_BUFFER_NUMBER;
    pars->file.skip = TIME_STEP_SKIPS;

    pars->file.bins_powspec = POWER_SPECTRUM_BINS;

    RUNTIME_INFO(puts("Initialized parameters.\n"));
}

/*
allocate memory and initialize all external variables
*/
void allocate_external(parameters_t *pars) {
    size_t Nx   = pars->x.N;
    size_t Ny   = pars->y.N;
    size_t Nz   = pars->z.N;
    size_t Ntot = pars->Ntot;
    size_t buf_size = pars->file.buf_size;
    size_t bins = pars->file.bins_powspec;

    //grid points
    grid = malloc((Nx + Ny + Nz) * sizeof *grid);

    //solutions for the field and the temporal derivative
    field      = fftw_malloc(2 * Ntot * sizeof *field);
    field_new  = fftw_malloc(2 * Ntot * sizeof *field_new);
    dfield     = fftw_malloc(2 * Ntot * sizeof *dfield);
    dfield_new = fftw_malloc(2 * Ntot * sizeof *dfield_new);
    // write out buffer for the field phi
    field_buf = calloc(buf_size * 2 * Ntot, sizeof *field_buf);

    // write out buffer for a
    time_buf = calloc(buf_size, sizeof *time_buf);

    // write out buffer for a
    f_a_buf = calloc(buf_size, sizeof *f_a_buf);

    // write out buffer for a
    rho_buf = calloc(buf_size, sizeof *rho_buf);

    // power spectrum and write out buffer
    pow_spec = calloc(bins, sizeof *pow_spec);
    pow_spec_buf = calloc(buf_size * bins, sizeof *pow_spec_buf);

    // default arrays to save coefficients of real to complex transforms
    size_t ncz = Nz / 2 + 1;
    cfftw_tmp   = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp);
    cfftw_tmp_x = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_x);
    cfftw_tmp_y = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_y);
    cfftw_tmp_z = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_z);

    // general purpose double memory blocks for temporary use
    dtmp_x = fftw_malloc(2 * Ntot * sizeof *dtmp_x);
    dtmp_y = fftw_malloc(2 * Ntot * sizeof *dtmp_y);
    dtmp_z = fftw_malloc(2 * Ntot * sizeof *dtmp_z);
    dtmp_grad2 = fftw_malloc(Ntot * sizeof *dtmp_grad2);
    dtmp_lap = fftw_malloc(Ntot * sizeof *dtmp_lap);

    if (!(grid && field && field_new && dfield && dfield_new && field_buf &&
        time_buf && f_a_buf && rho_buf && pow_spec && pow_spec_buf &&
        cfftw_tmp && cfftw_tmp_x && cfftw_tmp_y && cfftw_tmp_z &&
        dtmp_x && dtmp_y && dtmp_z && dtmp_grad2 && dtmp_lap))
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }
    RUNTIME_INFO(puts("Allocated memory for external variables.\n"));
}

/*
	make the N fourier spectral gridpoints for the interval [a, b],
	and the spectral operators (N \times N matrices) for first and second order
	derivatives; x, D1 and D2 need to be initialized and memory needs to be
	allocated before calling this function
	see mathematica notebook spectral_operators for further info
*/
void mk_grid(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    double ax = pars->x.a;
    double bx = pars->x.b;
    double ay = pars->y.a;
    double by = pars->y.b;
    double az = pars->z.a;
    double bz = pars->z.b;

	if (Nx <= 0 || Ny <= 0 || Nz <= 0)
	{
		fputs("Need positive number of gridpoints", stderr);
	}

	// Grid points
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

    // Console output for debugging
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

/*
create all fftw plans, IMPORTANT: create BEFORE initializing arrays!
*/
void mk_fftw_plans(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;

    #ifdef SHOW_TIMING_INFO
    double start =  get_wall_time();
    #endif
    p_fw_3d = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, cfftw_tmp,
                                            FFTW_DEFAULT_FLAG);
    p_bw_3d = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, cfftw_tmp, field,
                                            FFTW_DEFAULT_FLAG);
    #ifdef SHOW_TIMING_INFO
    double end =  get_wall_time();
    fftw_time_plan += end - start;
    #endif
    RUNTIME_INFO(puts("Created fftw plans.\n"));
}

/*
setup initial conditions for the field
*/
void mk_initial_conditions(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;
    size_t osx, osy;
    double x, y, z;

    // random phases
    srand(1113);
    double *theta = calloc(6, sizeof *theta);
    for (size_t i = 0; i < 6; ++i)
    {
        theta[i] = 2.0 * PI * (double)rand() / (double)RAND_MAX;
    }

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
                field[Ntot + osy + k] = dphi_init(x, y, z);
            }
        }
    }

    f_a = 1.0;

    // Console output for debugging
    #ifdef DEBUG
    	puts("phi");
        print_vector(field, Ntot);
        puts("\ndphi");
        print_vector(field + Ntot, Ntot);
        puts("\n");
    #endif

    free(theta);
    RUNTIME_INFO(puts("Initialized the field and its temporal derivative.\n"));
}

/*
example functions for initial conditions
those need to be periodic in the spatial domain
*/
double phi_init(double x, double y, double z, double *phases) {
    double frac = 0.4;
    double phi0 = 0.73 * frac;
    double deltaphi = phi0 / frac;
	return phi0 + deltaphi *
                (cos(1.0 * x + phases[0]) + cos(-1.0 * x + phases[1]) +
                 cos(1.0 * y + phases[2]) + cos(-1.0 * y + phases[3]) +
                 cos(1.0 * z + phases[4]) + cos(-1.0 * z + phases[5]));
	// return tanh(pow(x, 8));
}

double dphi_init(double x, double y, double z) {
	// return -sin(x) * sin(y) * sin(z);
	return 0.0;
}

void free_and_destroy_all(parameters_t *pars) {
    destroy_fftw_plans();
    free_all_external(pars);
}

/*
destroy all fftw plans
*/
void destroy_fftw_plans() {
    fftw_destroy_plan(p_fw_3d);
    fftw_destroy_plan(p_bw_3d);
    RUNTIME_INFO(puts("Destroyed fftw plans.\n"));
}

/*
free all allocated memory
*/
void free_all_external(parameters_t *pars) {
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

//****************************** printing functions
/*
for debugging mostly
*/
void print_matrix(double *matrix, size_t N){
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            printf("%f\t", matrix[j*N+i]);
        }
        printf("\n");
    }
}

void print_vector(double *vector, size_t N) {
    for (size_t i = 0; i < N; i++)
    {
        printf("%f\n", vector[i]);
    }
}