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
	initialize_parameters(pars);
    allocate_external(pars);
    mk_grid(pars);
    mk_fftw_plans(pars);
    mk_initial_conditions(pars);
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

    pars->file.datapath = calloc(strlen(DATAPATH) + 1,
                                    sizeof *pars->file.datapath);
    pars->file.name_field = calloc(pars->file.filename_buf,
                                    sizeof *pars->file.name_field);
    pars->file.name_powspec = calloc(pars->file.filename_buf,
                                    sizeof *pars->file.name_field);;

    if (!pars->file.datapath || !pars->file.name_field || !pars->file.name_powspec)
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    strcpy(pars->file.datapath, DATAPATH);
    pars->file.filename_buf = FILE_NAME_BUFFER_SIZE;

    pars->file.mode_field = FIELD_MODE;
    switch (pars->file.mode_field)
    {
        case 0:
            pars->file.num_field = 0;
            pars->file.skip_field = pars->t.Nt + 1;
            RUNTIME_INFO(puts("Not writing field to disc."));
            break;
        case 1:
            pars->file.num_field = 0;
            pars->file.skip_field = pars->t.Nt - 1;
            RUNTIME_INFO(puts("Writing field only at last timeslice to disc."));
            break;
        case 2:
            pars->file.num_field = FIELD_NUMBER;
            pars->file.skip_field = pars->t.Nt / FIELD_NUMBER;
            RUNTIME_INFO(puts("WARNING: Only part of the evolution is written "
                    "to disc. Some timesteps in the end might be missing!"));
            break;
        case 3:
            pars->file.num_field = pars->t.Nt;
            pars->file.skip_field = 1;
            RUNTIME_INFO(puts("Writing field at every timeslice to disc."));
            break;
    }
    strcpy(pars->file.name_field, FIELD_NAME);

    pars->file.mode_powspec = POWER_SPECTRUM_MODE;
    switch (pars->file.mode_powspec)
    {
        case 0:
            pars->file.num_powspec = 0;
            pars->file.skip_powspec = pars->t.Nt + 1;
            RUNTIME_INFO(puts("Not writing power spectrum to disc."));
            break;
        case 1:
            pars->file.num_powspec = 0;
            pars->file.skip_powspec = pars->t.Nt - 1;
            RUNTIME_INFO(puts("Writing power spectrum only at last timeslice "
                                "to disc."));
            break;
        case 2:
            pars->file.num_powspec = POWER_SPECTRUM_NUMBER;
            pars->file.skip_powspec = pars->t.Nt / POWER_SPECTRUM_NUMBER;
            RUNTIME_INFO(puts("WARNING: Only part of the power spectrum is "
            "written to disc. Some timesteps in the end might be missing!"));
            break;
        case 3:
            pars->file.num_powspec = pars->t.Nt;
            pars->file.skip_powspec = 1;
            RUNTIME_INFO(puts("Writing power spectrum at every timeslice "
                                "to disc."));
            break;
    }
    pars->file.bins_powspec = POWER_SPECTRUM_BINS;
    strcpy(pars->file.name_powspec, POWER_SPECTRUM_NAME);

    RUNTIME_INFO(puts("Initialized parameters.\n"));
}

/*
allocate memory and initialize all external variables
*/
void allocate_external(parameters_t *pars) {
    size_t Nx   = pars->x.N;
    size_t Ny   = pars->y.N;
    size_t Nz   = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;
    size_t Nt   = pars->t.Nt;

    //grid points
    grid = malloc((Nx + Ny + Nz) * sizeof *grid);

    //solutions for the field and the temporal derivative (we are saving each
    //timestep: 2 * Nx * Nt space)
    field      = fftw_malloc(2 * Ntot * sizeof *field);
    field_new  = fftw_malloc(2 * Ntot * sizeof *field_new);
    dfield     = fftw_malloc(2 * Ntot * sizeof *dfield);
    dfield_new = fftw_malloc(2 * Ntot * sizeof *dfield_new);

    // solution for the scale parameter a: Nt space
    frw_a = calloc(Nt, sizeof *frw_a);

    // T00 of the scalar field over time: Nt space
    rho = calloc(Nt, sizeof *rho);

    // power spectrum
    pow_spec = calloc(pars->file.bins_powspec, sizeof *pow_spec);

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

    if (!(grid && field && frw_a && rho && pow_spec &&
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
create all fftw plans, IMPORTANT: create BEFORE initializsing arrays!
*/
void mk_fftw_plans(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;

    double start =  get_wall_time();
    p_fw_3d = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, cfftw_tmp,
                                            FFTW_DEFAULT_FLAG);
    p_bw_3d = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, cfftw_tmp, field,
                                            FFTW_DEFAULT_FLAG);
    double end =  get_wall_time();
    fftw_time_plan += end - start;
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

    frw_a[0] = 1.0;
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
    free(pars->file.datapath);
    free(pars->file.name_field);
    free(pars->file.name_powspec);
    free(grid);
    fftw_free(field);
    fftw_free(field_new);
    fftw_free(dfield);
    fftw_free(dfield_new);
    free(frw_a);
    free(rho);
    free(pow_spec);
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