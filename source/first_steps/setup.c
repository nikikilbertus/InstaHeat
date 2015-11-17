#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
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
	pars->dt  = DELTA_T > 0 ? DELTA_T : pars->dt;
	pars->ti  = INITIAL_TIME;
	pars->tf  = FINAL_TIME;
	pars->x.a = SPATIAL_LOWER_BOUND_X;
	pars->x.b = SPATIAL_UPPER_BOUND_X;
    pars->y.a = SPATIAL_LOWER_BOUND_Y;
    pars->y.b = SPATIAL_UPPER_BOUND_Y;
    pars->z.a = SPATIAL_LOWER_BOUND_Z;
    pars->z.b = SPATIAL_UPPER_BOUND_Z;
	pars->Nt = ceil((pars->tf - pars->ti) / pars->dt);
	pars->nt = 0;
	pars->t  = 0.0;
    pars->cutoff_fraction = CUTOFF_FRACTION;
    #ifndef WRITE_OUT_LAST_ONLY
    if (WRITE_OUT_SIZE < 0)
    {
        RUNTIME_INFO(puts("WARNING: Only part of the evolution is written to "
                        "disc. Some timesteps in the end might be missing!"));
    }
    pars->file_row_skip = WRITE_OUT_SIZE < 0 ? 1 : (pars->Nt / WRITE_OUT_SIZE);
    #endif
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
    size_t Nt   = pars->Nt;

    // some space for the file names
    pars->field_name = calloc(FILE_NAME_BUFFER_SIZE, sizeof *pars->field_name);

    //grid points
    grid = malloc((Nx + Ny + Nz) * sizeof *grid);
    if (!grid)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    //solutions for the field and the temporal derivative (we are saving each
    //timestep: 2 * Nx * Nt space)
    field = fftw_malloc(2 * Ntot * sizeof *field);
    if (!field)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    // solution for the scale parameter a: Nt space
    frw_a = calloc(Nt, sizeof *frw_a);
    if (!frw_a)
    {
    	fputs("Allocating memory failed.", stderr);
    	exit(EXIT_FAILURE);
    }

    // T00 of the scalar field over time: Nt space
    rho = calloc(Nt, sizeof *rho);
    if (!rho)
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    // default arrays to save coefficients of real to complex transforms
    size_t ncx = Nx / 2 + 1;
    size_t ncy = Ny / 2 + 1;
    size_t ncz = Nz / 2 + 1;

    cfftw_tmp_x = fftw_malloc(ncx * Ny * Nz * sizeof *cfftw_tmp_x);
    cfftw_tmp_y = fftw_malloc(ncy * Nx * Nz * sizeof *cfftw_tmp_y);
    cfftw_tmp_z = fftw_malloc(ncz * Nx * Ny * sizeof *cfftw_tmp_z);
    if (!cfftw_tmp_x || !cfftw_tmp_y || !cfftw_tmp_z)
    {
        fputs("Allocating memory failed.", stderr);
        exit(EXIT_FAILURE);
    }

    // general purpose double memory blocks for temporary use
    dtmp_x = fftw_malloc(2 * Ntot * sizeof *dtmp_x);
    dtmp_y = fftw_malloc(2 * Ntot * sizeof *dtmp_y);
    dtmp_z = fftw_malloc(2 * Ntot * sizeof *dtmp_z);
    if (!dtmp_x || !dtmp_y || !dtmp_z)
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

    clock_t start = clock();
    p_fw_laplacian = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, field, cfftw_tmp_z,
                                            FFTW_DEFAULT_FLAG);
    p_bw_laplacian = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, cfftw_tmp_z, field,
                                            FFTW_DEFAULT_FLAG);
    int rank = 1;
    int nx[] = {Nx};
    int howmany = Ny * Nz;
    int *inembed = nx;
    int istride = Ny * Nz;
    int idist = 1;
    int *onembed = nx;
    int ostride = istride;
    int odist = 1;
    p_fw_Dx = fftw_plan_many_dft_r2c(rank, nx, howmany, field, inembed, istride,
                idist, cfftw_tmp_x, onembed, ostride, odist, FFTW_DEFAULT_FLAG);
    p_bw_Dx = fftw_plan_many_dft_c2r(rank, nx, howmany, cfftw_tmp_x, inembed,
            istride, idist, field, onembed, ostride, odist, FFTW_DEFAULT_FLAG);

    int ny[] = {Ny};
    howmany = Nz;
    inembed = ny;
    istride = Nz;
    idist   = 1;
    onembed = ny;
    ostride = istride;
    odist   = 1;
    p_fw_Dy = fftw_plan_many_dft_r2c(rank, ny, howmany, field, inembed, istride,
                idist, cfftw_tmp_y, onembed, ostride, odist, FFTW_DEFAULT_FLAG);
    p_bw_Dy = fftw_plan_many_dft_c2r(rank, ny, howmany, cfftw_tmp_y, inembed,
            istride, idist, field, onembed, ostride, odist, FFTW_DEFAULT_FLAG);

    int nz[] = {Nz};
    howmany = Nx * Ny;
    inembed = nz;
    istride = 1;
    idist = Nz;
    onembed = nz;
    ostride = istride;
    odist = Nz / 2 + 1;
    p_fw_Dz = fftw_plan_many_dft_r2c(rank, nz, howmany, field, inembed, istride,
                idist, cfftw_tmp_z, onembed, ostride, odist, FFTW_DEFAULT_FLAG);
    idist   = odist;
    odist   = Nz;
    p_bw_Dz = fftw_plan_many_dft_c2r(rank, nz, howmany, cfftw_tmp_z, inembed,
            istride, idist, field, onembed, ostride, odist, FFTW_DEFAULT_FLAG);

    clock_t end = clock();
    fftw_time_plan += (double)(end - start) / CLOCKS_PER_SEC;
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
                field[osy + k] = phi_init(x, y, z);
                field[Ntot + osy + k] = dphi_init(x, y, z);
            }
        }
    }

    frw_a[0] = 1.0;

    // Console output for debugging
    #ifdef DEBUG
    	puts("phi");
        print_vector(field, Ntot);
        puts("\ndphi");
        print_vector(field + Ntot, Ntot);
        puts("\n");
    #endif
    RUNTIME_INFO(puts("Initialized the field and its temporal derivative.\n"));
}

/*
example functions for initial conditions
those need to be periodic in the spatial domain
*/
double phi_init(double x, double y, double z) {
	return sin(x) * sin(y) * sin(z);
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
    fftw_destroy_plan(p_fw_laplacian);
    fftw_destroy_plan(p_bw_laplacian);
    fftw_destroy_plan(p_fw_Dx);
    fftw_destroy_plan(p_bw_Dx);
    fftw_destroy_plan(p_fw_Dy);
    fftw_destroy_plan(p_bw_Dy);
    fftw_destroy_plan(p_fw_Dz);
    fftw_destroy_plan(p_bw_Dz);
    RUNTIME_INFO(puts("Destroyed fftw plans.\n"));
}

/*
free all allocated memory
*/
void free_all_external(parameters_t *pars) {
    free(pars->field_name);
    free(grid);
    fftw_free(field);
    free(frw_a);
    free(rho);
    fftw_free(cfftw_tmp_x);
    fftw_free(cfftw_tmp_y);
    fftw_free(cfftw_tmp_z);
    fftw_free(dtmp_x);
    fftw_free(dtmp_y);
    fftw_free(dtmp_z);
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