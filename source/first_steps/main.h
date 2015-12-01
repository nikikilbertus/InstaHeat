#ifndef __MAIN__
#define __MAIN__

#include <stddef.h>
#include <complex.h>
#include <fftw3.h>

double get_wall_time();

/*
macros for debugging, testing and printing additional information during
execution
*/
#define SHOW_RUNTIME_INFO
// #define RUN_TESTS_ONLY
// #define DEBUG
// #define ENABLE_PROFILER

#ifdef SHOW_RUNTIME_INFO
#define RUNTIME_INFO(f) do {\
		(f); \
	} while (0)
#else
#define RUNTIME_INFO(f)
#endif

/*
how many threads to use for openmp parallelization
*/
#define THREAD_NUMBER			(4)

/*
check for NaNs during time evolution
*/
// #define CHECK_FOR_NAN

/*
apply a frequency cutoff filter during the time evolution
*/
// #define ENABLE_FFT_FILTER
// cutoff fraction used in spectral filtering
#define CUTOFF_FRACTION 		(1.0/3.0)

/*
the plan flag used for fftw plans
*/
#define FFTW_DEFAULT_FLAG 		FFTW_ESTIMATE

/*
should the power spectrum be written to disc
*/
#define POWER_SPECTRUM_MODE		(2) // default: 2
// how many bins for |k| are used in the computation of the power spectrum
#define POWER_SPECTRUM_BINS		(30)
// number of timeslices written to disc
#define POWER_SPECTRUM_NUMBER	(100)
// name of the file for the power spectrum
#define POWER_SPECTRUM_NAME		("pow_spec_000")

/*
file handling and write to disc parameters
*/
#define FIELD_MODE				(2) // default: 2
// number of timeslices written to disc
#define FIELD_NUMBER			(100)
// name of the file for the field
#define FIELD_NAME				("field_000")

// where to write the files
#define DATAPATH				("../../../data/")
// maximal length of file names
#define FILE_NAME_BUFFER_SIZE	(64)

/*
mathematical constants and macros
*/
#define PI (3.141592653589793238462643383279)
#define MAX(x, y)				((x) > (y) ? (x) : (y))
#define MIN(x, y)				((x) > (y) ? (y) : (x))

/*
simulation parameters
*/
// spatial
#define GRIDPOINTS_X  			(32)
#define GRIDPOINTS_Y  			(32)
#define GRIDPOINTS_Z  			(32)
#define GRIDPOINTS_TOTAL		((GRIDPOINTS_X)*(GRIDPOINTS_Y)*(GRIDPOINTS_Z))
#define SPATIAL_LOWER_BOUND_X	(-PI)
#define SPATIAL_UPPER_BOUND_X 	(PI)
#define SPATIAL_LOWER_BOUND_Y	(-PI)
#define SPATIAL_UPPER_BOUND_Y 	(PI)
#define SPATIAL_LOWER_BOUND_Z	(-PI)
#define SPATIAL_UPPER_BOUND_Z 	(PI)
// temporal
#define DELTA_T					(-0.1) // negative for manual adjustment
#define INITIAL_TIME 			(0.0)
#define FINAL_TIME	 			(250.0)
// potential
#define MASS 					(1.0)
#define COUPLING 				(1.0) // coupling in a phi4 potential
#define LAMBDA					(1.876e-4) // "cosmological constant"

// conversion form 3D indices to 1D index
// #define idx(i,j,k) 	((k) + GRIDPOINTS_Z * ( (j) + GRIDPOINTS_Y * (i)))
// #define idx_z(lin) 	( (lin) % GRIDPOINTS_Z )
// #define idx_y(lin) 	( (lin - (idx_z(lin)) / GRIDPOINTS_Z) % GRIDPOINTS_Y )
// #define idx_x(lin) 	( ( (lin - (idx_z(lin)) / GRIDPOINTS_Z) - (idx_y(lin)) ) / GRIDPOINTS_Y )

/*
representing one dimension of a multi dimensional grid
*/
typedef struct {
	size_t N;
	double a;
	double b;
}grid_dimension_t;

/*
encapsulate timing related parameters
*/
typedef struct {
	size_t Nt; // Number of timesteps
	double dt; // size of timestep delta t
	double ti; // initial time
	double tf; // final time
}timing_t;

/*
file handling parameters
mode: 0 - write nothing, 1 - write last only, 2 - write according to num/skip
3- write every timeslice
*/
typedef struct {
	uint8_t mode_field;
	size_t num_field;
	size_t skip_field;
	uint8_t mode_powspec;
	size_t num_powspec;
	size_t skip_powspec;
	size_t bins_powspec;
	size_t filename_buf;
	char *name_field;
	char *name_powspec;
	char *datapath;
}file_parameters_t;

/*
simulation parameters struct
*/
typedef struct {
	grid_dimension_t x;
	grid_dimension_t y;
	grid_dimension_t z;
	timing_t t;
	double cutoff_fraction; // used in spectral filtering during time evolution
	file_parameters_t file;
}parameters_t;

extern parameters_t pars;

// spatial gridpoints
extern double *grid;

// solutions for the scalar field we evolve and its temporal derivative
extern double *field, *dfield;
extern double *field_new, *dfield_new;

// solution for the FRW euqations
extern double *frw_a;
extern double f_a, df_a;
extern double f_a_new, df_a_new;

// T^{00} component of the field
extern double *rho;

// power spectrum
extern double *pow_spec;

// default arrays for real to complex dfts
extern complex *cfftw_tmp;
extern complex *cfftw_tmp_x;
extern complex *cfftw_tmp_y;
extern complex *cfftw_tmp_z;

// general purpose memory block for temporary use (eg for gradient)
extern double *dtmp_x;
extern double *dtmp_y;
extern double *dtmp_z;
extern double *dtmp_grad2;
extern double *dtmp_lap;

// fftw plans
extern fftw_plan p_fw_3d;
extern fftw_plan p_bw_3d;

// monitoring the time taken by certain parts
extern double fftw_time_exe;
extern double fftw_time_plan;

#endif