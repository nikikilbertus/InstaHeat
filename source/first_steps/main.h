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
#define WRITE_OUT_POWER_SPECTRUM
// how many bins for |k| are used in the computation of the power spectrum
#define POWER_SPECTRUM_SHELLS	(70)
// number of timeslices written to disc ( < 0 --> write out all)
#define WRITE_OUT_SIZE_POW_SPEC	(100)

/*
file handling and write to disc parameters
*/
// number of timeslices written to disc ( < 0 --> write out all)
#define WRITE_OUT_SIZE			(100)
// ifdef: write only last timeslice to disc regardless of WRITE_OUT_SIZE
// #define WRITE_OUT_LAST_ONLY
#define DATAPATH				"../../../data/" // where to write the files
#define FILE_NAME_BUFFER_SIZE	(64) // maximal length of file names

/*
mathematical constants
*/
#define PI (3.141592653589793238462643383279)

/*
simulation parameters
*/
// spatial
#define GRIDPOINTS_X  			(16)
#define GRIDPOINTS_Y  			(16)
#define GRIDPOINTS_Z  			(16)
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
*/

/*
simulation parameters struct
*/
typedef struct {
	grid_dimension_t x;
	grid_dimension_t y;
	grid_dimension_t z;
	timing_t t;
	double cutoff_fraction; // used in spectral filtering during time evolution
	char *field_name;
	size_t file_write_size;
	size_t pow_spec_shells;
	size_t file_write_size_pow_spec;
}parameters_t;

extern parameters_t pars;

/*
spatial gridpoints
*/
extern double *grid;

// solutions for the scalar field we evolve and its temporal derivative
extern double *field;

// solution for the FRW euqations
extern double *frw_a;

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