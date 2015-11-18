#ifndef __MAIN__
#define __MAIN__

#include <stddef.h>
#include <complex.h>
#include <fftw3.h>

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
apply a frequency cutoff based filter during the time evolution
*/
// #define ENABLE_FFT_FILTER
// #define ENABLE_ADAPTIVE_FILTER

/*
the plan flag used for fftw plans
*/
#define FFTW_DEFAULT_FLAG 		FFTW_MEASURE

/*
paths to files
*/
#define DATAPATH				"../../../data/"
#define FILE_NAME_BUFFER_SIZE	(64)

/*
mathematical constants
*/
#define PI (3.141592653589793238462643383279)

/*
simulation parameters
*/
// Number of gridpoints (needs to be odd for fourier grid)
#define GRIDPOINTS_X  			(16)
#define GRIDPOINTS_Y  			(16)
#define GRIDPOINTS_Z  			(16)
#define GRIDPOINTS_TOTAL		((GRIDPOINTS_X)*(GRIDPOINTS_Y)*(GRIDPOINTS_Z))
// timestep
#define DELTA_T					(-0.1)
// starttime (initial time)
#define INITIAL_TIME 			(0.0)
// final time
#define FINAL_TIME	 			(100.0)
// boundaries of spatial region
#define SPATIAL_LOWER_BOUND_X	(-PI)
#define SPATIAL_UPPER_BOUND_X 	(PI)
#define SPATIAL_LOWER_BOUND_Y	(-PI)
#define SPATIAL_UPPER_BOUND_Y 	(PI)
#define SPATIAL_LOWER_BOUND_Z	(-PI)
#define SPATIAL_UPPER_BOUND_Z 	(PI)
// what fraction of the evolution should be written to disc
#define WRITE_OUT_SIZE			(-1)
// if defined, only the last timeslice is written to file
#define WRITE_OUT_LAST_ONLY
// mass parameter
#define MASS 					(1.0)
// coupling in a phi4 potential
#define COUPLING 				(1.0)
// cutoff fraction used in spectral filtering during time evolution (disabled
// when ENABLE_ADAPTIVE_FILTER)
#define CUTOFF_FRACTION 		(0.5)
// in adaptive filtering, amount of energy required within cutoff (cumulative
// power spectral density)
#define CPSD_FRACTION 			(0.99)

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
simulation parameters struct
*/
typedef struct {
	grid_dimension_t x;
	grid_dimension_t y;
	grid_dimension_t z;
	size_t Nt; // Number of timesteps
	size_t nt; // counter for timesteps (not currently used)
	double t;  // current time (not currently used)
	double dt; // size of timestep delta t
	double ti; // initial time
	double tf; // final time
	double a;  // lower boundary of spatial domain
	double b;  // upper boudnary of spatial domain
	double cutoff_fraction; // used in spectral filtering during time evolution
	char *field_name;
	size_t file_row_skip;
}parameters_t;

extern parameters_t pars;

/*
spatial gridpoints
*/
extern double *grid;

/*
GENERAL REMARKS
- Dx/y/z as suffix indicates spatial derivatives
- d as prefix indicates a temporal derivative
- in the variable 'field' we save phi and dphi
*/

// solutions for the field we evolve and its temporal derivative
extern double *field;

// solution for the FRW euqations
extern double *frw_a;

// total energy of the field
extern double *rho;

// default arrays for real to complex dfts
extern complex *cfftw_tmp_x;
extern complex *cfftw_tmp_y;
extern complex *cfftw_tmp_z;
extern complex *cfftw_tmp_zx;
extern complex *cfftw_tmp_zy;
extern complex *cfftw_tmp_zz;

// general purpose memory block for temporary use (eg for gradient)
extern double *dtmp_x;
extern double *dtmp_y;
extern double *dtmp_z;
extern double *dtmp_grad2;
extern double *dtmp_lap;

// frequently used fftw plans
extern fftw_plan p_fw_laplacian;
extern fftw_plan p_bw_laplacian;
extern fftw_plan p_fw_Dx;
extern fftw_plan p_bw_Dx;
extern fftw_plan p_fw_Dy;
extern fftw_plan p_bw_Dy;
extern fftw_plan p_fw_Dz;
extern fftw_plan p_bw_Dz;

extern double fftw_time_exe;
extern double fftw_time_plan;

#endif