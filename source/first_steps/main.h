#ifndef __MAIN__
#define __MAIN__

#include <stddef.h>
#include <complex.h>

/*
macros for debugging, testing and printing additional information during
execution
*/
#define SHOW_RUNTIME_INFO
// #define RUN_TESTS_ONLY
// #define DEBUG

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
mathematical constants
*/
#define PI (3.141592653589793238462643383279)

/*
simulation parameters
*/
// Number of gridpoints (needs to be odd for fourier grid)
#define GRIDPOINTS_X  			(11)
#define GRIDPOINTS_Y  			(11)
#define GRIDPOINTS_Z  			(11)
#define GRIDPOINTS_TOTAL		((GRIDPOINTS_X)*(GRIDPOINTS_Y)*(GRIDPOINTS_Z))
// timestep
#define DELTA_T					(-0.1)
// starttime (initial time)
#define INITIAL_TIME 			(0.0)
// final time
#define FINAL_TIME	 			(1.1)
// boundaries of spatial region
#define SPATIAL_LOWER_BOUND_X	(-PI)
#define SPATIAL_UPPER_BOUND_X 	(PI)
#define SPATIAL_LOWER_BOUND_Y	(-PI)
#define SPATIAL_UPPER_BOUND_Y 	(PI)
#define SPATIAL_LOWER_BOUND_Z	(-PI)
#define SPATIAL_UPPER_BOUND_Z 	(PI)
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
#define idx(i,j,k) 	((k) + GRIDPOINTS_Z * ( (j) + GRIDPOINTS_Y * (i)))
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

// solutions for the field we evolve its temporal derivative
extern double *field;

// solution for the FRW euqations
extern double *frw_a;

// total energy of the field
extern double *rho;

// default arrays for real to complex dfts
extern complex *cfftw_tmp_x;
extern complex *cfftw_tmp_y;
extern complex *cfftw_tmp_z;

// general purpose memory block for temporary use
extern double *dmisc_tmp_x;
extern double *dmisc_tmp_y;
extern double *dmisc_tmp_z;

#endif