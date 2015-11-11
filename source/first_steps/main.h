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
#define GRIDPOINTS_SPATIAL 		(51)
// timestep
#define DELTA_T					(-0.1)
// starttime (initial time)
#define INITIAL_TIME 			(0.0)
// final time
#define FINAL_TIME	 			(10)
// boundaries of spatial region
#define SPATIAL_LOWER_BOUND		(-PI)
#define SPATIAL_UPPER_BOUND 	(PI)
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

/*
simulation parameters struct
*/
typedef struct {
	size_t Nx; // Number of gridpoints in spatial direction (odd for fourier)
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
extern double *x;

/*
GENERAL REMARKS
- D1,D2 as suffix indicates 1st and 2nd order spatial derivatives
- d as prefix indicates a temporal derivative
- in the variable 'field' we save phi and dphi
*/

// solutions for the field we evolve its temporal derivative
extern double *field;

// solution for the FRW euqations
extern double *frw_a;

// total energy of the field
extern double *rho;

// default array and plans for real to complex dfts
extern complex *cfftw_tmp;

// general purpose 2*N double memory block for temporary use
extern double *dmisc_tmp;

#endif