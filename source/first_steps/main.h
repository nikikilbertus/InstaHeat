#ifndef __MAIN__
#define __MAIN__

#include <stddef.h>

/*
macros for debugging, testing and printing additional information during
execution
*/
#define DEBUG_MODE
// #define RUN_TESTS_ONLY
// #define PRINT_SPECTRAL_OPERATORS
// #define PRINT_INITIAL_CONDITIONS

#ifdef DEBUG_MODE
#define DEBUG(f) do {\
		(f); \
	} while (0)
#else
#define DEBUG(f)
#endif

/*
choose the differentiation procedure by defining either one of the two
*/
#define SPECTRAL_OPERATOR_DERIVATIVE
// #define FFT_DERIVATIVE

/*
check for NaNs during time evolution
*/
// #define CHECK_FOR_NAN

/*
mathematical constants
*/
#define PI (3.141592653589793238462643383279)

/*
simulation parameters
*/
// Number of gridpoints (needs to be odd for fourier grid)
#define GRIDPOINTS_SPATIAL 		(31)
// timestep
#define DELTA_T					(0.1)
// starttime (initial time)
#define INITIAL_TIME 			(0.0)
// final time
#define FINAL_TIME	 			(10.0)
// boundaries of spatial region
#define SPATIAL_LOWER_BOUND		(-PI)
#define SPATIAL_UPPER_BOUND 	(PI)

/*
simulation parameters struct
*/
typedef struct {
	size_t Nx; // Number of gridpoints in spatial direction (odd for fourier)
	size_t Nt; // Number of timesteps
	size_t nt; // counter for timesteps
	double dt; // size of timestep delta t
	double ti; // initial time
	double tf; // final time
	double a;
	double b;
}parameters_t;

extern parameters_t pars;

/*
gridpoints: x; spectral operators (matrices for first and second order
derivatives: D1, D2)
*/
extern double *x, *D1, *D2;

/*
variable for lapack routines (everything is passed by reference!)
*/
// currently empty

/*
GENERAL REMARKS
D1,D2 as suffix indicates 1st and 2nd order spatial derivatives

in the variable field we save phi and dphi
*/

//solutions for the field we evolve and its spatial derivatives
extern double *field, *fieldD1, *fieldD2;

#endif