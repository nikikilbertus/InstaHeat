#ifndef __MAIN__
#define __MAIN__

#include <stddef.h>

/*
macros for debugging, testing and printing additional information during
execution
*/
// #define RUN_TESTS_ONLY
// #define PRINT_SPECTRAL_OPERATORS

/*
mathematical constants
*/
#define PI (3.141592653589793238462643383279)

/*
simulation parameters
*/
// Number of gridpoints (needs to be odd for fourier grid)
#define NOGP 		(31)
// maximal number of timesteps
#define MAXTS 		(10000)
// timestep
#define DT			(0.01)
// starttime (initial time)
#define TI 			(0.0)
// final time
#define TF 			(PI*100)
// boundaries of spatial region
#define LOW_BND		(-PI)
#define UP_BND 		(PI)

/*
counter for the timesteps and number of timesteps
*/
extern size_t nt;
extern size_t TS;

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
d as a prefix indicates temporal derivatives EXAMPLE: dphiD2 is
\dot{\sigma}^{\prime \prime} (1 temporal and 2 spatial derivatives)

We might not want to keep all derivatives through out the whole evolution in
memory later.
*/

//solutions for phi and its spatial derivatives
extern double *phi, *phiD1, *phiD2;

//solution for dphi and its spatial derivatives
extern double *dphi, *dphiD1, *dphiD2;

#endif