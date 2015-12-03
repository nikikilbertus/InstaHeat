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
#define FFTW_DEFAULT_FLAG 		(FFTW_ESTIMATE)

/*
should the power spectrum be written to disk
*/
// how many bins for |k| are used in the computation of the power spectrum
#define POWER_SPECTRUM_BINS		(30)

/*
file handling and write to disk parameters
*/
// file name
#define DATAPATH				("../../../data/run.h5")
// how many timeslices to keep in memory before write out
#define WRITE_OUT_BUFFER_NUMBER	(20)
// how many timeslices to skip in between (1 to write each)
#define TIME_SLICE_SKIP_NUMBER  (1)

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
#define DELTA_T					(0.1) // negative for manual adjustment
#define INITIAL_TIME 			(0.0)
#define FINAL_TIME	 			(10.0)
// potential
#define MASS 					(1.0)
#define COUPLING 				(1.0)      // coupling in a phi4 potential
#define LAMBDA					(1.876e-4) // "cosmological constant"

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
	double dt; // size of (initial) timestep delta t
	double ti; // initial time
	double tf; // final time
	double t;  // current time
}timing_t;

/*
file handling parameters
*/
typedef struct {
	size_t id;			// h5 file id of the output file
	size_t dset_phi;	// h5 data set id of the field phi
	size_t dset_time;	// h5 data set id of the time
	size_t dset_a;		// h5 data set id of the scaling parameter a
	size_t dset_rho;	// h5 data set id of the energy density rho
	size_t index;		// current index within the buffers
	size_t buf_size;	// size of the buffer before dump to disk
	size_t skip;		// how many timesteps to skip in between write out
	size_t bins_powspec; // how many bins are used for the power spectrum
}file_parameters_t;

/*
simulation parameters struct
*/
typedef struct {
	grid_dimension_t x;
	grid_dimension_t y;
	grid_dimension_t z;
	size_t Ntot;
	timing_t t;
	double cutoff_fraction; // used in spectral filtering during time evolution
	file_parameters_t file;
}parameters_t;

extern parameters_t pars;

// spatial gridpoints
extern double *grid;

// time slices buffer
extern double *time_buf;

// solutions for the scalar field we evolve and its temporal derivative
extern double *field, *dfield;
extern double *field_new, *dfield_new;
extern double *field_buf;

// solution for the FRW euqations
extern double f_a, df_a;
extern double f_a_new, df_a_new;
extern double *f_a_buf;

// T^{00} component of the field
extern double rho;
extern double *rho_buf;

// power spectrum
extern double *pow_spec;
extern double *pow_spec_buf;

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