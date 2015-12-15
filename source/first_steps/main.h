#ifndef __MAIN__
#define __MAIN__

#include <complex.h>
#include <fftw3.h>

// -------------------mathematical constants and macros-------------------------
#define PI                      (3.141592653589793238462643383279)
#define MAX(x, y)               ((x) > (y) ? (x) : (y))
#define MIN(x, y)               ((x) > (y) ? (y) : (x))
// seed for random number generator in creation of initial conditions
#define SEED                    (1113)

/*
compiler switches for debugging, testing, profiling and additional information
during execution
*/
#define SHOW_RUNTIME_INFO // recommended
#define SHOW_TIMING_INFO // recommended
// #define CHECK_FOR_NAN // not recommended (performance)
// #define ENABLE_PROFILER // only recommended for debugging
// #define DEBUG // only recommended for debugging (huge output!)
// #define RUN_TESTS_ONLY // testing only (tests.c)

// always leave this here uncommented
#ifdef SHOW_RUNTIME_INFO
#define RUNTIME_INFO(f) do {\
        (f); \
    } while (0)
#else
#define RUNTIME_INFO(f)
#endif

// ------------file handling parameters for writing to disk---------------------
// the output is bundled in one .h5 file, enter path here
#define DATAPATH                ("../../../data/run.h5")
/**
 *  write out data is buffered to not access hard drive too frequently (actually
 *  there are several levels of buffering, since hdf5 does it's own buffering 
 *  too), this parameter determines how many time slices are buffered before
 *  writing them to disk, beware of the memory consumption of large buffers!
 */
#define WRITE_OUT_BUFFER_NUMBER (20)
// how many timeslices to skip in between writing to file (1: write out all)
#define TIME_STEP_SKIPS         (30)
/**
 *  there is a (very crude and biased!) estimation of the power spectrum to
 *  track stability, therefore we sum up fourier coefficients into bins
 *  depending on the norm of their k vector, this gives the number of bins used
 */
#define POWER_SPECTRUM_BINS     (50)

// -----------------------simulation preferences--------------------------------
// how many threads to use for openmp parallelization (also used by fftw)
// if <= 0, the return value of omp_get_max_threads() is used
#define THREAD_NUMBER           (0)
// the plan flag used for fftw plans
#define FFTW_DEFAULT_FLAG       (FFTW_ESTIMATE)
/** 
 *  apply a frequency cutoff filter at each time step during the time evolution
 *  (compiler switch) the specific cutoff (e.g. step function with certain
 *  fraction (2/3 rule), or fourier smoothing) is determined by the function
 *  filter_window_function in evolution_toolkit.c, this should not be necessary
 *  for smaller grids, and simple (not too nonlinear) scenarios rather try to
 *  adjust tolerances and see what happens to the power spectrum
*/
#define ENABLE_FFT_FILTER

// ------------------computational domain---------------------------------------
// spatial
// TODO currently one needs at least 2, better 4 points in each direction,
// because only the 3d fftw is used. 2d and 1d will be properly implemented
// soonish
#define GRIDPOINTS_X            (64)
#define GRIDPOINTS_Y            (64)
#define GRIDPOINTS_Z            (64)
#define SPATIAL_LOWER_BOUND_X   (-PI)
#define SPATIAL_UPPER_BOUND_X   (PI)
#define SPATIAL_LOWER_BOUND_Y   (-PI)
#define SPATIAL_UPPER_BOUND_Y   (PI)
#define SPATIAL_LOWER_BOUND_Z   (-PI)
#define SPATIAL_UPPER_BOUND_Z   (PI)
// temporal
// initial step size for adaptive stepping (dopri853) or fixed step size (RK4) 
#define DELTA_T                 (0.001)
#define INITIAL_TIME            (0.0)
#define FINAL_TIME              (0.5)
#define MAX_STEPS               (1e6)
#define MINIMAL_DELTA_T         (1.0e-5)

// ----------------parameters used in the potential-----------------------------
#define MASS                    (1.0)
#define COUPLING                (1.0)      // coupling in a phi4 potential
#define LAMBDA                  (7.8e-2) // "cosmological constant"

// -------------------additional parameters for dopri853------------------------
// maximal/minimal rescaling of dt per step (don't change)
#define SMALLEST_SCALING        (0.333)
#define LARGEST_SCALING         (6.0)
// internal parameters for determining the rescaling of dt (don't change)
#define BETA                    (0.0) // ALPHA = 1.0/8.0 - BETA * 0.2
#define SAFE                    (0.9)
// error tolerancees, those can be changed (typical: between 1e-10 and 1e-3)
#define RELATIVE_TOLERANCE      (1.0e-9)
#define ABSOLUTE_TOLERANCE      (1.0e-9)

// ------------------------typedefs---------------------------------------------
// representing one spatial dimension of a multi dimensional grid
typedef struct {
    size_t N; // number of gridpoints
    double a; // lower bound of interval
    double b; // upper bound of interval
    double L; // lenght of the interval
}grid_dimension_t;

// encapsulate timing related parameters
typedef struct {
    size_t Nt; // number of timesteps (only relevant for fixed step size)
    double dt; // size of (initial) timestep delta t
    double ti; // initial time
    double tf; // final time
    double t;  // current time
}timing_t;

//file handling parameters
typedef struct {
    size_t id;           // h5 file id of the output file
    size_t dset_phi;     // h5 data set id of the scalar field
    size_t dset_powspec; // h5 data set id of the power spectrum
    size_t dset_time;    // h5 data set id of the time
    size_t dset_a;       // h5 data set id of the scaling parameter a
    size_t dset_rho;     // h5 data set id of the energy density rho
    size_t index;        // current index within the buffers
    size_t buf_size;     // size of the buffer
    size_t skip;         // how many timesteps to skip in between write out
    size_t bins_powspec; // how many bins are used for the power spectrum
}file_parameters_t;

/**
 *  simulation parameters
 *  throughout all files holds:
 *  Nx = number of grid points in the x direction
 *  Ny = number of grid points in the y direction
 *  Nz = number of grid points in the z direction
 *  N  = number of gridpoints for the whole spatial grid = Nx * Ny * Nz
 *  N2 = 2 * N
 *  Ntot = number of scalar equations = 2 * N + 1 (N for scalar field, N for its
 *  first temporal derivative, 1 for the FRW scaling parameter a)
 */
typedef struct {
    grid_dimension_t x;
    grid_dimension_t y;
    grid_dimension_t z;
    size_t N;
    timing_t t;
    file_parameters_t file;
}parameters_t;

// --------------------------global variables-----------------------------------
// we are using rather many global variables; that has the advantage of central
// one time allocation/initialization and deallocation;
// it also saves a lot of typing

// simulation parameters
extern parameters_t pars;

// spatial gridpoints
extern double *grid;

// time slices buffer
extern double *time_buf;

// the scalar field plus scaling parameter a we evolve, the temporal derivatives
// and the buffer for the scalar field 
extern double *field, *dfield;
extern double *field_new, *dfield_new;
extern double *field_buf;

// the scalar perturbation psi
extern double *psi;

// buffer for scaling parameter a
extern double *f_a_buf;

// energy density rho  = T^{00}_{\phi} and the buffer
extern double rho_avg;
extern double *rho;
extern double *rho_buf;

// power spectrum and the buffer
extern double *pow_spec;
extern double *pow_spec_buf;

// default arrays for real to complex dfts
extern complex *cfftw_tmp;
extern complex *cfftw_tmp_x;
extern complex *cfftw_tmp_y;
extern complex *cfftw_tmp_z;

// general purpose memory blocks for temporary use (e.g. for gradient)
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
extern double filter_time;
extern double h5_time_write;

// for timing information during execution
#ifdef SHOW_TIMING_INFO
double get_wall_time();
#endif

#endif
