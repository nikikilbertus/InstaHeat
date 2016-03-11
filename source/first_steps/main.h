#ifndef __MAIN__
#define __MAIN__

#include <complex.h>
#include <fftw3.h>
#include "hdf5.h"

// -------------------mathematical constants and macros-------------------------
#define PI                      (3.141592653589793238462643383279)
#define TWOPI                   (6.283185307179586476925286766559)
#define MAX(x, y)               (((x) > (y)) ? (x) : (y))
#define MIN(x, y)               (((x) > (y)) ? (y) : (x))

// seed for random number generator in creation of initial conditions
#define SEED                    (1139)

/*
compiler switches for debugging, testing, profiling and additional information
during execution
*/
#define SHOW_RUNTIME_INFO // recommended
#define SHOW_TIMING_INFO // recommended
/* #define CHECK_FOR_NAN // not recommended (performance) */
/* #define ENABLE_PROFILER // only recommended for debugging */
/* #define DEBUG // only recommended for debugging (huge output!) */
/* #define RUN_TESTS_ONLY // testing only (tests.c) */

// always leave this here uncommented
#ifdef SHOW_RUNTIME_INFO
#define INFO(f) do {\
        (f); \
    } while (0)
#else
#define INFO(f)
#endif

// ------------writing current commit hash to disk (hg, git, none)--------------
#define VERSION_CONTROL_NONE    (0)
#define VERSION_CONTROL_HG      (1)
#define VERSION_CONTROL_GIT     (2)

#define VERSION_CONTROL         VERSION_CONTROL_HG

//-------------which files to write to disk-------------------------------------
#define OUTPUT_PHI
/* #define OUTPUT_DPHI */
#define OUTPUT_PSI
/* #define OUTPUT_DPSI */
#define OUTPUT_RHO

#define OUTPUT_POWER_SPECTRUM

#define OUTPUT_PHI_MEAN
#define OUTPUT_DPHI_MEAN
#define OUTPUT_PSI_MEAN
#define OUTPUT_DPSI_MEAN
#define OUTPUT_RHO_MEAN

#define OUTPUT_PHI_VARIANCE
#define OUTPUT_DPHI_VARIANCE
#define OUTPUT_PSI_VARIANCE
#define OUTPUT_DPSI_VARIANCE
#define OUTPUT_RHO_VARIANCE

// ------------file handling parameters for writing to disk---------------------
// where to get the initial conditions from
#define IC_FROM_INTERNAL_FUNCTION   (0)
#define IC_FROM_H5_FILE             (1)
#define IC_FROM_DAT_FILE            (2)
#define IC_FROM_BUNCH_DAVIES        (3)
#define INITIAL_CONDITIONS          IC_FROM_DAT_FILE

// the output is bundled in one .h5 file, enter path here
#define DATAPATH                ("../../../data/valtest.h5")
/* #define DATAPATH                ("bla.h5") */
#define INITIAL_DATAPATH        ("../../../data/karsten/data_64.dat")


/**
 *  write out data is buffered to not access hard drive too frequently (actually
 *  there are several levels of buffering, since hdf5 does it's own buffering
 *  too), this parameter determines how many time slices are buffered before
 *  writing them to disk, beware of the memory consumption of large buffers!
 */
#define WRITE_OUT_BUFFER_NUMBER (5)

/**
 *  there is a (very crude and biased!) estimation of the power spectrum to
 *  track stability, therefore we sum up fourier coefficients into bins
 *  depending on the norm of their k vector, this gives the number of bins used
 */
#define POWER_SPECTRUM_BINS     (30)

// how many timeslices to skip in between writing to file (1: write out all)
#define TIME_STEP_SKIPS         (1)

// spatial output strides
#define STRIDE_X                (1)
#define STRIDE_Y                (1)
#define STRIDE_Z                (1)

// -----------------------simulation preferences--------------------------------
// how many threads to use for openmp parallelization (also used by fftw)
// if <= 0, the return value of omp_get_max_threads() is used
#define THREAD_NUMBER           (0)

// the plan flag used for fftw plans (ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE)
#define FFTW_DEFAULT_FLAG       (FFTW_PATIENT)

/**
 *  apply a frequency cutoff filter at each time step during the time evolution
 *  (compiler switch) the specific cutoff (e.g. step function with certain
 *  fraction (2/3 rule), or fourier smoothing) is determined by the function
 *  filter_window_function in evolution_toolkit.c, this should not be necessary
 *  for smaller grids, and simple (not too nonlinear) scenarios rather try to
 *  adjust tolerances and see what happens to the power spectrum
 */
#define ENABLE_FFT_FILTER

// choose between the parabolic or the hyperbolic equation to evolve psi
#define PSI_ELLIPTIC            (0)
#define PSI_HYPERBOLIC          (1)
#define PSI_PARABOLIC           (2)
#define PSI_METHOD              (PSI_HYPERBOLIC)

// ------------------computational domain---------------------------------------
// spatial (order is important! use y=z=1 for 1D; use z=1 for 2D)
#define GRIDPOINTS_X            (64)
#define GRIDPOINTS_Y            (64)
#define GRIDPOINTS_Z            (64)
#define SPATIAL_LOWER_BOUND_X   (0.0)
#define SPATIAL_UPPER_BOUND_X   (1.0)
#define SPATIAL_LOWER_BOUND_Y   (0.0)
#define SPATIAL_UPPER_BOUND_Y   (1.0)
#define SPATIAL_LOWER_BOUND_Z   (0.0)
#define SPATIAL_UPPER_BOUND_Z   (1.0)

// temporal
// initial step size for adaptive stepping (dopri853) or fixed step size (RK4)
#define DELTA_T                 (1.0e-5)
#define INITIAL_TIME            (0.0)
#define FINAL_TIME              (1.0e1)
#define MAX_STEPS               (1e15)
#define MINIMAL_DELTA_T         (1.0e-10)

// ----------------parameters used in the potential-----------------------------
/* #define MASS                    (0.11026) // for 50 e-fold hom inflation */
#define MASS                    (1.0e-2) // exactly karste
/* #define MASS                    (0.002000003836216) // compare_2 compare_psi */
/* #define MASS                    (6.0) // for compare.dat */
#define MASS_KARSTEN            (1.0e-2)
#define COUPLING                (1.0) // coupling in a phi4 potential
#define LAMBDA                  (4.721e-5) // "cosmological constant"
/* #define A_INITIAL               (1.05249e3) // compare_2, compare_psi, 5500 */
/* #define A_INITIAL               (4.09376e3) // compare_2, compare_psi, 6000 */
/* #define A_INITIAL               (6.1625e2) // compare_2, compare_psi, 5450 */
#define A_INITIAL               (1.370074629050061e+05) // data_64
/* #define A_INITIAL               (1.0) */
// for notch potential test: LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5

// -------------------additional parameters for dopri853------------------------
// maximal/minimal rescaling of dt per step (don't change)
#define SMALLEST_SCALING        (0.333)
#define LARGEST_SCALING         (6.0)

// internal parameters for determining the rescaling of dt (don't change)
#define BETA                    (0.0) // ALPHA = 1.0/8.0 - BETA * 0.2
#define SAFE                    (0.9)

// error tolerancees, those can be changed (typical: between 1e-10 and 1e-3)
#define RELATIVE_TOLERANCE      (1.0e-9)
#define ABSOLUTE_TOLERANCE      (1.0e-14)
// the timestep is limited from above by a fraction of the hubble time 1/H
#define MIN_DT_FRACTION         (1.0e-3)

// ------------------------typedefs---------------------------------------------
// representing one spatial dimension of a multi dimensional grid
typedef struct {
    size_t N; // number of gridpoints
    // depending on dimension, different upper bounds in for loops, see
    // initialization in setup.c for more information
    size_t M;
    double a; // lower bound of interval
    double b; // upper bound of interval
    complex k; // a factor used in computing k vectors: k = 2 pi I / L
    double k2; // k2 = k*k = -4 pi^2 / L^2
    size_t stride; // strides for output
    size_t outN; // number of output points in this dimension
}grid_dimension_t;

// encapsulate timing related parameters
typedef struct {
    size_t Nt; // number of timesteps (only relevant for fixed step size)
    double dt; // size of (initial) timestep delta t
    double ti; // initial time
    double tf; // final time
    double t;  // current time
}timing_t;

// bundle data set identifiers for hdf5 output
typedef struct {
    hsize_t field;
    hsize_t mean;
    hsize_t var;
    hsize_t dfield;
    hsize_t dmean;
    hsize_t dvar;
}datasets_t;

//file handling parameters
typedef struct {
    hsize_t id;             // h5 file id of the output file
    datasets_t dset_phi;    // h5 data set ids for the scalar field phi
    datasets_t dset_psi;    // h5 data set ids for the perturbation psi
    datasets_t dset_rho;    // h5 data set ids for the energy density rho
    hsize_t dset_powspec;   // h5 data set id of the power spectrum
    hsize_t dset_time;      // h5 data set id of the time
    hsize_t dset_a;         // h5 data set id of the scaling parameter a
    size_t index;           // current index within the buffers
    size_t buf_size;        // size of the buffer
    size_t skip;            // how many timesteps to skip in between write out
    size_t bins_powspec;    // how many bins are used for the power spectrum
}file_parameters_t;

/**
 *  simulation parameters
 *  throughout all files holds:
 *  Nx = number of grid points in the x direction
 *  Ny = number of grid points in the y direction
 *  Nz = number of grid points in the z direction
 *  N  = number of gridpoints for the whole spatial grid = Nx * Ny * Nz
 *  N2 = 2 * N, N3 = 3 * N, ...
 *  Ntot = number of scalar equations; depends on the used method
 *  elliptic: 2 * N + 1 (order: phi, dphi, a)
 *  parabolic: 3 * N + 1 (oder: phi, dphi, psi, a)
 *  hyperbolic: 4 * N + 1 (oder: phi, dphi, psi, dpsi, a)
 */
typedef struct {
    grid_dimension_t x;
    grid_dimension_t y;
    grid_dimension_t z;
    size_t N;
    size_t Ntot;
    size_t Nall;
    size_t outN; // total number of spatial gridpoints for output (with strides)
    size_t dim;
    timing_t t;
    file_parameters_t file;
}parameters_t;

// temporal veriables used for computations of gradients, ffts and the like
typedef struct {
    double  *xphi;
    double  *yphi;
    double  *zphi;
    complex *phic;
    complex *xphic;
    complex *yphic;
    complex *zphic;
    double  *grad;
    double  *lap;
    double  *deltarho;
    double  *f;
    complex *deltarhoc;
    complex *fc;
    complex *psic;
    complex *dpsic;
}temporary_t;

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

// contain phi, dphi, a, psi, dpsi (only first Ntot entries are integrated)
extern double *field;
extern double *dfield;
extern double *field_new;
extern double *dfield_new;

//TODO: create structure to bundle buffers (equal for phi and psi)
// buffers for the scalar inflaton field and the scalar metric perturbation
extern double *phi_buf;
extern double *dphi_buf;
extern double *psi_buf;
extern double *dpsi_buf;

extern double phi_mean;
extern double phi_var;
extern double *phi_mean_buf;
extern double *phi_var_buf;

extern double dphi_mean;
extern double dphi_var;
extern double *dphi_mean_buf;
extern double *dphi_var_buf;

extern double psi_mean;
extern double psi_var;
extern double *psi_mean_buf;
extern double *psi_var_buf;

extern double dpsi_mean;
extern double dpsi_var;
extern double *dpsi_mean_buf;
extern double *dpsi_var_buf;

// buffer for scaling parameter a
extern double *a_buf;

// energy density rho  = T^{00}_{\phi} and the buffer
extern double *rho;
extern double *rho_buf;
extern double rho_mean;
extern double rho_var;
extern double *rho_mean_buf;
extern double *rho_var_buf;

// pressure (only needed for PSI_METHOD == PSI_HYPERBOLIC, but always defined)
extern double *pressure;
extern double pressure_mean;

// power spectrum and the buffer
extern double *pow_spec;
extern double *pow_spec_buf;

// default arrays with temporary memory for real to complex dfts
extern temporary_t tmp;

// fftw plans
extern fftw_plan p_fw;
extern fftw_plan p_bw;

// monitoring the time taken by certain parts
extern double fftw_time_exe;
extern double fftw_time_plan;
extern double filter_time;
extern double poisson_time;
extern double h5_time_write;

// for timing information during execution
#ifdef SHOW_TIMING_INFO
double get_wall_time();
#endif

#endif
