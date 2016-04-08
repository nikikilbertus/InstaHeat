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
#define SEED                    (_SEED_)

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
#define RK4                 (0)
#define DOPRI853            (1)
#define INTEGRATION_METHOD  (_IM_)

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

#define VERSION_CONTROL         (_VC_)

//-------------which files to write to disk-------------------------------------
#define _OPHI_
#define _ODPHI_
#define _OPSI_
#define _ODPSI_
#define _ORHO_

#define _OPOWSPEC_

#define _OPHIM_
#define _ODPHIM_
#define _OPSIM_
#define _ODPSIM_
#define _ORHOM_

#define _OPHIV_
#define _ODPHIV_
#define _OPSIV_
#define _ODPSIV_
#define _ORHOV_

// ------------file handling parameters for writing to disk---------------------
// where to get the initial conditions from
#define IC_FROM_INTERNAL_FUNCTION   (0)
#define IC_FROM_H5_FILE             (1)
#define IC_FROM_DAT_FILE            (2)
#define IC_FROM_BUNCH_DAVIES        (3)
#define INITIAL_CONDITIONS          _IC_

// the output is bundled in one .h5 file, enter path here
#define DATAPATH                ("_PATH_")
/* #define DATAPATH                ("datanewlong.h5") */
#define INITIAL_DATAPATH        ("_IPATH_")


/**
 *  write out data is buffered to not access hard drive too frequently (actually
 *  there are several levels of buffering, since hdf5 does it's own buffering
 *  too), this parameter determines how many time slices are buffered before
 *  writing them to disk, beware of the memory consumption of large buffers!
 */
#define WRITE_OUT_BUFFER_NUMBER (_BUF_)

/**
 *  there is a (very crude and biased!) estimation of the power spectrum to
 *  track stability, therefore we sum up fourier coefficients into bins
 *  depending on the norm of their k vector, this gives the number of bins used
 */
#define POWER_SPECTRUM_BINS     (_BINS_)

// how many timeslices to skip in between writing to file (1: write out all)
#define TIME_STEP_SKIPS         (_TSKIP_)

// spatial output strides
#define STRIDE_X                (_XSKIP_)
#define STRIDE_Y                (_YSKIP_)
#define STRIDE_Z                (_ZSKIP_)

// -----------------------simulation preferences--------------------------------
// how many threads to use for openmp parallelization (also used by fftw)
// if <= 0, the return value of omp_get_max_threads() is used
#define THREAD_NUMBER           (_THREADS_)

// the plan flag used for fftw plans (ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE)
#define FFTW_DEFAULT_FLAG       (_FFTW_)

/**
 *  apply a frequency cutoff filter at each time step during the time evolution
 *  (compiler switch) the specific cutoff (e.g. step function with certain
 *  fraction (2/3 rule), or fourier smoothing) is determined by the function
 *  filter_window_function in evolution_toolkit.c, this should not be necessary
 *  for smaller grids, and simple (not too nonlinear) scenarios rather try to
 *  adjust tolerances and see what happens to the power spectrum
 */
#define _FILTER_

// choose between the parabolic or the hyperbolic equation to evolve psi
#define PSI_ELLIPTIC            (0)
#define PSI_HYPERBOLIC          (1)
#define PSI_PARABOLIC           (2)
#define PSI_METHOD              (_PM_)

// ------------------computational domain---------------------------------------
// spatial (order is important! use y=z=1 for 1D; use z=1 for 2D)
#define GRIDPOINTS_X            (_GPX_)
#define GRIDPOINTS_Y            (_GPY_)
#define GRIDPOINTS_Z            (_GPZ_)
#define SPATIAL_LOWER_BOUND_X   (_LX_)
#define SPATIAL_UPPER_BOUND_X   (_UX_)
#define SPATIAL_LOWER_BOUND_Y   (_LY_)
#define SPATIAL_UPPER_BOUND_Y   (_UY_)
#define SPATIAL_LOWER_BOUND_Z   (_LZ_)
#define SPATIAL_UPPER_BOUND_Z   (_UZ_)

// temporal
// initial step size for adaptive stepping (dopri853) or fixed step size (RK4)
#define DELTA_T                 (_DELT_)
#define INITIAL_TIME            (_TI_)
#define FINAL_TIME              (_TF_)
#define MAX_STEPS               (_MAXSTEPS_)
#define MINIMAL_DELTA_T         (_MINDELT_)

// ----------------parameters used in the potential-----------------------------
/* #define MASS                    (0.11026) // for 50 e-fold hom inflation */
#define MASS                    (_M_) // data_64
/* #define MASS                    (0.002000003836216) // compare_2 compare_psi */
/* #define MASS                    (16.418149637955437) // defrost */
/* #define MASS                    (9.973557010035818e-7) // pspectre defrost */
#define MASS_KARSTEN            (_MKARSTEN_)
#define MASS_PLANCK             (_MPLANCK_) // for bunch davies
#define COUPLING                (_COUPLING_) // coupling in a phi4 potential
#define LAMBDA                  (_LAMBDA_) // "cosmological constant"
/* #define A_INITIAL               (1.05249e3) // compare_2, compare_psi, 5500 */
/* #define A_INITIAL               (4.09376e3) // compare_2, compare_psi, 6000 */
/* #define A_INITIAL               (6.1625e2) // compare_2, compare_psi, 5450 */
/* #define A_INITIAL               (1.370074629050061e5) // data_64_0 */
/* #define A_INITIAL               (6.227758258677358e4) // data_64psi_1 */
/* #define A_INITIAL               (1.868327477603207e4) // data_64psi_2 */
/* #define A_INITIAL               (6.227611966276128e3) // data_64psi_3 */
#define A_INITIAL               (_A_)
// for notch potential test: LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5

// -------------------additional parameters for dopri853------------------------
// maximal/minimal rescaling of dt per step (don't change)
#define SMALLEST_SCALING        (_MINSCAL_)
#define LARGEST_SCALING         (_MAXSCAL_)

// internal parameters for determining the rescaling of dt (don't change)
#define BETA                    (_BETA_) // ALPHA = 1.0/8.0 - BETA * 0.2
#define SAFE                    (_SAFE_)

// error tolerancees, those can be changed (typical: between 1e-10 and 1e-3)
#define RELATIVE_TOLERANCE      (_RELTOL_)
#define ABSOLUTE_TOLERANCE      (_ABSTOL_)
// the timestep is limited from above by this fraction of the hubble time 1/H
#define MAX_DT_HUBBLE_FRACTION  (_HFRAC_)

// ------------------------struct definitions-----------------------------------
// representing one spatial dimension of a multi dimensional grid
struct grid_dimension
{
    size_t N; // number of gridpoints
    // depending on dimension, different upper bounds in for loops, see
    // initialization in setup.c for more information
    size_t M;
    double a; // lower bound of interval
    double b; // upper bound of interval
    double k; // a factor used in computing k vectors: k = 2 pi I / L
    double k2; // k2 = k*k = -4 pi^2 / L^2
    size_t stride; // strides for output
    size_t outN; // number of output points in this dimension
};

// encapsulate timing related parameters
struct timing
{
    size_t Nt; // number of timesteps (only relevant for fixed step size)
    double dt; // size of (initial) timestep delta t
    double ti; // initial time
    double tf; // final time
    double t;  // current time
};

// bundle data set identifiers for hdf5 output
struct datasets
{
    hsize_t field;
    hsize_t mean;
    hsize_t var;
    hsize_t dfield;
    hsize_t dmean;
    hsize_t dvar;
};

//file handling parameters
struct file_parameters
{
    hsize_t id;             // h5 file id of the output file
    struct datasets dset_phi;    // h5 data set ids for the scalar field phi
    struct datasets dset_psi;    // h5 data set ids for the perturbation psi
    struct datasets dset_rho;    // h5 data set ids for the energy density rho
    hsize_t dset_powspec;   // h5 data set id of the power spectrum
    hsize_t dset_time;      // h5 data set id of the time
    hsize_t dset_a;         // h5 data set id of the scaling parameter a
    size_t index;           // current index within the buffers
    size_t buf_size;        // size of the buffer
    size_t skip;            // how many timesteps to skip in between write out
    size_t bins_powspec;    // how many bins are used for the power spectrum
};

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
struct parameters
{
    struct grid_dimension x;
    struct grid_dimension y;
    struct grid_dimension z;
    size_t N;
    size_t Ntot;
    size_t Nall;
    size_t outN; // total number of spatial gridpoints for output (with strides)
    size_t M;
    size_t dim;
    struct timing t;
    struct file_parameters file;
};

// temporal veriables used for computations of gradients, ffts and the like
struct temporary
{
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
};

// grids for ksq and kx, ky, kz
struct k_grid
{
    double *sq;
    double *x;
    double *y;
    double *z;
};

// --------------------------global variables-----------------------------------
// we are using rather many global variables; that has the advantage of central
// one time allocation/initialization and deallocation;
// it also saves a lot of typing

// simulation parameters
extern struct parameters pars;

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

// filter mask for fourier filtering
extern double *filter;

// default arrays with temporary memory for real to complex dfts
extern struct k_grid kvec;

// default arrays with temporary memory for real to complex dfts
extern struct temporary tmp;

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
