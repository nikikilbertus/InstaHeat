#ifndef __MAIN__
#define __MAIN__

#include <complex.h>
#include <fftw3.h>
#include "hdf5.h"

/**
 * @file main.h
 * @brief Contains most parameters and compiler switches (for program flow and
 * options) as preprocessor defines, definitions of global structures,
 * declaration of global variables and function declarations for main.c.
 *
 * The main.h file is created from a shell script that is called from the
 * Makefile __before__ compilation. There is a main_template.h file which
 * contains placeholder that are filled according to the values in
 * parameters.sh. To keep this file short, we do not comment and document the
 * preprocessor directives. An extensive documentation of the parameter.sh file
 * is found in a separate documentation file.
 *
 * @see doc_parameters.md
 */

#define PI                      (3.141592653589793238462643383279)
#define TWOPI                   (6.283185307179586476925286766559)
#define MAX(x, y)               (((x) > (y)) ? (x) : (y))
#define MIN(x, y)               (((x) > (y)) ? (y) : (x))

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

#define RK4                 (0)
#define DOPRI853            (1)
#define INTEGRATION_METHOD  (_IM_)

#define VERSION_CONTROL_NONE    (0)
#define VERSION_CONTROL_HG      (1)
#define VERSION_CONTROL_GIT     (2)
#define VERSION_CONTROL         (_VC_)

#define PSI_ELLIPTIC            (0)
#define PSI_HYPERBOLIC          (1)
#define PSI_PARABOLIC           (2)
#define PSI_METHOD              (_PM_)

#define IC_FROM_INTERNAL_FUNCTION   (0)
#define IC_FROM_H5_FILE             (1)
#define IC_FROM_DAT_FILE            (2)
#define IC_FROM_BUNCH_DAVIES        (3)
#define INITIAL_CONDITIONS          (_IC_)

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

#define DATAPATH                ("_PATH_")
#define INITIAL_DATAPATH        ("_IPATH_")
#define WRITE_OUT_BUFFER_NUMBER (_BUF_)
#define POWER_SPECTRUM_BINS     (_BINS_)
#define TIME_STEP_SKIPS         (_TSKIP_)
#define STRIDE_X                (_XSKIP_)
#define STRIDE_Y                (_YSKIP_)
#define STRIDE_Z                (_ZSKIP_)
#define THREAD_NUMBER           (_THREADS_)
#define FFTW_DEFAULT_FLAG       (_FFTW_)
#define _FILTER_

#define GRIDPOINTS_X            (_GPX_)
#define GRIDPOINTS_Y            (_GPY_)
#define GRIDPOINTS_Z            (_GPZ_)
#define SPATIAL_LOWER_BOUND_X   (_LX_)
#define SPATIAL_UPPER_BOUND_X   (_UX_)
#define SPATIAL_LOWER_BOUND_Y   (_LY_)
#define SPATIAL_UPPER_BOUND_Y   (_UY_)
#define SPATIAL_LOWER_BOUND_Z   (_LZ_)
#define SPATIAL_UPPER_BOUND_Z   (_UZ_)

#define DELTA_T                 (_DELT_)
#define INITIAL_TIME            (_TI_)
#define FINAL_TIME              (_TF_)
#define MAX_STEPS               (_MAXSTEPS_)
#define MINIMAL_DELTA_T         (_MINDELT_)

#define SEED                    (_SEED_)

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

#define SMALLEST_SCALING        (_MINSCAL_)
#define LARGEST_SCALING         (_MAXSCAL_)
#define BETA                    (_BETA_) // ALPHA = 1.0/8.0 - BETA * 0.2
#define SAFE                    (_SAFE_)
#define RELATIVE_TOLERANCE      (_RELTOL_)
#define ABSOLUTE_TOLERANCE      (_ABSTOL_)
#define MAX_DT_HUBBLE_FRACTION  (_HFRAC_)

/**
 * @brief Holds the values necessary to describe one dimension of the grid.
 *
 * All values describe a single dimension of the simulation volume only.
 */
struct grid_dimension
{
    size_t N; ///< number of gridpoints in real space
    size_t M; ///< number of gridpoints in Fourier space
    double a; ///< lower bound of the interval
    double b; ///< upper bound of the interval
    double k; ///< used to compute k vectors: k = 2 pi I / (b - a)
    double k2; ///< used to compute k^2: k2 = k*k = -4 pi^2 / L^2
    size_t stride; ///< strides for output
    size_t outN; ///< number of output points in this dimension
};

/**
 * @brief Parameters related to the time evolution.
 */
struct timing
{
    size_t Nt; ///< number of timesteps (only relevant for fixed step size)
    double dt; ///< size of the (initial) timestep delta t
    double ti; ///< initial time
    double tf; ///< final time
    double t;  ///< current time (constantly updated during simulation)
};

/**
 * @brief A collection of dataset IDs for a single field.
 */
struct datasets
{
    hsize_t field;  ///< id for the field
    hsize_t mean;   ///< id for the mean of the field
    hsize_t var;    ///< id for the variance of the field
    hsize_t dfield; ///< id for the time derivative of the field
    hsize_t dmean;  ///< id for the mean of the time derivative of the field
    hsize_t dvar;   ///< id for the variance of the time derivative of the field
};

/**
 * @brief All parameters related to file IO operations
 */
struct file_parameters
{
    hsize_t id;               ///< h5 file id of the output file
    struct datasets dset_phi; ///< h5 dataset ids for the scalar field phi
    struct datasets dset_psi; ///< h5 dataset ids for the perturbation psi
    struct datasets dset_rho; ///< h5 dataset ids for the energy density rho
    hsize_t dset_powspec;     ///< h5 dataset id for the power spectrum
    hsize_t dset_time;        ///< h5 dataset id for the time parameter
    hsize_t dset_a;           ///< h5 dataset id for the scaling parameter a
    size_t index;             ///< current index within the buffers
    size_t buf_size;          ///< size of the buffer
    size_t skip;              ///< timesteps to skip in between write outs
    size_t bins_powspec;      ///< number of bins for the power spectrum
};

/**
 *  @brief Collection of all parameters.
 *
 *  @note Throughout the source code holds:
 *  Nx = number of grid points in the x direction
 *  Ny = number of grid points in the y direction
 *  Nz = number of grid points in the z direction
 *  N  = number of gridpoints for the whole spatial grid = Nx * Ny * Nz
 *  N2 = 2 * N, N3 = 3 * N, N2p = 2 * N + 2, N3p = 3 * N + 2
 *  Nall = size of field = 4 * N + 2
 *  Ntot = number of scalar equations; depends on the used method
 *  elliptic: Ntot = 2 * N + 1 (order: phi, dphi, a)
 *  parabolic: Ntot = 3 * N + 1 (oder: phi, dphi, psi, a)
 *  hyperbolic: Ntot = 4 * N + 1 (oder: phi, dphi, psi, dpsi, a)
 */
struct parameters
{
    struct grid_dimension x; ///< specification of the x direction
    struct grid_dimension y; ///< specification of the y direction
    struct grid_dimension z; ///< specification of the z direction
    size_t N; ///< number of spatial gridpoints: N=Nx*Ny*Nz
    size_t Ntot; ///< number of scalar equations evolved depending on PSI_METHOD
    size_t Nall; ///< size of field: Nall=4*N+2
    size_t outN; ///< number of spatial gridpoints for output
    size_t M; ///< number of gridpoints in Fourier space
    size_t dim; ///< dimensions of the simulation (1, 2 or 3)
    struct timing t; ///< time evolution parameters
    struct file_parameters file; ///< file IO parameters
};

/**
 * @brief Temporal arrays used for computations of gradients, ffts and the like.
 */
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

/**
 * @brief Grids for ksq and kx, ky, kz, i.e. the k vector and its square.
 */
struct k_grid
{
    double *sq;
    double *x;
    double *y;
    double *z;
};

extern struct parameters pars;
extern double *time_buf;

// contain phi, dphi, a, psi, dpsi (only first Ntot entries are integrated)
extern double *field;
extern double *dfield;
extern double *field_new;
extern double *dfield_new;

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

extern double *a_buf;

extern double *rho;
extern double *rho_buf;
extern double rho_mean;
extern double rho_var;
extern double *rho_mean_buf;
extern double *rho_var_buf;

// only needed for PSI_METHOD == PSI_HYPERBOLIC, but always defined
extern double *pressure;
extern double pressure_mean;

extern double *pow_spec;
extern double *pow_spec_buf;

extern double *filter;

extern struct k_grid kvec;

extern struct temporary tmp;

extern fftw_plan p_fw;
extern fftw_plan p_bw;

extern double fftw_time_exe;
extern double fftw_time_plan;
extern double filter_time;
extern double poisson_time;
extern double h5_time_write;

#ifdef SHOW_TIMING_INFO
double get_wall_time();
#endif

#endif
