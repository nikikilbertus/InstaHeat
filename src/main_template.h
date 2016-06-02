#ifndef __MAIN__
#define __MAIN__

#include <complex.h>
#include <fftw3.h>
#include "hdf5.h"

/**
 * @file main.h
 *
 * @brief Contains most parameters and compiler switches (for program flow and
 * options) as preprocessor defines, definitions of global structures,
 * declaration of global variables and function declarations for `main.c`.
 *
 * The `main.h` file is created by a shell script that is called from the
 * Makefile __before__ compilation. There is a main_template.h file which
 * contains placeholder that are replaced according to the values in
 * `parameters.sh`. An extensive documentation of the `parameter.sh` file is
 * found in a separate documentation file `doc_parameters.md` in the root
 * directory. One should not look for documentation and explanation here.
 *
 * @see `doc_parameters.md` in the root directory.
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

#ifdef SHOW_RUNTIME_INFO
#define INFO(f) do {\
        (f); \
    } while (0)
#else
#define INFO(f)
#endif

#ifdef SHOW_TIMING_INFO
#define TIME(f) do {\
        (f); \
    } while (0)
#else
#define TIME(f)
#endif

#define RK4                 (0)
#define DOPRI853            (1)
#define INTEGRATION_METHOD  (_IM_)

#define VERSION_CONTROL_NONE    (0)
#define VERSION_CONTROL_HG      (1)
#define VERSION_CONTROL_GIT     (2)
#define VERSION_CONTROL         (_VC_)

#define IC_FROM_INTERNAL_FUNCTION       (0)
#define IC_FROM_H5_FILE                 (1)
#define IC_FROM_DAT_FILE_WITH_PSI       (2)
#define IC_FROM_DAT_FILE_WITHOUT_PSI    (3)
#define IC_FROM_BUNCH_DAVIES            (4)
#define INITIAL_CONDITIONS              (_IC_)

#if INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITH_PSI \
    || INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITHOUT_PSI
    #define IC_FROM_DAT_FILE
#endif

#define _OPHI_
#define _ODPHI_
#define _OPSI_
#define _ODPSI_
#define _ORHO_
#define _OPHIPS_
#define _OPSIPS_
#define _ORHOPS_
#define _OPHIS_
#define _ODPHIS_
#define _OPSIS_
#define _ODPSIS_
#define _ORHOS_
#define _OCSTR_

#if defined(OUTPUT_PHI) || defined(OUTPUT_DPHI) || defined(OUTPUT_PSI) \
    || defined(OUTPUT_DPSI)|| defined(OUTPUT_RHO)
    #define LARGE_OUTPUT
#endif

#if defined(OUTPUT_PHI_PS) || defined(OUTPUT_DPHI_PS) || defined(OUTPUT_PSI_PS) \
    || defined(OUTPUT_DPSI_PS) || defined(OUTPUT_RHO_PS)
    #define OUTPUT_PS
#endif

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
#define MASS                    (_M_)
#define MASS_KARSTEN            (_MKARSTEN_)
#define INFLATON_MASS           (_MINFL_)
#define COUPLING                (_COUPLING_)
#define LAMBDA                  (_LAMBDA_)
#define A_INITIAL               (_A_)
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
 * @brief Specifies a dataset in the output h5 file.
 */
struct output
{
    size_t dim;     ///< length of output on a single time slice
    hsize_t id;     ///< id for the dataset in the h5 file
    double *buf;    ///< buffer WRITE_OUT_BUFFER_NUMBER many time slices
    double *tmp;    ///< the output values on the current time slice
};

/**
 * @brief All parameters related to file IO operations
 */
struct file_parameters
{
    hsize_t id;        ///< h5 file id of the output file
    size_t index;      ///< current index within the buffers
    size_t buf_size;   ///< size of the buffer
    size_t skip;       ///< timesteps to skip in between write outs
};

/**
 *  @brief Collection of all parameters.
 *
 *  @note Throughout the source code holds:
 *  Nx = number of grid points in the x direction
 *  Ny = number of grid points in the y direction
 *  Nz = number of grid points in the z direction
 *  N  = number of gridpoints for the whole spatial grid = Nx * Ny * Nz
 *  Next  = (extended) number of doubles in complex grid = Nx * Ny * (Nz + 2)
 *  Ntot = number of scalar equations in the integration: 4 * N + 4 * Next + 1
 *  (order: phi, dphi, psi, dpsi, h1, h2, dh1, dh2, a)
 */
struct parameters
{
    struct grid_dimension x; ///< specification of the x direction
    struct grid_dimension y; ///< specification of the y direction
    struct grid_dimension z; ///< specification of the z direction
    size_t N; ///< number of spatial gridpoints: N=Nx*Ny*Nz
    size_t Next; ///< number of doubles in complex grid: N=Nx*Ny*(Nz+2)
    size_t Ntot; ///< number of scalar equations
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
    double  *xphi; ///< The x derivative of $$\phi$$ in real space
    double  *yphi; ///< The y derivative of $$\phi$$ in real space
    double  *zphi; ///< The z derivative of $$\phi$$ in real space
    complex *phic; ///< The inflaton field $$\phi$$ in real space
    complex *xphic; ///< The x derivative of $$\phi$$ in Fourier space
    complex *yphic; ///< The y derivative of $$\phi$$ in Fourier space
    complex *zphic; ///< The z derivative of $$\phi$$ in Fourier space
    double  *grad; ///< The _squared_ gradient of $$\phi$$ in real space
    double  *lap; ///< The Laplacian of $$\phi$$ in real space
    double  *deltarho; ///< $$\delta \rho = \rho - <\rho>$$ in real space
    double  *f; ///< Various purposes (real space)
    complex *deltarhoc; ///< $$\delta \rho = \rho - <\rho>$$ in Fourier space
    complex *fc; ///< Various purposes (Fourier space)
    complex *psic; ///< The metric perturbation $$\psi$$ in Fourier space
    complex *dpsic; ///< The derivative $$\dot{\psi}$$ in Fourier space
};

/**
 * @brief Grids for ksq and kx, ky, kz, i.e. the k vector and its square.
 */
struct k_grid
{
    double *sq; ///< The sqaure of the k vector
    double *x; ///< The x direction of the k vector (with zeros at N/2)
    double *y; ///< The y direction of the k vector (with zeros at N/2)
    double *z; ///< The z direction of the k vector (with zeros at N/2)
    double *xf; ///< The x direction of the k vector
    double *yf; ///< The y direction of the k vector
    double *zf; ///< The z direction of the k vector
};

/**
 * @brief Collection of variables monitoring timing and function calls.
 */
struct monitor
{
    /**
     * @brief Number of calls to `mk_rhs(const double t, double *f, double
     * *result)`
     */
    size_t calls_rhs;
    double fftw_time_exe; ///< Total wall clock time for fft execution
    double fftw_time_plan; ///< Total wall clock time for fftw planning
    double filter_time; ///< Total wall clock time for filtering
    double poisson_time; ///< Total wall clock time for `mk_psi(double *f)`
    double h5_time_write; ///< Total wall clock time for write out
    double copy_buffer_time; ///< Total wall clock time for copying buffers
    double cstr_time; ///< Total wall clock time for computing constraints
    double smry_time; ///< Total wall clock time for computing summaries
    double stt_time; ///< Total wall clock time for computing S_{ij}^{TT}
};

extern struct parameters pars; ///< Only instance of the parameters_t struct

// contain phi, dphi, psi, dpsi, h1, h2, a
extern double *field;
extern double *dfield;
extern double *field_new;
extern double *dfield_new;
extern struct output phi;
extern struct output dphi;
extern struct output psi;
extern struct output dpsi;
#define SUMMARY_VALUES      (4) ///< mean, variance, minimum, maximum
extern struct output phi_smry;
extern struct output dphi_smry;
extern struct output psi_smry;
extern struct output dpsi_smry;
extern struct output h1_smry;
extern struct output h2_smry;
extern struct output t_out;
extern struct output a_out;
extern double *rho;
extern double rho_mean;
extern struct output rho_out;
extern struct output rho_smry;
extern double *pressure;
extern double pressure_mean;
extern struct output phi_ps;
extern struct output psi_ps;
extern struct output rho_ps;
#define NUMBER_CONSTRAINTS  (4) ///< Hamiltonian and momentum: l2 and l\infty
extern struct output cstr;
extern struct output gw;
extern double *filter;
extern struct k_grid kvec;
extern struct temporary tmp;
extern fftw_plan p_fw;
extern fftw_plan p_bw;
extern struct monitor mon;

#ifdef SHOW_TIMING_INFO
double get_wall_time();
#endif

#endif
