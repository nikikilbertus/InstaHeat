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
 * contains placeholders that are replaced according to the values in
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
/* #define CHECK_FOR_NAN // not recommended (performance) */
/* #define ENABLE_PROFILER // only recommended for debugging */
/* #define DEBUG // only recommended for debugging (huge output!) */

#ifdef SHOW_RUNTIME_INFO
#define INFO(f) do {\
        (f); \
    } while (0)
#else
#define INFO(f)
#endif

#ifdef ENABLE_TIMING
#define TIME(f) do {\
        (f); \
    } while (0)
#else
#define TIME(f)
#endif

#define RK4                 (0)
#define DOPRI853            (1)
#define INTEGRATION_METHOD  (__INTEGRATION_METHOD__)

#define VERSION_CONTROL_NONE    (0)
#define VERSION_CONTROL_HG      (1)
#define VERSION_CONTROL_GIT     (2)
#define VERSION_CONTROL         (__VERSION_CONTROL__)

#define IC_FROM_INTERNAL_FUNCTION       (0)
#define IC_FROM_H5_FILE                 (1)
#define IC_FROM_DAT_FILE_WITH_PSI       (2)
#define IC_FROM_DAT_FILE_WITHOUT_PSI    (3)
#define IC_FROM_BUNCH_DAVIES            (4)
#define INITIAL_CONDITIONS              (__INITIAL_CONDITIONS__)

#if INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITH_PSI \
    || INITIAL_CONDITIONS == IC_FROM_DAT_FILE_WITHOUT_PSI
    #define IC_FROM_DAT_FILE
#endif

#define __OUTPUT_PHI__
#define __OUTPUT_DPHI__
#define __OUTPUT_PSI__
#define __OUTPUT_DPSI__
#define __OUTPUT_RHO__
#define __OUTPUT_PHI_PS__
#define __OUTPUT_PSI_PS__
#define __OUTPUT_RHO_PS__
#define __OUTPUT_PHI_SMRY__
#define __OUTPUT_DPHI_SMRY__
#define __OUTPUT_PSI_SMRY__
#define __OUTPUT_DPSI_SMRY__
#define __OUTPUT_RHO_SMRY__
#define __OUTPUT_PRESSURE_SMRY__
#define __OUTPUT_H1_SMRY__
#define __OUTPUT_H2_SMRY__
#define __OUTPUT_CONSTRAINTS__

#if defined(OUTPUT_PHI) || defined(OUTPUT_DPHI) || \
    defined(OUTPUT_PSI) || defined(OUTPUT_DPSI) || \
    defined(OUTPUT_RHO)
    #define LARGE_OUTPUT
#endif

#if defined(OUTPUT_PHI_PS) || defined(OUTPUT_DPHI_PS) || \
    defined(OUTPUT_PSI_PS) || defined(OUTPUT_DPSI_PS) || \
    defined(OUTPUT_RHO_PS)
    #define OUTPUT_PS
#endif

#if defined(OUTPUT_PHI_SMRY) || defined(OUTPUT_DPHI_SMRY) || \
    defined(OUTPUT_PSI_SMRY) || defined(OUTPUT_DPSI_SMRY) || \
    defined(OUTPUT_RHO_SMRY)
    #define OUTPUT_SMRY
#endif

#define DATAPATH                ("__DATAPATH__")
#define INITIAL_DATAPATH        ("__INITIAL_DATAPATH__")
#define WRITE_OUT_BUFFER_NUMBER (__WRITE_OUT_BUFFER_NUMBER__)
#define POWER_SPECTRUM_BINS     (__POWER_SPECTRUM_BINS__)
#define TIME_STEP_SKIPS         (__TIME_STEP_SKIPS__)
#define STRIDE_X                (__STRIDE_X__)
#define STRIDE_Y                (__STRIDE_Y__)
#define STRIDE_Z                (__STRIDE_Z__)
#define THREAD_NUMBER           (__THREAD_NUMBER__)
#define FFTW_DEFAULT_FLAG       (__FFTW_DEFAULT_FLAG__)
#define GRIDPOINTS_X            (__GRIDPOINTS_X__)
#define GRIDPOINTS_Y            (__GRIDPOINTS_Y__)
#define GRIDPOINTS_Z            (__GRIDPOINTS_Z__)
#define SPATIAL_LOWER_BOUND_X   (__SPATIAL_LOWER_BOUND_X__)
#define SPATIAL_UPPER_BOUND_X   (__SPATIAL_UPPER_BOUND_X__)
#define SPATIAL_LOWER_BOUND_Y   (__SPATIAL_LOWER_BOUND_Y__)
#define SPATIAL_UPPER_BOUND_Y   (__SPATIAL_UPPER_BOUND_Y__)
#define SPATIAL_LOWER_BOUND_Z   (__SPATIAL_LOWER_BOUND_Z__)
#define SPATIAL_UPPER_BOUND_Z   (__SPATIAL_UPPER_BOUND_Z__)
#define DELTA_T                 (__DELTA_T__)
#define INITIAL_TIME            (__INITIAL_TIME__)
#define FINAL_TIME              (__FINAL_TIME__)
#define MAX_STEPS               (__MAX_STEPS__)
#define MINIMAL_DELTA_T         (__MINIMAL_DELTA_T__)
#define SEED                    (__SEED__)
#define MASS                    (__MASS__)
#define MASS_KARSTEN            (__MASS_KARSTEN__)
#define INFLATON_MASS           (__INFLATON_MASS__)
#define COUPLING                (__COUPLING__)
#define LAMBDA                  (__LAMBDA__)
#define A_INITIAL               (__A_INITIAL__)
#define SMALLEST_SCALING        (__SMALLEST_SCALING__)
#define LARGEST_SCALING         (__LARGEST_SCALING__)
#define BETA                    (__BETA__) // ALPHA = 1.0/8.0 - BETA * 0.2
#define SAFE                    (__SAFE__)
#define RELATIVE_TOLERANCE      (__RELATIVE_TOLERANCE__)
#define ABSOLUTE_TOLERANCE      (__ABSOLUTE_TOLERANCE__)
#define MAX_DT_HUBBLE_FRACTION  (__MAX_DT_HUBBLE_FRACTION__)
#define BUNCH_DAVIES_CUTOFF     (__BUNCH_DAVIES_CUTOFF__)
#define MAX_RUNTIME             (__MAX_RUNTIME__)
#define __ENABLE_FFT_FILTER__
#define __ENABLE_GW__
#define __ENABLE_FOLLOWUP__
#define __ENABLE_TIMING__

#ifndef ENABLE_GW
    #undef OUTPUT_H1_SMRY
    #undef OUTPUT_H2_SMRY
#endif

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
    double k; ///< used to compute k vectors: k = 2 pi / (b - a)
    double k2; ///< used to compute k^2: k2 = k*k = -4 pi^2 / L^2
    size_t stride; ///< strides for output
    size_t outN; ///< number of output points in this dimension
};

/**
 * @brief Holds all parameters related to the time evolution.
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
    double *buf;    ///< buffer for `WRITE_OUT_BUFFER_NUMBER` many time slices
    double *tmp;    ///< the output values on the current time slice
};

/**
 * @brief Holds all parameters related to file IO operations.
 */
struct file_parameters
{
    hsize_t id;        ///< h5 file id of the output file
    size_t index;      ///< current index within the buffers for time evolution
    size_t buf_size;   ///< size of the buffer
    size_t skip;       ///< timesteps to skip in between write outs
};

/**
 *  @brief An overall collection of all simulation parameters.
 *
 *  @note Throughout the source code holds:
 *  Nx = number of grid points in the x direction
 *  Ny = number of grid points in the y direction
 *  Nz = number of grid points in the z direction
 *  N  = number of gridpoints for the whole spatial grid = Nx * Ny * Nz
 *  Next  = (extended) number of doubles in complex grid = Nx * Ny * (Nz + 2)
 *  Ntot = number of scalar equations in the integration: 4 * N + 4 * Next + 1
 *  (order: \f$\phi\f$ (N), \f$\dot{\phi}\f$ (N), \f$\psi\f$ (N),
 *  \f$\dot{\psi}\f$ (N), \f$h_1\f$ (Next), \f$h_2\f$ (Next), \f$\dot{h}_1\f$
 *  (Next), \f$\dot{h}_2\f$ (Next), \f$a\f$ (1))
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
    size_t bunch_davies_cutoff; ///< cutoff of the initial Bunch Davies spectrum
    struct timing t; ///< time evolution parameters
    struct file_parameters file; ///< file IO parameters
    size_t max_runtime; ///< The maximal overall runtime of the program
};

/**
 * @brief Memory for DFTs and temporary quantities.
 */
struct temporary
{
    double  *xphi; ///< the x derivative of \f$\phi\f$ in real space
    double  *yphi; ///< the y derivative of \f$\phi\f$ in real space
    double  *zphi; ///< the z derivative of \f$\phi\f$ in real space
    complex *phic; ///< the inflaton field \f$\phi\f$ in real space
    complex *xphic; ///< the x derivative of \f$\phi\f$ in Fourier space
    complex *yphic; ///< the y derivative of \f$\phi\f$ in Fourier space
    complex *zphic; ///< the z derivative of \f$\phi\f$ in Fourier space
    double  *grad; ///< the _squared_ gradient of \f$\phi\f$ in real space
    double  *lap; ///< the Laplacian of \f$\phi\f$ in real space
    double  *deltarho; ///< \f$\delta \rho = \rho - \langle \rho \rangle\f$ in real space
    double  *f; ///< various purposes (real space)
    complex *deltarhoc; ///< \f$\delta \rho = \rho - \langle \rho \rangle\f$ in Fourier space
    complex *fc; ///< various purposes (Fourier space)
    complex *psic; ///< the metric perturbation \f$\psi\f$ in Fourier space
    complex *dpsic; ///< the derivative \f$\dot{\psi}\f$ in Fourier space
};

/**
 * @brief Grids for ksq and kx, ky, kz, i.e. the k vector and its square.
 */
struct k_grid
{
    double *sq; ///< the sqaure of the k vector
    double *x; ///< the x component of the k vector (with zeros at N/2)
    double *y; ///< the y component of the k vector (with zeros at N/2)
    double *z; ///< the z component of the k vector (with zeros at N/2)
    double *xf; ///< the x component of the k vector
    double *yf; ///< the y component of the k vector
    double *zf; ///< the z component of the k vector
    double k2_max; ///< The maximal square k vector on the grid.
    double k_max; ///< The sqrt of the maximal square k vector on the grid.
};

/**
 * @brief Collection of variables monitoring timing and function calls.
 */
struct monitor
{
    /**
     * @brief Number of calls to `mk_rhs(const double t, double *f, double
     * *result)` in `toolbox.c`
     */
    size_t calls_rhs;
    double total; ///< Total wall clock time for the overall program execution
    double fftw_exe; ///< Total wall clock time for fft execution
    double fftw_plan; ///< Total wall clock time for fftw planning
    double filter; ///< Total wall clock time for filtering
    double elliptic; ///< Total wall clock time for `mk_psi(double *f)`
    double integration; ///< Total wall clock time for the integration
    double h5_write; ///< Total wall clock time for write out to disk
    double cpy_buffers; ///< Total wall clock time for copying buffers
    double cstr; ///< Total wall clock time for computing constraints
    double smry; ///< Total wall clock time for computing summaries
    double gw_sources; ///< Total wall clock time for computing \f$S_{ij}^{TT}\f$
};

extern struct parameters pars; ///< Only instance of the struct `parameters`

// contain phi, dphi, psi, dpsi, h1, h2, a
extern double *field; ///< An array bundling all fields that are evolved
extern double *dfield; ///< An array bundling temporal derivatives of all fields
extern double *field_new; ///< An array just like `field` for copying
extern double *dfield_new; ///< An array just like 'dfield' for copying
extern struct output phi; ///< The output struct for \f$\phi\f$
extern struct output dphi; ///< The output struct for \f$\dot{\phi}\f$
extern struct output psi; ///< The output struct for \f$\psi\f$
extern struct output dpsi; ///< The output struct for \f$\dot{\psi}\f$
#define SUMMARY_VALUES (4) ///< Summaries contain mean, variance, minimum, maximum
extern struct output phi_smry; ///< The output struct for the summary of \f$\phi\f$
extern struct output dphi_smry; ///< The output struct for the summary of \f$\dot{\phi}\f$
extern struct output psi_smry; ///< The output struct for the summary of \f$\psi\f$
extern struct output dpsi_smry; ///< The output struct for the summary of \f$\dot{\psi}\f$
extern struct output h1_smry; ///< The output struct for the summary of \f$h_1\f$
extern struct output h2_smry; ///< The output struct for the summary of \f$h_2\f$
extern struct output t_out; ///< The output struct for the summary of \f$t\f$
extern struct output a_out; ///< The output struct for the summary of \f$a\f$
extern double *rho; ///< An array for the energy density \f$\rho\f$
extern double rho_mean; ///< The mean value of the energy density \f$\rho\f$
extern struct output rho_out; ///< The output struct for \f$\rho\f$
extern struct output rho_smry; ///< The output struct for the summary of \f$\rho\f$
extern double *pressure; ///< An array for the pressure \f$p\f$
extern double pressure_mean; ///< The mean value of the pressure \f$p\f$
extern struct output p_smry; ///< The output struct for the summary of \f$p\f$
extern struct output phi_ps; ///< The output struct for the power spectrum of \f$\phi\f$
extern struct output psi_ps; ///< The output struct for the power spectrum of \f$\psi\f$
extern struct output rho_ps; ///< The output struct for the power spectrum of \f$\rho\f$
/**
 * @brief We compute the Hamiltonian and the momentum constraint in \f$l_2\f$
 * and \f$l_{\infty}\f$ norm
 */
#define NUMBER_CONSTRAINTS (4)
extern struct output cstr; ///< The output struct for the constraints
extern struct output gw; ///< The output struct for the gravitational wave power spectrum
extern double *filter; ///< A grid with the filter values for frequency filtering
extern struct k_grid kvec; ///< The only instance of the struct `k_grid`
extern struct temporary tmp; ///< The only instance of the struct `temporary`
extern fftw_plan p_fw; ///< FFTw3 plan for the Fourier transforms
extern fftw_plan p_bw; ///< FFTw3 plan for the inverse Fourier transforms
extern struct monitor mon; ///< The only instance of the struct `monitor`

double get_wall_time();

#endif
