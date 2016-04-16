#ifndef __DOPRI853_STEPPER__
#define __DOPRI853_STEPPER__

#include "main.h"
#if INTEGRATION_METHOD == DOPRI853

/**
 * @file dopri853.h
 * @brief Typedefs and function declarations for `dopri853.c` and for
 * `dopri853_constants.c`.
 * @see <a href="http://numerical.recipes">Numerical Recipes</a>
 */

/**
 * @brief Holds the Butcher tableaux for the Dormand Prince integration routine.
 * @note This is just a huge list of constant double values. Do not change
 * anything here.
 */
struct dopri853_constants
{
    const double c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c14,c15,c16,
    b1,b6,b7,b8,b9,b10,b11,b12,
    bhh1,bhh2,bhh3,
    er1,er6,er7,er8,er9,er10,er11,er12,
    a21,
    a31,a32,
    a41,a43,
    a51,a53,a54,
    a61,a64,a65,
    a71,a74,a75,a76,
    a81,a84,a85,a86,a87,
    a91,a94,a95,a96,a97,a98,
    a101,a104,a105,a106,a107,a108,a109,
    a111,a114,a115,a116,a117,a118,a119,a1110,
    a121,a124,a125,a126,a127,a128,a129,a1210,a1211,
    a141,a147,a148,a149,a1410,a1411,a1412,a1413,
    a151,a156,a157,a158,a1511,a1512,a1513,a1514,
    a161,a166,a167,a168,a169,a1613,a1614,a1615,
    d41,d46,d47,d48,d49,d410,d411,d412,d413,d414,d415,d416,
    d51,d56,d57,d58,d59,d510,d511,d512,d513,d514,d515,d516,
    d61,d66,d67,d68,d69,d610,d611,d612,d613,d614,d615,d616,
    d71,d76,d77,d78,d79,d710,d711,d712,d713,d714,d715,d716;
};

/**
 * @brief Holds intermediate evaluations of the right hand side and errors for
 * the Dormand Prince integrator.
 *
 * Holds pointers to memory blocks for the intermediate evaluations of
 * the right hand side of the partial differential equation (as determined in
 * mk_rhs(const double t, double *f, double *result) in toolbox.c as
 * well as memory blocks for the error estimates (5th and 3rd order).
 */
struct dopri853_values
{
        double *k2, *k3, *k4, *k5, *k6, *k7, *k8, *k9, *k10, *k_tmp;
        double *yerr, *yerr2;
};

/**
 * @brief Holds parameters for the Dormand Prince integrator.
 *
 * Most of these parameters should be self explanatory by their names. For more
 * information see <a href="http://numerical.recipes">Numerical Recipes</a>.
 */
struct dopri853_control
{
    double t; ///< The current time
    double t_old; ///< The previous time (on last time slice)
    double ti; ///< The initial time
    double tf; ///< The final time
    double dt; ///< The time step size
    double dt_did; ///< The previously used time step size
    double dt_next; ///< The proposed next time step size
    double dt_min; ///< The minimal permissible time step size
    size_t max_steps; ///< The maximal number of steps
    int n_stp; ///< The number of performed steps
    int n_ok; ///< The number of successful steps
    int n_bad; ///< The number of unsuccessful steps
    double beta; ///< Internal parameter for the error estimates
    double alpha; ///< Internal parameter for the error estimates
    double safe; ///< Internal parameter for the error estimates
    double minscale; ///< Minimal permissible rescaling of the time step size
    double maxscale; ///< Maximal permissible rescaling of the time step size
    double a_tol; ///< Absolute tolerance
    double r_tol; ///< Relative tolerance
    double err_old; ///< The previous error (on the last time slice)
    int reject; ///< Flag whether time step is rejected or accepted
    double eps; ///< Epsilon value for comparisons
};

extern struct dopri853_constants dpc; ///< Dormand Prince Butcher tableaux constants
extern struct dopri853_values dpv; ///< Intermediate fields and errors
extern struct dopri853_control dp; ///< Dormand Prince parameters

void initialize_dopri853();
void run_dopri853();
int perform_step(const double dt_try);
void try_step(const double dt);
double error(const double dt);
int success(const double err, double *dt);
void allocate_dopri853_values();
void free_dopri853_values();

#endif
#endif
