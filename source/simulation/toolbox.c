#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "toolbox.h"
#include "main.h"
#include "io.h"
#include "setup.h"

/**
 * @file toolbox.c
 * @brief A collection of routines to compute the right hand side of the pde as
 * well as the desired output.
 */

static void assemble_gradient_squared();
static double potential(const double f);
static double potential_prime(const double f);
static void mk_constraints();
static void mk_power_spectrum(const fftw_complex *in, struct output out);
static void apply_filter_real(double *inout);
static void apply_filter_fourier(fftw_complex *phi_io, fftw_complex *dphi_io,
        fftw_complex *psi_io, fftw_complex *dpsi_io);
static void center(double *f, const size_t N);
static double mean(const double *f, const size_t N);
static void mean_var_min_max(const double *f, double *smry);
static double variance(const double mean, const double *f, const size_t N);
static void contains_nan(const double *f, const size_t N);
static void contains_nanc(const complex *f, const size_t N);

struct evolution_flags evo_flags = {.filter = 0,
                                    .compute_pow_spec = 0,
                                    .compute_cstr = 0};

/**
 * @brief Compute the right hand side of the pde, i.e. all first order temporal
 * derivatives.
 *
 * Depending on `PSI_METHOD` we are computing the temporal derivatives of different fields, and hence also evovle different fields in the integration routine:
 *
 * @param[in] t The current time
 * @param[in] All necessary fields bundled in one array
 * @param[out] The right hand side of the partial differential equation, i.e.
 * the first temporal derivatives of the fields
 *
 * - For `PSI_METHOD=PSI_ELLIPTIC` we are only evolving $$\phi$$,
 *   $$\dot{\phi}$$ and $$a$$ withthe integration routine and use an elliptic
 *   constraint equation to compute $$\psi$$ and $$\dot{\psi}$$ separately on
 *   each timeslice.
 * - For 'PSI_METHOD=PSI_PARABOLIC` we are evolving $$\phi$$, $$\dot{\phi}$$,
 *   $$\psi$$ and $$a$$ according to a parabolic constraint equation with the
 *   integration routine. $$\dot{\psi}$$ is computed separately on each
 *   timeslice.
 * - For 'PSI_METHOD=PSI_HYPERBOLIC` we are evolving $$\phi$$, $$\dot{\phi}$$,
 *   $$\psi$$, $$\dot{\psi}$$ and $$a$$ according to a hyperbolic constraint
 *   equation with the integration routine.
 */
void mk_rhs(const double t, double *f, double *result)
{
    mon.calls_rhs += 1;
    const size_t N = pars.N;
    const size_t N2 = 2 * N;
    const size_t N2p = N2 + 2;
    const size_t N3p = 3 * N + 2;
    const double a = f[N2];
    const double a2 = a * a;

    mk_gradient_squared_and_laplacian(f);
    mk_rho(f);
    #ifdef OUTPUT_CONSTRAINTS
    mk_constraints();
    #endif
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;

    // copy dphi in all cases
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        result[i] = f[N + i];
    }

    #ifndef EVOLVE_WITHOUT_PSI
    // parabolic: equation for dpsi; hyperbolic: copy dpsi and equation for ddpsi
    #if PSI_METHOD != PSI_ELLIPTIC
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        #if PSI_METHOD == PSI_PARABOLIC
        result[N2p + i] = -hubble * f[N2p + i] - 0.5 * (rho[i] - rho_mean) / h3
            + tmp.f[i] / (h3 * a2);
        #elif PSI_METHOD == PSI_HYPERBOLIC
        result[N2p + i] = f[N3p + i];
        result[N3p + i] = 0.5 * pressure[i] + (f[N2p + i] - 0.5) * pressure_mean
            - 4.0 * hubble * f[N3p + i];
        #endif
    }
    #else
    mk_psi(f);
    #endif
    #endif

    // equation for ddphi in all cases (psi & dpsi have to be provided first)
    double df, p, dp;
    #pragma omp parallel for private(df, p, dp)
    for (size_t i = 0; i < N; ++i) {
        df = f[N + i];
        p = f[N2p + i];
        #if PSI_METHOD != PSI_PARABOLIC
        dp = f[N3p + i];
        #else
        dp = result[N2p + i];
        #endif
        result[N + i] = (1.0 + 4.0 * p) * tmp.lap[i] / a2 -
            (h3 - 4.0 * dp) * df - (1.0 + 2.0 * p) * potential_prime(f[i]);
    }

    // update da
    result[N2] = a * hubble;
}

/**
 * @brief Compute the Laplacian and the sqaure gradient.
 *
 * @param[in] in An array containing the fields, in particular $$\phi$$ and
 * $$\dot{\phi}$$
 *
 * The Laplacian and the squared gradient as well as the partial derivatives of
 * $$\phi$$ are stored in global variables for reuse in other functions. When
 * the function returns:
 * `tmp.xphi` contains $$\partial_x \phi$$
 * `tmp.yphi` contains $$\partial_y \phi$$
 * `tmp.zphi` contains $$\partial_z \phi$$
 * `tmp.lap` contains the Laplacian $$\sum_{i=1}^3 \partial_i^2 \phi$$
 * `tmp.grad` contains the _squared_ gradient $$(\nabla \phi)^2$$
 * All the above values persist until the next call of `mk_rhs(const double t,
 * double *f, double *result)`
 * If `PSI_METHOD=PSI_PARABOLIC` or `OUTPUT_CONSTRAINTS` is defined,
 * additionally, `tmp.f` contains $$\Delta \psi$$, the Lagrangian of $$\psi$$.
 */
void mk_gradient_squared_and_laplacian(double *in)
{
    const size_t N = pars.N;
    const size_t M = pars.M;

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, in, tmp.phic);
        #if PSI_METHOD == PSI_PARABOLIC || defined(OUTPUT_CONSTRAINTS)
        const size_t N2p = 2 * N + 2;
        fftw_execute_dft_r2c(p_fw, in + N2p, tmp.psic);
        #endif
    TIME(mon.fftw_time_exe += get_wall_time());

    if (evo_flags.compute_pow_spec == 1) {
        mk_power_spectrum(tmp.phic, phi_ps);
    }

    complex pre;
    #pragma omp parallel for private(pre)
    for (size_t i = 0; i < M; ++i) {
        pre = tmp.phic[i] * I / N;
        tmp.xphic[i] = pre * kvec.x[i];
        tmp.yphic[i] = pre * kvec.y[i];
        tmp.zphic[i] = pre * kvec.z[i];
        tmp.phic[i] *= kvec.sq[i] / N;
        #if PSI_METHOD == PSI_PARABOLIC || defined(OUTPUT_CONSTRAINTS)
        tmp.psic[i] *= kvec.sq[i] / N;
        #endif
    }

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, tmp.xphic, tmp.xphi);
    if (pars.dim > 1) {
        fftw_execute_dft_c2r(p_bw, tmp.yphic, tmp.yphi);
        if (pars.dim > 2) {
            fftw_execute_dft_c2r(p_bw, tmp.zphic, tmp.zphi);
        }
    }
    fftw_execute_dft_c2r(p_bw, tmp.phic, tmp.lap);
    #if PSI_METHOD == PSI_PARABOLIC || defined(OUTPUT_CONSTRAINTS)
    fftw_execute_dft_c2r(p_bw, tmp.psic, tmp.f);
    #endif
    TIME(mon.fftw_time_exe += get_wall_time());
    assemble_gradient_squared();
}

/**
 * @brief Constructs the squared gradient from the (up to) three partial
 * spatial derivatives.
 *
 * The partial derivatives of $$\phi$$ computed in
 * `mk_gradient_squared_and_laplacian(double *in)` are indiviually squared,
 * added up and the result is stored in `tmp.grad`
 */
static void assemble_gradient_squared()
{
    const size_t N = pars.N;
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        tmp.grad[i] = tmp.xphi[i] * tmp.xphi[i];
        if (pars.dim > 1) {
            tmp.grad[i] += tmp.yphi[i] * tmp.yphi[i];
            if (pars.dim > 2) {
                tmp.grad[i] += tmp.zphi[i] * tmp.zphi[i];
            }
        }
    }
}

/**
 * @brief Compute the energy density $$\rho$$ and its average value. For
 * `PSI_METHOD=PSI_HYPERBOLIC` also compute the pressure and its average
 * value.
 *
 * @param[in] f An array containing the fields
 *
 * Everything is stored in global variables for reuse in other functions. When
 * the function returns:
 * `rho` contains the energy density $$\rho$$
 * `rho_mean` contains the average energy density $$< \rho >$$
 * If `PSI_METHOD=PSI_HYPERBOLIC`
 * `pressure` contains the pressure $$p$$
 * `pressure_mean` contains the average pressure $$< p >$$
 * All the above values persist until the next call of `mk_rhs(const double t,
 * double *f, double *result)`
 */
void mk_rho(const double *f)
{
    const size_t N = pars.N;
    const size_t N2 = 2 * N;
    const size_t N2p = N2 + 2;
    const double a = f[N2];
    const double a2 = a * a;
    rho_mean = 0.0;
    #if PSI_METHOD == PSI_HYPERBOLIC
    pressure_mean = 0.0;
    #endif

    double df, p, t1, t2;
    #pragma omp parallel for private(df, p, t1, t2) \
                                reduction(+: rho_mean, pressure_mean)
    for (size_t i = 0; i < N; ++i) {
        df = f[N + i];
        #if PSI_METHOD != PSI_ELLIPTIC
        p = f[N2p + i];
        t1 = (0.5 - p) * df * df;
        t2 = (0.5 + p) * tmp.grad[i] / a2;
        rho[i] = t1 + t2 + potential(f[i]);
            #if PSI_METHOD == PSI_HYPERBOLIC
            pressure[i] = t1 - t2 / 3.0 - potential(f[i]);
            pressure_mean += pressure[i];
            #endif
        #else
        rho[i] = (df * df + tmp.grad[i] / a2) / 2.0 + potential(f[i]);
        #endif
        rho_mean += rho[i];
    }
    rho_mean /= N;
    #if PSI_METHOD == PSI_HYPERBOLIC
    pressure_mean /= N;
    #endif
}

/**
 * @brief The potential of the scalar inflaton field $$\phi$$
 *
 * @param[in] f The field value where to evaluate the potential
 * @return The potential value at the given input
 */
static double potential(const double f)
{
    // higgs metastability potential
    /* double l = LAMBDA / 1.0e-10; */
    /* double a = 0.01 * l, b = 1.0 * l; */
    /* return f == 0.0 ? LAMBDA : */
    /*     (LAMBDA + a * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) * pow(f, 4) + */
    /*     b * pow(f, 6)); */

    // notch or step potential
    // LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5
    /* double lambda = 100.0; */
    /* return LAMBDA / (1.0 + exp(-lambda * f)); */

    // standard f squared potential
    return MASS * MASS * f * f / 2.0;
}

/**
 * @brief The derivative of the potential of the scalar inflaton field $$\phi$$
 *
 * @param[in] f The field value where to evaluate the derivative of the
 * potential
 * @return The value of the derivative of the potential at given input
 */
static double potential_prime(const double f)
{
    // higgs metastability potential
    /* double l = LAMBDA / 1.0e-10; */
    /* double a = 0.01 * l, b = 1.0 * l; */
    /* return f == 0.0 ? 0 : */
    /*     (4.0 * a * pow(f, 3) * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) - */
    /*     (0.1 * a * pow(f, 4) * ((f > 0.0) - (f < 0.0))) / (fabs(f) * log(10.0)) */
    /*     + 6.0 * b * pow(f, 5)); */

    // notch or step potential
    // LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5
    /* double lambda = 100.0; */
    /* double tmp = exp(lambda * f); */
    /* return LAMBDA * lambda * tmp / ((1.0 + tmp) * (1.0 + tmp)); */

    // standard f squared potential
    return MASS * MASS * f;
}

/**
 * @brief Monitor how well the Hamiltonian and momentum constraint are fulfilled.
 *
 * Computes the Hamiltonian and momentum constraint in a form such they should
 * give 0, if fulfilled exactly. Of these combinations save the l2 Norm as well
 * as the maximum asbolute value to `cstr.tmp` for output.
 *
 * @note So far only the Hamiltonian constraint is implemented! (4/14/2016)
 */
static void mk_constraints()
{
    TIME(mon.cstr_time -= get_wall_time());
    const size_t N = pars.N;
    const size_t N2 = 2 * N;
    const size_t N2p = 2 * N + 2;
    const size_t N3p = 3 * N + 2;
    const double a = field[N2];
    const double a2 = a * a;
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;
    const double phi_mean = mean(field, N);
    const double dphi_mean = mean(field + N, N);
    double ham, ham_l2 = 0.0, ham_max = 0.0;
    double mom, mom_l2 = 0.0, mom_max = 0.0;
    double tmp1;

    #pragma omp parallel for private(tmp1, ham, mom) \
        reduction(max: ham_max, mom_max) reduction(+: ham_l2, mom_l2)
    for (size_t i = 0; i < N; ++i) {
        tmp1 = hubble * field[N2p + i] + field[N3p + i];
        ham = tmp.f[i] / a2 - h3 * tmp1 - 0.5 * (rho[i] - rho_mean);
        mom = tmp1 - 0.5 * dphi_mean * (field[i] - phi_mean);
        ham_l2 += ham * ham;
        mom_l2 += mom * mom;
        ham_max = MAX(ham_max, fabs(ham));
        mom_max = MAX(mom_max, fabs(mom));
    }
    cstr.tmp[0] = ham_l2;
    cstr.tmp[1] = ham_max;
    cstr.tmp[2] = mom_l2;
    cstr.tmp[3] = mom_max;
    TIME(mon.cstr_time += get_wall_time());
}

/**
 * @brief Compute $$\psi$$ and $$\dot{\psi}$$
 *
 * @param[in] f An array containing the fields. Expects $$\phi$$,
 * $$\dot{\phi}$$ and $$a$$ right after each other in @p f.
 *
 * We use an elliptic equation from the Hamiltonian constraint combined with
 * the momentum contraint to compute $$\psi$$ and $$\dot{\psi}$$ from given
 * $$\phi$$, $$\dot{\phi}$$ and $$a$$.
 */
void mk_psi(double *f)
{
    //TODO[performance]: optimize be saving often used values (h3, dphiextra...)
    const size_t N = pars.N;
    const size_t M = pars.M;
    const size_t N2p = 2 * N + 2;
    const size_t N3p = 3 * N + 2;
    const double a = f[2 * N];
    const double a2 = a * a;
    const double hubble = sqrt(rho_mean / 3.0);

    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_c2r(p_bw, tmp.psic, f + N2p); */
    /* fftw_execute_dft_c2r(p_bw, tmp.dpsic, f + N3p); */
    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe += get_wall_time(); */
    /* mon.poisson_time += get_wall_time(); */
    /* #endif */

    // sophisticated version square first then average plus gradient
    /* #ifdef SHOW_TIMING_INFO */
    /* mon.poisson_time -= get_wall_time(); */
    /* #endif */

    /* double extra1 = 0.0; */
    /* double extra2 = 0.0; */
    /* #pragma omp parallel for reduction(+: extra1, extra2) */
    /* for (size_t i = 0; i < N; ++i) { */
    /*     tmp.f[i] = f[N + i] * tmp.xphi[i]; */
    /*     tmp.deltarho[i] = rho[i] - rho_mean; */
    /*     extra1 += f[N + i] * f[N + i]; */
    /*     extra2 += tmp.grad[i]; */
    /* } */
    /* const double extra = 0.5 * (extra1 - extra2 / a2) / N; */

    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_r2c(p_fw, tmp.f, tmp.fc); */
    /* fftw_execute_dft_r2c(p_fw, tmp.deltarho, tmp.deltarhoc); */
    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe += get_wall_time(); */
    /* #endif */

    /* tmp.fc[0] = 0.0; */
    /* tmp.psic[0] = 0.0; */
    /* tmp.dpsic[0] = 0.0; */
    /* #pragma omp parallel for */
    /* for (size_t i = 1; i < M; ++i) { */
    /*     if (fabs(kvec.x[i]) < DBL_EPSILON) { */
    /*         tmp.fc[i] = 0.0; */
    /*     } else { */
    /*         tmp.fc[i] /= kvec.x[i] * I; */
    /*     } */
    /*     if (fabs(kvec.sq[i] / a2 + extra) < 1.0e-14) { */
    /*         tmp.psic[i] = 0.0; */
    /*     } else { */
    /*         tmp.psic[i] = 0.5 * (tmp.deltarhoc[i] + 3.0 * hubble * tmp.fc[i]) / */
    /*             ((kvec.sq[i] / a2 + extra) * N); */
    /*     } */
    /*     tmp.dpsic[i] = 0.5 * tmp.fc[i] / N - hubble * tmp.psic[i]; */
    /* } */

    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_c2r(p_bw, tmp.psic, f + N2p); */
    /* fftw_execute_dft_c2r(p_bw, tmp.dpsic, f + N3p); */
    /* #ifdef SHOW_TIMING_INFO */
    /* mon.fftw_time_exe += get_wall_time(); */
    /* mon.poisson_time += get_wall_time(); */
    /* #endif */

    // simpler version
    TIME(mon.poisson_time -= get_wall_time());
    const double phi_mean = mean(f, N);
    const double dphi_mean = mean(f + N, N);
    double extra1 = 0.0;
    double extra2 = 0.0;
    #pragma omp parallel for reduction(+: extra1, extra2)
    for (size_t i = 0; i < N; ++i) {
        tmp.deltarho[i] = rho[i] - rho_mean;
        tmp.f[i] = dphi_mean * (f[i] - phi_mean);
        extra1 += f[N + i] * f[N + i];
        extra2 += tmp.grad[i];
    }
    const double extra = 0.5 * (extra1 - extra2 / a2) / N;
    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, tmp.deltarho, tmp.deltarhoc);
    fftw_execute_dft_r2c(p_fw, tmp.f, tmp.fc);
    TIME(mon.fftw_time_exe += get_wall_time());

    for (size_t i = 1; i < M; ++i) {
        tmp.phic[i] = 0.5 * (tmp.deltarhoc[i] +
                3 * hubble * tmp.fc[i]) / ((kvec.sq[i] / a2 + extra) * N);
    }
    tmp.phic[0] = 0.0;
    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, tmp.phic, f + N2p);
    TIME(mon.fftw_time_exe += get_wall_time());
    for (size_t i = 0; i < N; ++i) {
        f[N3p + i] = 0.5 * tmp.f[i] - hubble * f[N2p + i];
    }
    TIME(mon.poisson_time += get_wall_time());
}

/**
 * @brief Compute the power spectrum
 *
 * @param[in] in The Fourier amplitudes of the field
 * @param[in, out] out The output struct providing information about and memory
 * for the power spectrum of the field.
 *
 * The power spectrum is constructed by binning the Fourier modes according to
 * the size of their wave vectors. The number of bins is determined by
 * `POWER_SPECTRUM_BINS`
 */
static void mk_power_spectrum(const fftw_complex *in, struct output out)
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    const size_t N = pars.N;
    const size_t M = pars.M;
    const size_t bins = out.dim;

    const double k2_max = pars.x.k2 * (Nx/2) * (Nx/2) +
                    pars.y.k2 * (Ny/2) * (Ny/2) + pars.z.k2 * (Nz/2) * (Nz/2);

    #pragma omp parallel for
    for (size_t i = 0; i < bins; ++i) {
        out.tmp[i] = 0.0;
    }

    double pow2_tmp = 0.0;
    size_t idx;
    for (size_t i = 0; i < M; ++i) {
        if (fabs(kvec.z[i]) < DBL_EPSILON) {
            pow2_tmp = in[i] * conj(in[i]);
        } else {
            pow2_tmp = 2.0 * in[i] * conj(in[i]);
        }
        idx = (int)trunc(bins * sqrt(kvec.sq[i] / k2_max) - 1.0e-14);
        out.tmp[idx] += pow2_tmp / N;
    }
}

/**
 * @brief Apply a Fourier filter to each field of a given input to cutoff high
 * frequency modes
 *
 * @param[in, out] inout The field which we want to filter. Expect $$\phi$$ at
 * index 0, $$\dot{\phi}$ at index N, $$\psi$$ at index 2*N+2, $$\dot{\psi}$$
 * at index 3*N+2. All four are overwritten by their filtered results.
 * The highest modes of each field are cut off according to `filter_window(const
 * double x)` in `setup.c`.
 */
static void apply_filter_real(double *inout)
{
    const size_t N = pars.N;
    const size_t N2p = 2 * N + 2;
    const size_t N3p = 3 * N + 2;

    TIME(mon.filter_time -= get_wall_time());
    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, inout, tmp.phic);
    fftw_execute_dft_r2c(p_fw, inout + N, tmp.xphic);
    #if PSI_METHOD != PSI_ELLIPTIC
    fftw_execute_dft_r2c(p_fw, inout + N2p, tmp.yphic);
        #if PSI_METHOD == PSI_HYPERBOLIC
        fftw_execute_dft_r2c(p_fw, inout + N3p, tmp.zphic);
        #endif
    #endif
    TIME(mon.fftw_time_exe += get_wall_time());

    apply_filter_fourier(tmp.phic, tmp.xphic, tmp.yphic, tmp.zphic);

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, tmp.phic, inout);
    fftw_execute_dft_c2r(p_bw, tmp.xphic, inout + N);
    #if PSI_METHOD != PSI_ELLIPTIC
    fftw_execute_dft_c2r(p_bw, tmp.yphic, inout + N2p);
        #if PSI_METHOD == PSI_HYPERBOLIC
        fftw_execute_dft_c2r(p_bw, tmp.zphic, inout + N3p);
        #endif
    #endif
    TIME(mon.fftw_time_exe += get_wall_time());
    TIME(mon.filter_time += get_wall_time());
}

/**
 * @brief Applying the filter mask to the complex fields
 *
 * @param[in, out] phi_io The field $$\phi$$ in Fourier space
 * @param[in, out] dphi_io The field $$\dot{\phi}$$ in Fourier space
 * @param[in, out] psi_io The field $$\psi$$ in Fourier space
 * @param[in, out] dpsi_io The field $$\dot{\psi}$$ in Fourier space
 *
 * The `filter` is constructed in `mk_filter_mask()` using the filter window
 * `filter_window(const double x)` in `setup.c`.
 */
static void apply_filter_fourier(fftw_complex *phi_io, fftw_complex *dphi_io,
        fftw_complex *psi_io, fftw_complex *dpsi_io)
{
    const size_t M = pars.M;
    double fil;
    #pragma omp parallel for private(fil)
    for (size_t i = 0; i < M; ++i) {
        fil = filter[i];
        phi_io[i] *= fil;
        dphi_io[i] *= fil;
        #if PSI_METHOD != PSI_ELLIPTIC
        psi_io[i] *= fil;
            #if PSI_METHOD == PSI_HYPERBOLIC
            dpsi_io[i] *= fil;
            #endif
        #endif
    }
}

/**
 * @brief Recompute all desired output quantities and save the current
 * timeslice to buffers
 */
void prepare_and_save_timeslice()
{
    #ifdef OUTPUT_PHI_PS
    evo_flags.compute_pow_spec = 1;
    #endif
    mk_gradient_squared_and_laplacian(field);
    evo_flags.compute_pow_spec = 0;
    mk_rho(field);
    #if PSI_METHOD == PSI_ELLIPTIC && !defined(EVOLVE_WITHOUT_PSI)
    mk_psi(field);
    #endif
    mk_summary();
    #ifdef OUTPUT_CONSTRAINTS
    mk_constraints();
    #endif
    save();
}

/**
 * @brief Center input vector around it's average.
 *
 * @param[in, out] f Any double vector of length @p N
 * @param[in] N The length of the vector @p f
 *
 * The vector @p f is overwritten by f - <f>
 */
static void center(double *f, const size_t N)
{
    double avg = mean(f, N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        f[i] -= avg;
    }
}

/**
 * @brief Save the summaries of the fields (containing the mean, variance,
 * minimum and maximum value at the current time).
 */
void mk_summary()
{
    //TODO[performance]: parallel sections instead of parallel loops here?
    TIME(mon.smry_time -= get_wall_time());
    #ifdef OUTPUT_PHI_SMRY
    mean_var_min_max(field, phi_smry.tmp);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    mean_var_min_max(field + pars.N, dphi_smry.tmp);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    mean_var_min_max(field + 2 * pars.N + 2, psi_smry.tmp);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
        #if PSI_METHOD == PSI_PARABOLIC
        mean_var_min_max(dfield + 2 * pars.N + 2, dpsi_smry.tmp);
        #else
        mean_var_min_max(field + 3 * pars.N + 2, dpsi_smry.tmp);
        #endif
    #endif
    #ifdef OUTPUT_RHO_SMRY
    mean_var_min_max(rho, rho_smry.tmp);
    #endif
    TIME(mon.smry_time += get_wall_time());
}

/**
 * @brief Compute the mean (or average) of a vector
 *
 * @param[in] f Any vector of length @p N
 * @param[in] N The length of the vector @p f
 * @return The mean value of @p f
 */
static double mean(const double *f, const size_t N)
{
    double mean = 0.0;
    #pragma omp parallel for reduction(+: mean)
    for (size_t i = 0; i < N; ++i) {
        mean += f[i];
    }
    return mean / (double)N;
}

/**
 * @brief Compute the summary of a vector, i.e. the mean, variance, minimum and
 * maximum value
 *
 * @param[in] f The input vector
 * @param[out] smry An array of size 4 which is filled with the summary: mean,
 * variance, min, max (in this order)
 *
 * @note The vector is implicitly assumed to have length `pars.N`
 */
static void mean_var_min_max(const double *f, double *smry)
{
    const size_t N = pars.N;
    double mean = 0.0;
    double min_val = f[0];
    double max_val = f[0];
    #pragma omp parallel for reduction(+: mean) reduction(max: max_val) \
        reduction(min: min_val)
    for (size_t i = 0; i < N; ++i) {
        mean += f[i];
        min_val = MIN(min_val, f[i]);
        max_val = MAX(max_val, f[i]);
    }
    smry[0] = mean / (double)N;
    smry[1] = variance(smry[0], f, N);
    smry[2] = min_val;
    smry[3] = max_val;
}

/**
 * @brief Compute the variance of a vector
 *
 * @param[in] mean The mean of the vector @p f
 * @param[in] f Any vector of length @p N
 * @param[in] N The length of the vector @p f
 * @return The variance value of @p f
 */
static double variance(const double mean, const double *f, const size_t N)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    double tmp;
    #pragma omp parallel for private(tmp) reduction(+: sum1, sum2)
    for (size_t i = 0; i < N; ++i) {
        tmp = f[i] - mean;
        sum1 += tmp * tmp;
        sum2 += tmp;
    }
    return (sum1 - sum2 * sum2 / (double)N) / (double)(N - 1);
}

/**
 * @brief Check and print whether a vector contains NaNs __(debugging only)__
 *
 * @param[in] f Any vector of length @p N
 * @param[in] N The length of the vector @p f
 */
static void contains_nan(const double *f, const size_t N)
{
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        if (isnan(f[i])) {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}

/**
 * @brief Check and print whether a vector contains NaNs __(debugging only)__
 *
 * @param[in] f Any vector of length @p N
 * @param[in] N The length of the vector @p f
 */
static void contains_nanc(const complex *f, const size_t N)
{
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        if (isnan(creal(f[i])) || isnan(cimag(f[i]))) {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}
