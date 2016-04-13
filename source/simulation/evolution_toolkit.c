#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"
#include "filehandling.h"
#include "setup.h"

/**
 * @brief A collection of routines to compute the right hand side of the pde as
 * well as the desired output.
 */

struct evolution_flags evo_flags = {.filter = 0, .compute_pow_spec = 0};

// compute right hand side of the pde, i.e. all first order temporal derivatives
// which fields are contained depends on PSI_METHOD
void mk_rhs(const double t, double *f, double *result)
{
    const size_t N = pars.N;
    const size_t N2 = 2 * N;
    const size_t N2p = N2 + 2;
    const size_t N3p = 3 * N + 2;
    const double a = f[N2];
    const double a2 = a * a;

    mk_gradient_squared_and_laplacian(f);
    mk_rho(f);
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

// compute the laplacian and the squared gradient of the input and store them
void mk_gradient_squared_and_laplacian(double *in)
{
    const size_t N = pars.N;
    const size_t M = pars.M;

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw, in, tmp.phic);
        #if PSI_METHOD == PSI_PARABOLIC
        const size_t N2p = 2 * N + 2;
        fftw_execute_dft_r2c(p_fw, in + N2p, tmp.psic);
        #endif
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

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
        #if PSI_METHOD == PSI_PARABOLIC
        tmp.psic[i] *= kvec.sq[i] / N;
        #endif
    }

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw, tmp.xphic, tmp.xphi);
    if (pars.dim > 1) {
        fftw_execute_dft_c2r(p_bw, tmp.yphic, tmp.yphi);
        if (pars.dim > 2) {
            fftw_execute_dft_c2r(p_bw, tmp.zphic, tmp.zphi);
        }
    }
    fftw_execute_dft_c2r(p_bw, tmp.phic, tmp.lap);
    #if PSI_METHOD == PSI_PARABOLIC
    fftw_execute_dft_c2r(p_bw, tmp.psic, tmp.f);
    #endif
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    assemble_gradient_squared();
}

void assemble_gradient_squared()
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

// compute energy density rho & average value
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

// A selection of potentials one can try, make sure to set the corresponding
// potential_prime, the derivative is not computed automatically
inline double potential(const double f)
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

inline double potential_prime(const double f)
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

// solve poisson like equation for scalar perturbation and its derivative
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
    /* fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_c2r(p_bw, tmp.psic, f + N2p); */
    /* fftw_execute_dft_c2r(p_bw, tmp.dpsic, f + N3p); */
    /* #ifdef SHOW_TIMING_INFO */
    /* fftw_time_exe += get_wall_time(); */
    /* poisson_time += get_wall_time(); */
    /* #endif */

    // sophisticated version square first then average plus gradient
    /* #ifdef SHOW_TIMING_INFO */
    /* poisson_time -= get_wall_time(); */
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
    /* fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_r2c(p_fw, tmp.f, tmp.fc); */
    /* fftw_execute_dft_r2c(p_fw, tmp.deltarho, tmp.deltarhoc); */
    /* #ifdef SHOW_TIMING_INFO */
    /* fftw_time_exe += get_wall_time(); */
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
    /* fftw_time_exe -= get_wall_time(); */
    /* #endif */
    /* fftw_execute_dft_c2r(p_bw, tmp.psic, f + N2p); */
    /* fftw_execute_dft_c2r(p_bw, tmp.dpsic, f + N3p); */
    /* #ifdef SHOW_TIMING_INFO */
    /* fftw_time_exe += get_wall_time(); */
    /* poisson_time += get_wall_time(); */
    /* #endif */

    // simpler version
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
    fftw_execute_dft_r2c(p_fw, tmp.deltarho, tmp.deltarhoc);
    fftw_execute_dft_r2c(p_fw, tmp.f, tmp.fc);
    for (size_t i = 1; i < M; ++i) {
        tmp.phic[i] = 0.5 * (tmp.deltarhoc[i] +
                3 * hubble * tmp.fc[i]) / ((kvec.sq[i] / a2 + extra) * N);
    }
    tmp.phic[0] = 0.0;
    fftw_execute_dft_c2r(p_bw, tmp.phic, f + N2p);
    for (size_t i = 0; i < N; ++i) {
        f[N3p + i] = 0.5 * tmp.f[i] - hubble * f[N2p + i];
    }
}

// computes a crude estimation of the power spectrum, more info in main.h
void mk_power_spectrum(const fftw_complex *in, struct output out)
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

// filter the real input field, input gets overwritten with filtered data
void apply_filter_real(double *inout)
{
    const size_t N = pars.N;
    const size_t N2p = 2 * N + 2;
    const size_t N3p = 3 * N + 2;

    #ifdef SHOW_TIMING_INFO
    filter_time -= get_wall_time();
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw, inout, tmp.phic);
    fftw_execute_dft_r2c(p_fw, inout + N, tmp.xphic);
    #if PSI_METHOD != PSI_ELLIPTIC
    fftw_execute_dft_r2c(p_fw, inout + N2p, tmp.yphic);
        #if PSI_METHOD == PSI_HYPERBOLIC
        fftw_execute_dft_r2c(p_fw, inout + N3p, tmp.zphic);
        #endif
    #endif
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    apply_filter_fourier(tmp.phic, tmp.xphic, tmp.yphic, tmp.zphic);

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw, tmp.phic, inout);
    fftw_execute_dft_c2r(p_bw, tmp.xphic, inout + N);
    #if PSI_METHOD != PSI_ELLIPTIC
    fftw_execute_dft_c2r(p_bw, tmp.yphic, inout + N2p);
        #if PSI_METHOD == PSI_HYPERBOLIC
        fftw_execute_dft_c2r(p_bw, tmp.zphic, inout + N3p);
        #endif
    #endif
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    filter_time += get_wall_time();
    #endif
}

// filtering in fourier domain for phi, dphi, psi simultaneously
void apply_filter_fourier(fftw_complex *phi_io, fftw_complex *dphi_io,
        fftw_complex *psi_io, fftw_complex *dpsi_io)
{
    const size_t M = pars.M;
    double fil;
    #pragma omp parallel for private(fil)
    for (size_t i = 0; i < M; ++i) {
        fil = filter[i];
        phi_io[i]  *= fil;
        dphi_io[i] *= fil;
        #if PSI_METHOD != PSI_ELLIPTIC
        psi_io[i]  *= fil;
            #if PSI_METHOD == PSI_HYPERBOLIC
            dpsi_io[i]  *= fil;
            #endif
        #endif
    }
}

// recompute current power spectrum and rho and save current timeslice to buffer
void prepare_and_save_timeslice()
{
    evo_flags.compute_pow_spec = 1;
    mk_gradient_squared_and_laplacian(field);
    evo_flags.compute_pow_spec = 0;
    mk_rho(field);
    #if PSI_METHOD == PSI_ELLIPTIC && !defined(EVOLVE_WITHOUT_PSI)
    mk_psi(field);
    #endif
    mk_summary();
    save();
}

void center(double *f, const size_t N)
{
    double avg = mean(f, N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        f[i] -= avg;
    }
}

void mk_summary()
{
    //TODO[performance]: parallel sections instead of parallel loops here?
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
}

inline double mean(const double *f, const size_t N)
{
    double mean = 0.0;
    #pragma omp parallel for reduction(+: mean)
    for (size_t i = 0; i < N; ++i) {
        mean += f[i];
    }
    return mean / (double)N;
}

void mean_var_min_max(const double *f, double *smry)
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

double variance(const double mean, const double *f, const size_t N)
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

void contains_nan(const double *f, const size_t N)
{
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        if (isnan(f[i])) {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}

void contains_nanc(const complex *f, const size_t N)
{
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        if (isnan(creal(f[i])) || isnan(cimag(f[i]))) {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}
