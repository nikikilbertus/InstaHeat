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
static void mk_stt(const double *f, complex **fsij);
static double potential(const double f);
static double potential_prime(const double f);
static void mk_gw_spectrum(double *f);
#ifdef OUTPUT_CONSTRAINTS
static void mk_constraints(double *f);
#endif
#ifdef OUTPUT_PS
static void mk_power_spectrum(const fftw_complex *in, struct output out);
#endif
#ifdef ENABLE_FFT_FILTER
static void apply_filter_real(double *inout);
static void apply_filter_fourier(fftw_complex *phi_io, fftw_complex *dphi_io,
        fftw_complex *psi_io, fftw_complex *dpsi_io);
#endif
/* static void center(double *f, const size_t N); */
static double mean(const double *f, const size_t N);
static void mean_var_min_max(const double *f, double *smry);
static double variance(const double mean, const double *f, const size_t N);
#ifdef CHECK_FOR_NAN
static void contains_nan(const double *f, const size_t N);
static void contains_nanc(const complex *f, const size_t N);
#endif

struct evolution_flags evo_flags = {.filter = 0,
                                    .output = 0};

/**
 * @brief Compute the right hand side of the pde, i.e. all first order temporal
 * derivatives.
 *
 * @param[in] t The current time
 * @param[in] f All necessary fields bundled in one array
 * @param[out] result The right hand side of the partial differential equation,
 * i.e. the first temporal derivatives of the fields
 *
 * We evolve $$\phi$$, $$\dot{\phi}$$, $$\psi$$, $$\dot{\psi}$$ and $$a$$
 * according to a hyperbolic constraint equation with the integration routine.
 */
void mk_rhs(const double t, double *f, double *result)
{
    mon.calls_rhs += 1;
    const size_t N = pars.N;
    const size_t Next = pars.Next;
    const size_t N2 = 2 * N;
    const size_t N3 = 3 * N;
    const size_t Nh1 = 4 * N;
    const size_t Nh2 = Nh1 + Next;
    const size_t Ndh1 = Nh2 + Next;
    const size_t Ndh2 = Ndh1 + Next;
    const double a = f[pars.Ntot - 1];
    const double a2 = a * a;

    mk_gradient_squared_and_laplacian(f);
    mk_rho(f);
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;

    if (evo_flags.output == 1) {
        // power spectrum of rho here, because need rho first
        #ifdef OUTPUT_RHO_PS
        TIME(mon.fftw_time_exe -= get_wall_time());
        fftw_execute_dft_r2c(p_fw, rho, tmp.phic);
        TIME(mon.fftw_time_exe += get_wall_time());
        mk_power_spectrum(tmp.phic, rho_ps);
        #endif
        #ifdef OUTPUT_CONSTRAINTS
        mk_constraints(f);
        #endif
        //TODO: is this correct use and position?
        mk_gw_spectrum(f);
    }

    double p, dp;
    #pragma omp parallel for private(p, dp)
    for (size_t i = 0; i < N; ++i) {
        p = f[N2 + i];
        dp = f[N3 + i];
        result[i] = f[N + i]; // copy dphi
        result[N + i] = (1.0 + 4.0 * p) * tmp.lap[i] / a2 -
            (h3 - 4.0 * dp) * f[N + i] -
            (1.0 + 2.0 * p) * potential_prime(f[i]); // eq. for ddphi
        result[N2 + i] = dp; // copy dpsi
        result[N3 + i] = 0.5 * pressure[i] + (p - 0.5) * pressure_mean
            - 4.0 * hubble * p; // eq. for ddpsi
    }

    #pragma omp parallel for
    for (size_t i = 0; i < 2 * Next; ++i) {
        result[Nh1 + i] = f[Ndh1 + i]; // copy dhijs
    }

    const size_t len = 6;
    complex **stt = malloc(len * sizeof *stt);
    for (size_t i = 0; i < len; ++i) {
        stt[i] = fftw_malloc(pars.M * sizeof *stt[i]);
    }
    mk_stt(f, stt);

    size_t i1, i2;
    double t1, t2;
    #pragma omp parallel for private(i1, i2, t1, t2)
    for (size_t i = 0; i < pars.M; ++i) {
        i1 = 2 * i;
        i2 = i1 + 1;
        t1 = - (2.0 * pressure_mean + kvec.sq[i] / a);
        t2 = 2.0 / a;
        result[Ndh1 + i1] = t1 * f[Nh1 + i1] - h3 * f[Ndh1 + i1] +
                            t2 * creal(stt[0][i]);
        result[Ndh1 + i2] = t1 * f[Nh1 + i2] - h3 * f[Ndh1 + i2] +
                            t2 * cimag(stt[0][i]);
        result[Ndh2 + i1] = t1 * f[Nh2 + i1] - h3 * f[Ndh2 + i1] +
                            t2 * creal(stt[1][i]);
        result[Ndh2 + i2] = t1 * f[Nh2 + i2] - h3 * f[Ndh2 + i2] +
                            t2 * cimag(stt[1][i]);
    }

    for (size_t i = 0; i < len; ++i) {
        fftw_free(stt[i]);
    }
    free(stt);

    // update da
    result[pars.Ntot - 1] = a * hubble;
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
 * If `OUTPUT_CONSTRAINTS` is defined, additionally, `tmp.f` contains $$\Delta
 * \psi$$, the Lagrangian of $$\psi$$.
 */
void mk_gradient_squared_and_laplacian(double *in)
{
    const size_t N = pars.N;
    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, in, tmp.phic);
    #if defined(OUTPUT_CONSTRAINTS) || defined(OUTPUT_PSI_PS)
    fftw_execute_dft_r2c(p_fw, in + 2 * N, tmp.psic);
    #endif
    TIME(mon.fftw_time_exe += get_wall_time());

    // good place for power spectrum of phi and psi, because fft exists
    if (evo_flags.output == 1) {
        #ifdef OUTPUT_PHI_PS
        mk_power_spectrum(tmp.phic, phi_ps);
        #endif
        #ifdef OUTPUT_PSI_PS
        mk_power_spectrum(tmp.psic, psi_ps);
        #endif
    }

    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        complex pre = tmp.phic[i] * I / N;
        tmp.xphic[i] = pre * kvec.x[i];
        tmp.yphic[i] = pre * kvec.y[i];
        tmp.zphic[i] = pre * kvec.z[i];
        tmp.phic[i] *= kvec.sq[i] / N;
        #ifdef OUTPUT_CONSTRAINTS
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
    #ifdef OUTPUT_CONSTRAINTS
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
    #pragma omp parallel for
    for (size_t i = 0; i < pars.N; ++i) {
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
 * @brief Compute the energy density $$\rho$$ and its average value as well as
 * the pressure and its average value.
 *
 * @param[in] f An array containing the fields
 *
 * Everything is stored in global variables for reuse in other functions. When
 * the function returns:
 * `rho` contains the energy density $$\rho$$
 * `rho_mean` contains the average energy density $$< \rho >$$
 * `pressure` contains the pressure $$p$$
 * `pressure_mean` contains the average pressure $$< p >$$
 * All the above values persist until the next call of `mk_rhs(const double t,
 * double *f, double *result)`
 */
void mk_rho(const double *f)
{
    const size_t N = pars.N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    rho_mean = 0.0;
    pressure_mean = 0.0;
    #pragma omp parallel for reduction(+: rho_mean, pressure_mean)
    for (size_t i = 0; i < N; ++i) {
        double t1 = (0.5 - f[2 * N + i]) * f[N + i] * f[N + i];
        double t2 = (0.5 + f[2 * N + i]) * tmp.grad[i] / a2;
        rho[i] = t1 + t2 + potential(f[i]);
        pressure[i] = t1 - t2 / 3.0 - potential(f[i]);
        pressure_mean += pressure[i];
        rho_mean += rho[i];
    }
    rho_mean /= N;
    pressure_mean /= N;
}

/**
 * @brief Compute the traceless transverse source terms of the tensor
 * perturbations needed in gravitational wave extraction.
 *
 * @param[in] f An array containing the fields
 * @param[in, out] fsij The traceless transverse part of the source terms in
 * the equation of motion of the tensor perturbations. We only need the first
 * two fields, because there are only two degrees of freedom.
 *
 * @see TODO[link thesis and papers]
 */
static void mk_stt(const double *f, complex **fsij)
{
    TIME(mon.stt_time -= get_wall_time());
    const double atmp = 3.0 * f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const size_t len = 6;

    double **sij = malloc(len * sizeof *sij);
    for (size_t i = 0; i < len; ++i) {
        sij[i] = fftw_malloc(pars.N * sizeof *sij[i]);
    }
    #pragma omp parallel for
    for (size_t i = 0; i < pars.N; ++i) {
        // TODO: do i include metric here?
        double gphi = - tmp.grad[i] * (2.0 + 4.0 * f[2 * pars.N + i]) / atmp;
        sij[0][i] = tmp.xphi[i] * tmp.xphi[i] + gphi;
        sij[1][i] = tmp.xphi[i] * tmp.yphi[i];
        sij[2][i] = tmp.xphi[i] * tmp.zphi[i];
        sij[3][i] = tmp.yphi[i] * tmp.yphi[i] + gphi;
        sij[4][i] = tmp.yphi[i] * tmp.zphi[i];
        sij[5][i] = tmp.zphi[i] * tmp.zphi[i] + gphi;
    }
    TIME(mon.fftw_time_exe -= get_wall_time());
    for (size_t i = 0; i < len; ++i) {
        fftw_execute_dft_r2c(p_fw, sij[i], fsij[i]);
        fftw_free(sij[i]);
    }
    TIME(mon.fftw_time_exe += get_wall_time());
    free(sij);

    #pragma omp parallel for
    for (size_t i = 1; i < pars.M; ++i) {
        double kx = kvec.xf[i];
        double ky = kvec.yf[i];
        double kz = kvec.zf[i];
        double fx = kx / kvec.sq[i];
        double fy = ky / kvec.sq[i];
        double fz = kz / kvec.sq[i];
        complex t1 = kx * fx * fsij[0][i] + 2.0 * kx * fy * fsij[1][i] +
             2.0 * kx * fz * fsij[2][i] + ky * fy * fsij[3][i] +
             2.0 * ky * fz * fsij[4][i] + kz * fz * fsij[5][i];
        complex t2 = fsij[0][i] + fsij[3][i] + fsij[5][i];
        complex s1 = t1 + t2;
        complex s2 = t1 - t2;

        // use s11 and s12
        if (fabs(kz) > DBL_EPSILON) {
            complex k1 = kx * fsij[0][i] + ky * fsij[1][i] + kz * fsij[2][i];
            complex k2 = kx * fsij[1][i] + ky * fsij[3][i] + kz * fsij[4][i];
            fsij[0][i] = fsij[0][i] - 2.0 * fx * k1 +
                0.5 * (fx * kx * s1 + s2);
            fsij[1][i] = fsij[1][i] - fx * k2 - fy * k1 +
                0.5 * fx * ky * s1;
        // use s11 and s13
        } else if (fabs(ky) > DBL_EPSILON) {
            complex k1 = kx * fsij[0][i] + ky * fsij[1][i] + kz * fsij[2][i];
            complex k3 = kx * fsij[2][i] + ky * fsij[4][i] + kz * fsij[5][i];
            fsij[0][i] = fsij[0][i] - 2.0 * fx * k1 +
                0.5 * (fx * kx * s1 + s2);
            fsij[1][i] = fsij[2][i] - fx * k3 - fz * k1 +
                0.5 * fx * kz * s1;
        // use s22 and s23
        } else {
            complex k2 = kx * fsij[1][i] + ky * fsij[3][i] + kz * fsij[4][i];
            complex k3 = kx * fsij[2][i] + ky * fsij[4][i] + kz * fsij[5][i];
            fsij[0][i] = fsij[3][i] - 2.0 * fy * k2 +
                0.5 * (fy * ky * s1 + s2);
            fsij[1][i] = fsij[4][i] - fy * k3 - fz * k2 +
                0.5 * fy * kz * s1;
        }
    }
    fsij[0][0] = 0.0;
    fsij[1][0] = 0.0;
    TIME(mon.stt_time += get_wall_time());
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
 * @brief Compute the power spectrum of the gravitational waves.
 *
 * Computes the power spectrum of the gravitational waves associated with the
 * tensor perturbations in @p f according to the prescription in TODO[link to
 * paper].
 *
 * @param[in] f An array with the current fields values
 */
static void mk_gw_spectrum(double *f)
{
    const size_t Ndh1 = 4 * pars.N + 2 * pars.Next;
    const size_t Ndh2 = Ndh1 + pars.Next;
    const double k2_max = pars.x.k2 * (pars.x.N/2) * (pars.x.N/2) +
                          pars.y.k2 * (pars.y.N/2) * (pars.y.N/2) +
                          pars.z.k2 * (pars.z.N/2) * (pars.z.N/2);

    #pragma omp parallel for
    for (size_t i = 0; i < gw.dim; ++i) {
        gw.tmp[i] = 0.0;
    }

    // cannot parallelize due to binning
    for (size_t i = 1; i < pars.M; ++i) {
        double kx = kvec.xf[i];
        double ky = kvec.yf[i];
        double kz = kvec.zf[i];
        double k2 = kvec.sq[i];
        double dh1r = f[Ndh1 + 2 * i];
        double dh2r = f[Ndh2 + 2 * i];
        double dh1i = f[Ndh1 + 2 * i + 1];
        double dh2i = f[Ndh2 + 2 * i + 1];
        double pow;
        // use h11 and h12
        if (fabs(kz) > DBL_EPSILON) {
            pow = 2.0 * k2 / (kz * kz * (ky * ky + kz * kz)) *
                ((kx * kx + kz * kz) * (dh1r * dh1r + dh1i * dh1i) +
                2.0 * kx * ky * (dh1r * dh2r + dh1i * dh2i) +
                (kx * kx + kz * kz) * (dh2r * dh2r + dh2i * dh2i));
        // use h11 and h13
        } else if (fabs(ky) > DBL_EPSILON) {
            pow = 2.0 * (kx * kx + ky * ky) / (ky * ky * ky * ky) *
                ((kx * kx + ky * ky) * (dh1r * dh1r + dh1i * dh1i) +
                ky * ky * (dh2r * dh2r + dh2i * dh2i));
        // use h22 and h23
        } else {
            pow = 2.0 * (dh1r * dh1r + dh1i * dh1i +
                         dh2r * dh2r + dh2i * dh2i);
        }
        if (fabs(kvec.z[i]) > DBL_EPSILON) {
            pow *= 2.0;
        }
        size_t idx = (int)trunc(gw.dim * sqrt(k2 / k2_max) - 1.0e-14);
        gw.tmp[idx] += pow / pars.N;
    }
}

#ifdef OUTPUT_CONSTRAINTS
/**
 * @brief Monitor how well the Hamiltonian and momentum constraint are fulfilled.
 *
 * Computes the Hamiltonian and momentum constraint in a form such they should
 * give 0, if fulfilled exactly. Of these combinations save the l2 Norm as well
 * as the maximum asbolute value to `cstr.tmp` for output.
 *
 * @param[in] f An array containing the fields
 * @note So far only the Hamiltonian constraint is implemented! (4/14/2016)
 */
static void mk_constraints(double *f)
{
    TIME(mon.cstr_time -= get_wall_time());
    const size_t N = pars.N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;
    const double phi_mean = mean(f, N);
    const double dphi_mean = mean(f + N, N);
    double ham_l2 = 0.0, ham_max = 0.0;
    double mom_l2 = 0.0, mom_max = 0.0;

    #pragma omp parallel for reduction(max: ham_max, mom_max) \
                             reduction(+: ham_l2, mom_l2)
    for (size_t i = 0; i < N; ++i) {
        double tmp1 = hubble * f[2 * N + i] + f[3 * N + i];
        double ham = tmp.f[i] / a2 - h3 * tmp1 - 0.5 * (rho[i] - rho_mean);
        double mom = tmp1 - 0.5 * dphi_mean * (f[i] - phi_mean);
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
#endif

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
    TIME(mon.poisson_time -= get_wall_time());
    const size_t N = pars.N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
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

    #pragma omp parallel for
    for (size_t i = 1; i < pars.M; ++i) {
        tmp.phic[i] = 0.5 * (tmp.deltarhoc[i] +
            3.0 * hubble * tmp.fc[i]) / ((kvec.sq[i] / a2 + extra) * N);
    }
    tmp.phic[0] = 0.0;

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, tmp.phic, f + 2 * N);
    TIME(mon.fftw_time_exe += get_wall_time());

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        f[3 * N + i] = 0.5 * tmp.f[i] - hubble * f[2 * N + i];
    }
    TIME(mon.poisson_time += get_wall_time());
}

#ifdef OUTPUT_PS
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
    const double k2_max = pars.x.k2 * (pars.x.N/2) * (pars.x.N/2) +
                          pars.y.k2 * (pars.y.N/2) * (pars.y.N/2) +
                          pars.z.k2 * (pars.z.N/2) * (pars.z.N/2);

    #pragma omp parallel for
    for (size_t i = 0; i < out.dim; ++i) {
        out.tmp[i] = 0.0;
    }

    // starting with 1 to explicitly exclude constant average
    for (size_t i = 1; i < pars.M; ++i) {
        double pow2_tmp = in[i] * conj(in[i]);
        if (fabs(kvec.z[i]) > DBL_EPSILON) {
            pow2_tmp *= 2.0;
        }
        size_t idx = (int)trunc(out.dim * sqrt(kvec.sq[i] / k2_max) - 1.0e-14);
        out.tmp[idx] += pow2_tmp / pars.N;
    }
}
#endif

#ifdef ENABLE_FFT_FILTER
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
    TIME(mon.filter_time -= get_wall_time());
    const size_t N = pars.N;
    const size_t N2 = 2 * N;
    const size_t N3 = 3 * N;

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, inout, tmp.phic);
    fftw_execute_dft_r2c(p_fw, inout + N, tmp.xphic);
    fftw_execute_dft_r2c(p_fw, inout + N2, tmp.yphic);
    fftw_execute_dft_r2c(p_fw, inout + N3, tmp.zphic);
    TIME(mon.fftw_time_exe += get_wall_time());

    apply_filter_fourier(tmp.phic, tmp.xphic, tmp.yphic, tmp.zphic);

    TIME(mon.fftw_time_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, tmp.phic, inout);
    fftw_execute_dft_c2r(p_bw, tmp.xphic, inout + N);
    fftw_execute_dft_c2r(p_bw, tmp.yphic, inout + N2);
    fftw_execute_dft_c2r(p_bw, tmp.zphic, inout + N3);
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
        psi_io[i] *= fil;
        dpsi_io[i] *= fil;
    }
}
#endif

/**
 * @brief Recompute all desired output quantities and save the current
 * timeslice to buffers
 */
void prepare_and_save_timeslice()
{
    evo_flags.output = 1;
    mk_rhs(pars.t.t, field, dfield);
    evo_flags.output = 0;
    mk_summary();
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
/* static void center(double *f, const size_t N) */
/* { */
/*     double avg = mean(f, N); */
/*     #pragma omp parallel for */
/*     for (size_t i = 0; i < N; ++i) { */
/*         f[i] -= avg; */
/*     } */
/* } */

/**
 * @brief Save the summaries of the fields (containing the mean, variance,
 * minimum and maximum value at the current time).
 */
void mk_summary()
{
    TIME(mon.smry_time -= get_wall_time());
    #ifdef OUTPUT_PHI_SMRY
    mean_var_min_max(field, phi_smry.tmp);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    mean_var_min_max(field + pars.N, dphi_smry.tmp);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    mean_var_min_max(field + 2 * pars.N, psi_smry.tmp);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    mean_var_min_max(field + 3 * pars.N, dpsi_smry.tmp);
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

#ifdef CHECK_FOR_NAN
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
#endif
