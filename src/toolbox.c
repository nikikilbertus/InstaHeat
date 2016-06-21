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

static void mk_ffts_and_filter(double *f);
static void assemble_gradient_squared();
static void output_all(double *f);
static void update_phi_psi(double *f, double *result);
#ifdef ENABLE_GW
static void update_h(double *f, double *result);
static void mk_gw_sources_tt(const double *f, complex **s);
static void mk_gw_sources(const double *f, complex **s);
static void mk_gw_spectrum(double *f);
#endif
static double potential(const double f);
static double potential_prime(const double f);
#ifdef OUTPUT_CONSTRAINTS
static void mk_constraints(double *f);
#endif
#ifdef OUTPUT_PS
static void mk_power_spectrum(const fftw_complex *in, struct output out);
#endif
#ifdef ENABLE_FFT_FILTER
static void apply_filter(double *f);
static void capply_filter(complex *in, double *out);
#endif
static void mk_summary();
static void mean_var_min_max(const double *f, double *smry, const size_t N);
#if defined(OUTPUT_H1_SMRY) || defined(OUTPUT_H2_SMRY)
static void fmean_var_min_max(const double *f, double *smry);
#endif
static double mean(const double *f, const size_t N);
static double variance(const double mean, const double *f, const size_t N);
static void real_to_complex(const double *in, complex *out);
static void complex_to_real(const complex *in, double *out);
static void fft(double *in, complex *out);
static void ifft(complex *in, double *out);
#ifdef CHECK_FOR_NAN
static void contains_nan(const double *f, const size_t N);
static void contains_nanc(const complex *f, const size_t N);
#endif

/**
 * @brief THe only instance of the `evolution_flags` structure.
 */
struct evolution_flags evo_flags = {.filter = 0, .output = 0};

/**
 * @brief Compute the right hand side of the pde, i.e. all first order temporal
 * derivatives of the fields.
 *
 * @param[in] t The current time.
 * @param[in] f The fields.
 * @param[out] result The right hand side of the partial differential equation,
 * i.e. the first temporal derivatives of the fields in @p f.
 *
 * This function has a whole lot of side effects. It will also trigger the
 * computation of all desired outputs. Depending on the simulation parameters
 * different functions will be called cascading from `mk_rhs(const double t,
 * double *f, double *result)`.
 *
 * @see TODO[link] for more on the evolved fields and the equations of motion.
 */
void mk_rhs(const double t, double *f, double *result)
{
    mon.calls_rhs += 1;
    mk_gradient_squared_and_laplacian(f);
    mk_rho_and_p(f);
    if (evo_flags.output == 1) {
        output_all(f);
    }
    update_phi_psi(f, result);
    #ifdef ENABLE_GW
    update_h(f, result);
    #endif
    result[pars.Ntot - 1] = f[pars.Ntot - 1] * sqrt(rho_mean / 3.0); //update da
}

/**
 * @brief Recompute all desired output quantities and copy the current
 * timeslice to buffers.
 *
 * @note This is only called once or twice throughout a whole simulation run.
 */
void prepare_and_save_timeslice()
{
    evo_flags.output = 1;
    mk_rhs(pars.t.t, field, dfield);
    evo_flags.output = 0;
    save();
}

/**
 * @brief Compute the Laplacian and the sqaured gradient. More computations
 * might be performed depending on flags in the `evolution_flags` structure
 * `evo_flags`.
 *
 * @param[in, out] in The fields, in particular \f$\phi\f$ and \f$\psi\f$. All
 * fields in @p in might be overwritten by their filtered versions.
 *
 * The Laplacian and the squared gradient as well as the partial derivatives of
 * \f$\phi\f$ are stored in global variables for reuse in other functions. When
 * `mk_gradient_squared_and_laplacian(double *in)` returns:
 * `tmp.xphi` contains \f$\partial_x \phi\f$.
 * `tmp.yphi` contains \f$\partial_y \phi\f$.
 * `tmp.zphi` contains \f$\partial_z \phi\f$.
 * `tmp.lap` contains the Laplacian \f$\sum_{i=1}^3 \partial_i^2 \phi\f$
 * `tmp.grad` contains the _squared_ gradient \f$(\nabla \phi)^2\f$
 * All the above values persist until the next call of `mk_rhs(const double t,
 * double *f, double *result)`.
 * If required, the Fourier transform of \f$\psi\f$ is computed and `tmp.f`
 * contains \f$\Delta \psi\f$, the Lagrangian of \f$\psi\f$.
 * If required, the fields \f$\phi\f$, \f$\psi\f$, \f$\dot{\phi}\f$,
 * \f$\dot{\psi}\f$ will be overwritten with their filtered versions.
 */
void mk_gradient_squared_and_laplacian(double *in)
{
    mk_ffts_and_filter(in);
    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        complex pre = tmp.phic[i] / pars.N;
        tmp.xphic[i] = pre * I * kvec.x[i];
        tmp.yphic[i] = pre * I * kvec.y[i];
        tmp.zphic[i] = pre * I * kvec.z[i];
        tmp.fc[i] = - pre * kvec.sq[i];
    }

    ifft(tmp.xphic, tmp.xphi);
    if (pars.dim > 1) {
        ifft(tmp.yphic, tmp.yphi);
        if (pars.dim > 2) {
            ifft(tmp.zphic, tmp.zphi);
        }
    }
    ifft(tmp.fc, tmp.lap);
    assemble_gradient_squared();
}

/**
 * @brief Computes the DFTs which are necessary for the right hand side or the
 * desired output. This will also trigger the filtering if required.
 *
 * @param[in, out] f The fields. They might be overwritten by their filtered
 * versions.
 *
 * Additionally, if required, the Laplacian of \f$\psi\f$ is computed here.
 */
static void mk_ffts_and_filter(double *f)
{
    const size_t N = pars.N;
    fft(f, tmp.phic);
    #if defined(ENABLE_FFT_FILTER)
    if (evo_flags.filter == 1) { // output == 1 => filter == 1
        fft(f + 2 * N, tmp.psic);
        capply_filter(tmp.phic, f); // fft already available
        capply_filter(tmp.psic, f + 2 * N); // fft already available
        apply_filter(f + pars.N); // fft not available
        apply_filter(f + 3 * pars.N); // fft not available
    }
    #elif defined(OUTPUT_CONSTRAINTS) || defined(OUTPUT_PSI_PS)
    if (evo_flags.output == 1) {
        fft(f + 2 * N, tmp.psic);
    }
    #endif
    #ifdef OUTPUT_CONSTRAINTS
    if (evo_flags.output == 1) {
        #pragma omp parallel for
        for (size_t i = 0; i < pars.M; ++i) {
            tmp.fc[i] = - tmp.psic[i] * kvec.sq[i] / N;
        }
        ifft(tmp.fc, tmp.f);
    }
    #endif
}

/**
 * @brief Constructs the squared gradient from the (up to) three partial
 * spatial derivatives (depending on the number of spatial dimensions).
 *
 * The partial derivatives of \f$\phi\f$ computed in
 * `mk_gradient_squared_and_laplacian(double *in)` are used to compute the
 * squared gradient of \f$\phi\f$, which is  stored in `tmp.grad`.
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
 * @brief Compute the energy density \f$\rho\f$ and its spatial average as well
 * as the pressure \f$p\f$ and its spatial average.
 *
 * @param[in] f The fields.
 *
 * This function will populate the following global variables:
 * `rho` contains the energy density \f$\rho\f$.
 * `rho_mean` contains the average energy density \f$\langle \rho \rangle\f$.
 * `pressure` contains the pressure \f$p\f$.
 * `pressure_mean` contains the average pressure \f$\langle p \rangle\f$.
 *
 * @note All the above values persist until the next call of `mk_rhs(const
 * double t, double *f, double *result)`.
 */
void mk_rho_and_p(const double *f)
{
    const size_t N = pars.N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    rho_mean = 0.0;
    pressure_mean = 0.0;
    #pragma omp parallel for reduction(+: rho_mean, pressure_mean)
    for (size_t i = 0; i < N; ++i) {
        double dph = f[N + i], ps = f[2 * N + i];
        double t1 = (0.5 - ps) * dph * dph;
        double t2 = (0.5 + ps) * tmp.grad[i] / a2;
        rho[i] = t1 + t2 + potential(f[i]);
        pressure[i] = t1 - t2 / 3.0 - potential(f[i]);
        pressure_mean += pressure[i];
        rho_mean += rho[i];
    }
    rho_mean /= N;
    pressure_mean /= N;
}

/**
 * @brief Checks which outputs are enabled on the current time slice, computes
 * them and fills the corresponding buffers.
 *
 * @param[in] f The fields.
 */
static void output_all(double *f)
{
    #ifdef OUTPUT_PHI_PS
    mk_power_spectrum(tmp.phic, phi_ps);
    #endif
    #ifdef OUTPUT_PSI_PS
    mk_power_spectrum(tmp.psic, psi_ps);
    #endif
    #ifdef OUTPUT_RHO_PS
    fft(rho, tmp.fc);
    mk_power_spectrum(tmp.fc, rho_ps);
    #endif
    #ifdef OUTPUT_CONSTRAINTS
    mk_constraints(f);
    #endif
    mk_summary();
    #ifdef ENABLE_GW
    mk_gw_spectrum(f);
    #endif
}

/**
 * @brief Computes and updates part of the right hand side of the pde, i.e. the
 * first order temporal derivatives of \f$\phi\f$, \f$\dot{phi}\f$, \f$\psi\f$,
 * \f$\dot{\psi}\f$.
 *
 * @param[in] f The fields.
 * @param[out] result The right hand side of the pde for the given fields.
 *
 * @see TODO[link] for more on the evolved fields and their equations of motion.
 */
static void update_phi_psi(double *f, double *result)
{
    const size_t N = pars.N, N2 = 2 * N, N3 = 3 * N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        double dph = f[N + i], ps = f[N2 + i], dps = f[N3 + i];
        result[i] = dph; // update dphi
        result[N + i] = (1.0 + 4.0 * ps) * tmp.lap[i] / a2 -
            (h3 - 4.0 * dps) * dph -
            (1.0 + 2.0 * ps) * potential_prime(f[i]); // update ddphi
        result[N2 + i] = dps; // update dpsi
        result[N3 + i] = 0.5 * pressure[i] + (ps - 0.5) * pressure_mean
            - 4.0 * hubble * dps; // update ddpsi
    }
}

#ifdef ENABLE_GW
/**
 * @brief Computes and updates part of the right hand side of the pde, i.e. the
 * first order temporal derivatives of the tensor perturbations \f$h\f$.
 *
 * @param[in] f The fields.
 * @param[out] result The right hand side of the pde for the given fields.
 *
 * @see TODO[link] for more on the evolved fields and their equations of motion.
 */
static void update_h(double *f, double *result)
{
    const size_t N = pars.N, Next = pars.Next, Nh1 = 4 * N, Nh2 = Nh1 + Next;
    const size_t Ndh1 = Nh2 + Next, Ndh2 = Ndh1 + Next;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
    const double h3 = 3.0 * hubble;

    #pragma omp parallel for
    for (size_t i = 0; i < 2 * Next; ++i) {
        result[Nh1 + i] = f[Ndh1 + i]; // update dhijs
    }
    const size_t len = 6;
    complex **stt = malloc(len * sizeof *stt);
    for (size_t i = 0; i < len; ++i) {
        stt[i] = fftw_malloc(pars.M * sizeof *stt[i]);
    }
    mk_gw_sources_tt(f, stt);

    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        size_t i1 = 2 * i, i2 = i1 + 1;
        double tmp = - (2.0 * pressure_mean + kvec.sq[i] / a2);
        result[Ndh1 + i1] = tmp * f[Nh1 + i1] - h3 * f[Ndh1 + i1] +
                            2.0 * creal(stt[0][i]) / a2;
        result[Ndh1 + i2] = tmp * f[Nh1 + i2] - h3 * f[Ndh1 + i2] +
                            2.0 * cimag(stt[0][i]) / a2;
        result[Ndh2 + i1] = tmp * f[Nh2 + i1] - h3 * f[Ndh2 + i1] +
                            2.0 * creal(stt[1][i]) / a2;
        result[Ndh2 + i2] = tmp * f[Nh2 + i2] - h3 * f[Ndh2 + i2] +
                            2.0 * cimag(stt[1][i]) / a2;
    }
    for (size_t i = 0; i < len; ++i) {
        fftw_free(stt[i]);
    }
    free(stt);
}

/**
 * @brief Compute the traceless transverse source terms for the equation of
 * motion of the tensor perturbations \f$h\f$.
 *
 * @param[in] f The fields.
 * @param[out] s The traceless transverse part of the source terms in
 * the equation of motion of the tensor perturbations.
 *
 * @note From the nine components of \f$h_{ij}\f$ only two are independent once
 * we impose symmetry, transversality and tracelessness. Thus we only keep two
 * of the components.
 *
 * @see The traceless transverse source term is the right hand side of equation
 * TODO[link] in the thesis.
 */
static void mk_gw_sources_tt(const double *f, complex **s)
{
    TIME(mon.gw_sources -= get_wall_time());
    mk_gw_sources(f, s);

    #pragma omp parallel for
    for (size_t i = 1; i < pars.M; ++i) {
        double kx = kvec.xf[i], ky = kvec.yf[i], kz = kvec.zf[i];
        double ksq = kvec.sq[i];
        double fx = kx / ksq, fy = ky / ksq, fz = kz / ksq;
        complex t1 = kx * fx * s[0][i] + ky * fy * s[3][i] +
                     kz * fz * s[5][i] + 2.0 * (kx * fy * s[1][i] +
                     kx * fz * s[2][i] + ky * fz * s[4][i]);
        complex t2 = s[0][i] + s[3][i] + s[5][i];
        complex s1 = t1 + t2, s2 = t1 - t2;

        if (fabs(kz) > DBL_EPSILON) { // use s11 and s12
            complex k1 = kx * s[0][i] + ky * s[1][i] + kz * s[2][i];
            complex k2 = kx * s[1][i] + ky * s[3][i] + kz * s[4][i];
            s[0][i] = s[0][i] - 2.0 * fx * k1 + 0.5 * (fx * kx * s1 + s2);
            s[1][i] = s[1][i] - fx * k2 - fy * k1 + 0.5 * fx * ky * s1;
        } else if (fabs(ky) > DBL_EPSILON) { // use s11 and s13
            complex k1 = kx * s[0][i] + ky * s[1][i] + kz * s[2][i];
            complex k3 = kx * s[2][i] + ky * s[4][i] + kz * s[5][i];
            s[0][i] = s[0][i] - 2.0 * fx * k1 + 0.5 * (fx * kx * s1 + s2);
            s[1][i] = s[2][i] - fx * k3 - fz * k1 + 0.5 * fx * kz * s1;
        } else { // use s22 and s23
            complex k2 = kx * s[1][i] + ky * s[3][i] + kz * s[4][i];
            complex k3 = kx * s[2][i] + ky * s[4][i] + kz * s[5][i];
            s[0][i] = s[3][i] - 2.0 * fy * k2 + 0.5 * (fy * ky * s1 + s2);
            s[1][i] = s[4][i] - fy * k3 - fz * k2 + 0.5 * fy * kz * s1;
        }
    }
    s[0][0] = 0.0;
    s[1][0] = 0.0;
    TIME(mon.gw_sources += get_wall_time());
}

/**
 * @brief Compute the source term, i.e. the right hand side, for the equation of
 * motion of the tensor perturbations \f$h_{ij}\f$.
 *
 * @param[in] f The fields.
 * @param[out] s Holds 6 arrays for the 6 components of the sources \f$S_{ij}\f$
 * in Fourier space after imposing symmetry.
 *
 * @see The source term is the right hand side of equation TODO[link] in the
 * thesis.
 */
static void mk_gw_sources(const double *f, complex **s)
{
    const size_t len = 6;
    double *t1[] = {tmp.xphi, tmp.xphi, tmp.xphi, tmp.yphi, tmp.yphi, tmp.zphi};
    double *t2[] = {tmp.xphi, tmp.yphi, tmp.zphi, tmp.yphi, tmp.zphi, tmp.zphi};
    double *s_tmp = fftw_malloc(pars.N * sizeof *s_tmp);
    for (size_t j = 0; j < len; ++j) {
        #pragma omp parallel for
        for (size_t i = 0; i < pars.N; ++i) {
            s_tmp[i] = t1[j][i] * t2[j][i];
        }
        fft(s_tmp, s[j]);
    }
    fftw_free(s_tmp);
}

/**
 * @brief Computes the power spectrum of the gravitational waves.
 *
 * @param[in] f The fields.
 *
 * Computes the power spectrum of the gravitational waves associated with the
 * tensor perturbations in the fields @p f according to the prescription in
 * TODO[link].
 */
static void mk_gw_spectrum(double *f)
{
    const size_t Ndh1 = 4 * pars.N + 2 * pars.Next;
    const size_t Ndh2 = Ndh1 + pars.Next;
    const double lx= pars.x.b - pars.x.a;
    const double ly= pars.y.b - pars.y.a;
    const double lz= pars.z.b - pars.z.a;
    const double hubble = sqrt(rho_mean / 3.0);
    // ratio dof of today to mater-radiation-equality to the 1/3
    const double rat = pow(0.01, 1.0/3.0);
    const double fac = PI * rat / (3.0 * hubble * hubble * lx * ly * lz);
    #pragma omp parallel for
    for (size_t i = 0; i < gw.dim; ++i) {
        gw.tmp[i] = 0.0;
    }
    // start with 1 to exclude constant offset, no omp!
    for (size_t i = 1; i < pars.M; ++i) {
        double kx = kvec.xf[i], ky = kvec.yf[i], kz = kvec.zf[i];
        double kx2 = kx * kx, ky2 = ky * ky, kz2 = kz * kz;
        double k2 = kvec.sq[i], k = sqrt(k2);
        double dh1r = f[Ndh1 + 2 * i], dh2r = f[Ndh2 + 2 * i];
        double dh1i = f[Ndh1 + 2 * i + 1], dh2i = f[Ndh2 + 2 * i + 1];
        double pow;
        // use h11 and h12
        if (fabs(kz) > DBL_EPSILON) {
            pow = 2.0 * k2 / (kz2 * (ky2 + kz2)) *
                ((kx2 + kz2) * (dh1r * dh1r + dh1i * dh1i) +
                2.0 * kx * ky * (dh1r * dh2r + dh1i * dh2i) +
                (ky2 + kz2) * (dh2r * dh2r + dh2i * dh2i));
        // use h11 and h13
        } else if (fabs(ky) > DBL_EPSILON) {
            pow = 2.0 * (kx2 + ky2) / (ky2 * ky2) *
                ((kx2 + ky2) * (dh1r * dh1r + dh1i * dh1i) +
                ky2 * (dh2r * dh2r + dh2i * dh2i));
        // use h22 and h23
        } else {
            pow = 2.0 * (dh1r * dh1r + dh1i * dh1i +
                         dh2r * dh2r + dh2i * dh2i);
        }
        if (fabs(kvec.z[i]) > DBL_EPSILON) {
            pow *= 2.0;
        }
        size_t idx = (int)trunc(gw.dim * k / kvec.k_max - 1.0e-14);
        gw.tmp[idx] += fac * k2 * k * pow;
    }
    #pragma omp parallel for
    for (size_t i = 0; i < gw.dim; ++i) {
        gw.tmp[i] /= pars.N;
    }
}
#endif

/**
 * @brief The potential \f$V\f$ of the scalar inflaton field \f$\phi\f$.
 *
 * @param[in] f The inflaton field value where to evaluate the potential.
 * @return The value of the potential at the given input @p f.
 */
static double potential(const double f)
{
    // standard phi squared potential
    return MASS * MASS * f * f / 2.0;
}

/**
 * @brief The derivative of the potential \f$V'\f$ of the scalar inflaton field
 * \f$\phi\f$.
 *
 * @param[in] f The value at which to evaluate the derivative of the potential.
 * @return The value of the derivative of the potential at given input @p f.
 */
static double potential_prime(const double f)
{
    // standard phi squared potential
    return MASS * MASS * f;
}

#ifdef OUTPUT_CONSTRAINTS
/**
 * @brief Computes how well the Hamiltonian and momentum constraint are
 * fulfilled for monitoring.
 *
 * @param[in] f The fields.
 *
 * Computes the Hamiltonian and momentum constraint as a sum of terms that
 * should cancel each other if fulfilled exactly. Of these combinations the
 * \f$\ell_2\f$ norm as well as the sup-norm, i.e. the maximum asbolute value,
 * are copied to the corresponding output buffer.
 */
static void mk_constraints(double *f)
{
    TIME(mon.cstr -= get_wall_time());
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
    cstr.tmp[0] = ham_l2 / N;
    cstr.tmp[1] = ham_max / N;
    cstr.tmp[2] = mom_l2 / N;
    cstr.tmp[3] = mom_max / N;
    TIME(mon.cstr += get_wall_time());
}
#endif

/**
 * @brief Computes \f$\psi\f$ and \f$\dot{\psi}\f$ from \f$\phi\f$ and
 * \f$\dot{\phi}\f$ on the initial timeslice.
 *
 * @param[in, out] f The fields. Expects \f$\phi\f$, \f$\dot{\phi}\f$ and \f$a\f$
 * to be given in @p f. Fills in \f$\psi\f$ and \f$\dot{\psi}\f$ in @p f.
 *
 * We use an elliptic equation from the Hamiltonian constraint combined with
 * the momentum contraint to compute \f$\psi\f$ and \f$\dot{\psi}\f$ from given
 * \f$\phi\f$, \f$\dot{\phi}\f$ and \f$a\f$. TODO[link]
 */
void mk_psi(double *f)
{
    TIME(mon.elliptic -= get_wall_time());
    const size_t N = pars.N;
    const double a2 = f[pars.Ntot - 1] * f[pars.Ntot - 1];
    const double hubble = sqrt(rho_mean / 3.0);
    const double phi_mean = mean(f, N);
    const double dphi_mean = mean(f + N, N);
    double extra1 = 0.0, extra2 = 0.0;

    #pragma omp parallel for reduction(+: extra1, extra2)
    for (size_t i = 0; i < N; ++i) {
        tmp.deltarho[i] = rho[i] - rho_mean;
        tmp.f[i] = dphi_mean * (f[i] - phi_mean);
        extra1 += f[N + i] * f[N + i];
        extra2 += tmp.grad[i];
    }
    const double extra = 0.5 * (extra1 - extra2 / a2) / N;

    fft(tmp.deltarho, tmp.deltarhoc);
    fft(tmp.f, tmp.fc);
    #pragma omp parallel for
    for (size_t i = 1; i < pars.M; ++i) {
        tmp.phic[i] = 0.5 * (tmp.deltarhoc[i] +
            3.0 * hubble * tmp.fc[i]) / ((- kvec.sq[i] / a2 + extra) * N);
    }
    tmp.phic[0] = 0.0;

    ifft(tmp.phic, f + 2 * N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        f[3 * N + i] = 0.5 * tmp.f[i] - hubble * f[2 * N + i];
    }
    TIME(mon.elliptic += get_wall_time());
}

#ifdef OUTPUT_PS
/**
 * @brief Computes the power spectrum of a field.
 *
 * @param[in] in The Fourier amplitudes of the field.
 * @param[in, out] out The `output` structure providing information about and
 * memory for the power spectrum of the field.
 *
 * The power spectrum is constructed by binning the Fourier modes according to
 * the length of their wave vectors \f$k\f$. The bin size is uniform in \f$k\f$
 * and the number of bins is determined by `POWER_SPECTRUM_BINS` in the
 * parameters.
 */
static void mk_power_spectrum(const fftw_complex *in, struct output out)
{
    #pragma omp parallel for
    for (size_t i = 0; i < out.dim; ++i) {
        out.tmp[i] = 0.0;
    }
    // starting with 1 to explicitly exclude constant average, no omp!
    for (size_t i = 1; i < pars.M; ++i) {
        double pow2_tmp = in[i] * conj(in[i]);
        if (fabs(kvec.z[i]) > DBL_EPSILON) {
            pow2_tmp *= 2.0;
        }
        size_t idx = (int)trunc(out.dim * sqrt(kvec.sq[i]) / kvec.k_max - 1.0e-14);
        out.tmp[idx] += pow2_tmp / pars.N;
    }
}
#endif

#ifdef ENABLE_FFT_FILTER
/**
 * @brief Applies a filter to the Fourier modes of the real space input to
 * cutoff high frequency modes.
 *
 * @param[in, out] f The field which we want to filter. It is overwritten by the
 * filtered version.
 *
 * @see The highest modes of the field are cut off according to
 * `filter_window(const double x)` in `setup.c`.
 */
static void apply_filter(double *f)
{
    TIME(mon.filter -= get_wall_time());
    fft(f, tmp.fc);
    TIME(mon.filter += get_wall_time());
    capply_filter(tmp.fc, f);
}

/**
 * @brief Applies a filter to the Fourier space input to cutoff high frequency
 * modes.
 *
 * @param[in, out] in The field in Fourier space which we want to filter. It is
 * overwritten by the filtered version.
 * @param[out] out The array that is populated by the filtered field in real
 * space generated from @p in.
 *
 * @see The highest modes of the field are cut off according to
 * `filter_window(const double x)` in `setup.c`.
 */
static void capply_filter(complex *in, double *out)
{
    TIME(mon.filter -= get_wall_time());
    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        in[i] *= filter[i];
        tmp.deltarhoc[i] = in[i] / pars.N;
    }
    ifft(tmp.deltarhoc, out);
    TIME(mon.filter += get_wall_time());
}
#endif

/**
 * @brief Constructs and saves summaries of the required fields (containing the
 * mean, variance, minimum and maximum value at the current time) to the
 * corresponding buffers.
 */
static void mk_summary()
{
    TIME(mon.smry -= get_wall_time());
    #ifdef OUTPUT_PHI_SMRY
    mean_var_min_max(field, phi_smry.tmp, pars.N);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    mean_var_min_max(field + pars.N, dphi_smry.tmp, pars.N);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    mean_var_min_max(field + 2 * pars.N, psi_smry.tmp, pars.N);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    mean_var_min_max(field + 3 * pars.N, dpsi_smry.tmp, pars.N);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    mean_var_min_max(rho, rho_smry.tmp, pars.N);
    #endif
    #ifdef OUTPUT_PRESSURE_SMRY
    mean_var_min_max(pressure, p_smry.tmp, pars.N);
    #endif
    // TODO: when to compute summary of h1 and h2, need it in real space
    #ifdef OUTPUT_H1_SMRY
    fmean_var_min_max(field + 4 * pars.N, h1_smry.tmp);
    #endif
    #ifdef OUTPUT_H2_SMRY
    fmean_var_min_max(field + 4 * pars.N + pars.Next, h2_smry.tmp);
    #endif
    TIME(mon.smry += get_wall_time());
}

/**
 * @brief Computes the summary of a field, i.e. the mean, variance, minimum and
 * maximum value.
 *
 * @param[in] f The input field.
 * @param[out] smry An array of size 4 which is filled with the summary: mean,
 * variance, min, max (in this order).
 * @param[in] N The length of the input field @p N.
 *
 * @note The field is implicitly assumed to have length `pars.N`.
 */
static void mean_var_min_max(const double *f, double *smry, const size_t N)
{
    double mean = 0.0, min_val = f[0], max_val = f[0];
    #pragma omp parallel for reduction(+: mean) reduction(max: max_val) \
        reduction(min: min_val)
    for (size_t i = 0; i < N; ++i) {
        mean += f[i];
        min_val = MIN(min_val, f[i]);
        max_val = MAX(max_val, f[i]);
    }
    smry[0] = mean / N;
    smry[1] = variance(smry[0], f, N);
    smry[2] = min_val;
    smry[3] = max_val;
}

#if defined(OUTPUT_H1_SMRY) || defined(OUTPUT_H2_SMRY)
/**
 * @brief Compute the summary of a vector, i.e. the mean, variance, minimum and
 * maximum value given its Fourier transform
 *
 * @param[in] f The input vector in the Fourier domain
 * @param[out] smry An array of size 4 which is filled with the summary: mean,
 * variance, min, max (in this order)
 *
 * @note The vector is implicitly assumed to have length `pars.Next` and is
 * given in the memory layout described in TODO[link]
 */
static void fmean_var_min_max(const double *f, double *smry)
{
    //TODO: this is still heavily flawed
    double var = 0.0;
    #pragma omp parallel for reduction(+: var)
    for (size_t i = 1; i < pars.Next; ++i) {
        var += f[i] * f[i];
    }
    smry[0] = f[0];
    smry[1] = var / pars.N;
    smry[2] = 0.0;
    smry[3] = 0.0;
}
#endif

/**
 * @brief Compute the mean (or average) of a 1D array.
 *
 * @param[in] f Any 1D array of length @p N.
 * @param[in] N The length of the 1D array @p f.
 * @return The mean value (average) of @p f.
 */
static double mean(const double *f, const size_t N)
{
    double mean = 0.0;
    #pragma omp parallel for reduction(+: mean)
    for (size_t i = 0; i < N; ++i) {
        mean += f[i];
    }
    return mean / N;
}

/**
 * @brief Compute the variance of a 1D array.
 *
 * @param[in] mean The mean of the 1D array @p f.
 * @param[in] f The 1D array of length @p N.
 * @param[in] N The length of the 1D array @p f.
 * @return The variance of @p f.
 */
static double variance(const double mean, const double *f, const size_t N)
{
    double sum1 = 0.0, sum2 = 0.0, tmp;
    #pragma omp parallel for private(tmp) reduction(+: sum1, sum2)
    for (size_t i = 0; i < N; ++i) {
        tmp = f[i] - mean;
        sum1 += tmp * tmp;
        sum2 += tmp;
    }
    return (sum1 - sum2 * sum2 / N) / (N - 1);
}

/**
 * @brief Copy complex field saved in the custom double memory layout to a
 * complex array.
 *
 * @param[in] in The double input array.
 * @param[out] out The complex output array.
 */
static void real_to_complex(const double *in, complex *out)
{
    #pragma omp parallel for
    for (size_t i = 0; i < pars.Next; i += 2) {
        out[i] = in[i] + I * in[i + 1];
    }
}

/**
 * @brief Copy complex array to the custom memory layout in a double layout.
 *
 * @param[in] in The complex input array.
 * @param[out] out The double output array.
 */
static void complex_to_real(const complex *in, double *out)
{
    #pragma omp parallel for
    for (size_t i = 0; i < pars.M; ++i) {
        size_t ii = 2 * i;
        out[ii] = creal(in[i]);
        out[ii + 1] = cimag(in[i]);
    }
}

/**
 * @brief Fourier transform
 *
 * @param[in] in Pointer to the function in real space
 * @param[out] out Pointer to an array for the Fourier transform
 */
static void fft(double *in, complex *out)
{
    TIME(mon.fftw_exe -= get_wall_time());
    fftw_execute_dft_r2c(p_fw, in, out);
    TIME(mon.fftw_exe += get_wall_time());
}

/**
 * @brief Inverse Fourier transform
 *
 * @param[in] in Pointer to the function in Fourier space
 * @param[out] out Pointer to an array for the inverse Fourier transform
 */
static void ifft(complex *in, double *out)
{
    TIME(mon.fftw_exe -= get_wall_time());
    fftw_execute_dft_c2r(p_bw, in, out);
    TIME(mon.fftw_exe += get_wall_time());
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
