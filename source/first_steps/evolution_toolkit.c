#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"
#include "filehandling.h"
#include "setup.h"

// see header file for more info
evolution_flags_t evo_flags = {.filter = 0, .compute_pow_spec = 0};

// compute the right hand side of the pde, i.e. the first temporal derivatives
// of all fields (scalar field, its first temporal derivative and a)
void mk_rhs(const double t, double *f, double *result) {
    size_t N = pars.N;
    size_t N2 = 2 * N;
    double a = f[N2];

    mk_rho(f);
    double hubble = sqrt(rho_avg / 3.0);

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        result[i] = f[N + i];
    }

    #ifdef INCLUDE_PSI
    mk_psi_and_dpsi(f);
    double df, p;
    #pragma omp parallel for private(df, p)
    for (size_t i = 0; i < N; ++i)
    {
        df = f[N + i];
        p = psi[i];
        result[N + i] = (1.0 + 4.0 * p) * tmp_phi.lap[i] / (a * a) -
            (3.0 * hubble - 4.0 * dpsi[i]) * df -
            (1.0 + 2.0 * p) * potential_prime(f[i]);
    }
    #else
    #pragma omp parallel for
    for (size_t i = N; i < N2; ++i)
    {
        result[i] = tmp_phi.lap[i - N] / (a * a);
        result[i] -= (3.0 * hubble * f[i] + potential_prime(f[i - N]));
    }
    #endif
    result[N2] = a * hubble;
}

// compute energy density rho, i.e. 00 of stress energy & save average value
void mk_rho(double *f) {
    size_t N = pars.N;
    double a = f[2 * N];
    rho_avg = 0.0;

    mk_gradient_squared_and_laplacian(f);

    double df, grad2;
    #pragma omp parallel for default(shared) private(df, grad2) \
        reduction(+: rho_avg)
    for (size_t i = 0; i < N; ++i)
    {
        df = f[N + i];
        grad2 = tmp_phi.grad[i] / (a * a);
        rho[i] = (df * df + grad2) / 2.0 + potential(f[i]);
        rho_avg += rho[i];
    }
    rho_avg /= N;
}

// compute the laplacian and the squared gradient of the input and store them
void mk_gradient_squared_and_laplacian(double *in) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t N  = pars.N;
    size_t Mx = pars.x.M;
    size_t My = pars.y.M;
    size_t Mz = pars.z.M;

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw, in, tmp_phi.c);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    if (evo_flags.compute_pow_spec == 1)
    {
        mk_power_spectrum(tmp_phi.c);
    }

    double k_sq;
    size_t osx, osy, id;
    #pragma omp parallel for private(osx, osy, id, k_sq)
    for (size_t i = 0; i < Mx; ++i)
    {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j)
        {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k)
            {
                id = osy + k;
                // for laplacian
                k_sq = pars.z.k2 * k * k;
                // x derivative
                if (i > Nx / 2)
                {
                    tmp_phi.cx[id] = tmp_phi.c[id] * pars.x.k
                        * ((int)i - (int)Nx) / N;
                    k_sq += pars.x.k2 * (Nx - i) * (Nx - i);
                }
                else if (2 * i == Nx)
                {
                    tmp_phi.cx[id] = 0.0;
                    k_sq += pars.x.k2 * i * i;
                }
                else
                {
                    tmp_phi.cx[id] = tmp_phi.c[id] * pars.x.k * i / N;
                    k_sq += pars.x.k2 * i * i;
                }
                // y derivative
                if (j > Ny / 2)
                {
                    tmp_phi.cy[id] = tmp_phi.c[id] * pars.y.k
                        * ((int)j - (int)Ny) / N;
                    k_sq += pars.y.k2 * (Ny - j) * (Ny - j);
                }
                else if (2 * j == Ny)
                {
                    tmp_phi.cy[id] = 0.0;
                    k_sq += pars.y.k2 * j * j;
                }
                else
                {
                    tmp_phi.cy[id] = tmp_phi.c[id] * pars.y.k * j / N;
                    k_sq += pars.y.k2 * j * j;
                }
                // z derivative
                if (k == Mz - 1)
                {
                    tmp_phi.cz[id] = 0.0;
                }
                else
                {
                    tmp_phi.cz[id] = tmp_phi.c[id] * pars.z.k * k / N;
                }
                // laplacian
                tmp_phi.c[id] *= k_sq / N;
            }
        }
    }

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw, tmp_phi.cx, tmp_phi.dx);
    if (pars.dim > 1)
    {
        fftw_execute_dft_c2r(p_bw, tmp_phi.cy, tmp_phi.dy);
    }
    if (pars.dim > 2)
    {
        fftw_execute_dft_c2r(p_bw, tmp_phi.cz, tmp_phi.dz);
    }
    fftw_execute_dft_c2r(p_bw, tmp_phi.c, tmp_phi.lap);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    // gradient squared
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        tmp_phi.grad[i] = tmp_phi.dx[i] * tmp_phi.dx[i];
        if (pars.dim > 1)
        {
            tmp_phi.grad[i] += tmp_phi.dy[i] * tmp_phi.dy[i];
        }
        if (pars.dim > 2)
        {
            tmp_phi.grad[i] += tmp_phi.dz[i] * tmp_phi.dz[i];
        }
    }
}

// A selection of potentials one can try, make sure to set the corresponding
// potential_prime, the derivative is not computed automatically
// TODO: change that?
inline double potential(const double f) {
    // higgs metastability potential
    /* double l = LAMBDA / 1.0e-10; */
    /* double a = 0.01 * l, b = 1.0 * l; */
    /* return f == 0.0 ? LAMBDA : */
    /*     (LAMBDA + a * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) * pow(f, 4) + */
    /*     b * pow(f, 6)); */

    // notch or step potential
    // LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5
    double lambda = 100.0;
    return LAMBDA / (1.0 + exp(-lambda * f));

    // standard f squared potential
    /* return MASS * MASS * f * f / 2.0; */

    // standard f to the fourth (with f squared) potential
    /* return MASS * MASS * f * f / 2.0 + COUPLING * f * f * f * f / 24.0; */

    /* return LAMBDA; */

    /* return 0.0; */
}

inline double potential_prime(const double f) {
    // higgs metastability potential
    /* double l = LAMBDA / 1.0e-10; */
    /* double a = 0.01 * l, b = 1.0 * l; */
    /* return f == 0.0 ? 0 : */
    /*     (4.0 * a * pow(f, 3) * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) - */
    /*     (0.1 * a * pow(f, 4) * ((f > 0.0) - (f < 0.0))) / (fabs(f) * log(10.0)) */
    /*     + 6.0 * b * pow(f, 5)); */

    // notch or step potential
    // LAMBDA = 3d: 1.876e-4, 2d: 4.721e-5, 1d: 4.1269e-5
    double lambda = 100.0;
    double tmp = exp(lambda * f);
    return LAMBDA * lambda * tmp / ((1.0 + tmp) * (1.0 + tmp));

    // standard f squared potential
    /* return MASS * MASS * f; */

    // standard f to the fourth (with f squared) potential
    /* return MASS * MASS * f + COUPLING * f * f * f / 6.0; */

    return 0.0;
}

// solve the poisson like equation for scalar perturbations and simultaneously
// solving for its derivative
void mk_psi_and_dpsi(double *f) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t N = pars.N;
    size_t Mx = pars.x.M;
    size_t My = pars.y.M;
    size_t Mz = pars.z.M;
    size_t M = Mx * My * Mz;
    double a = f[2 * N];
    double a2 = a * a;
    double hubble = sqrt(rho_avg / 3.0);

    #ifdef SHOW_TIMING_INFO
    poisson_time -= get_wall_time();
    #endif

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        tmp_psi.dx[i] = f[N + i] * tmp_phi.dx[i];
        tmp_psi.dy[i] = rho[i] - rho_avg;
    }

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw, tmp_psi.dx, tmp_psi.cx);
    fftw_execute_dft_r2c(p_fw, tmp_psi.dy, tmp_psi.cy);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    double k_sq;
    size_t osx, osy, id;
    #pragma omp parallel for private(k_sq, osx, osy, id)
    for (size_t i = 0; i < Mx; ++i)
    {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j)
        {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k)
            {
                id = osy + k;
                k_sq = pars.z.k2 * k * k;
                if (i > Nx / 2)
                {
                    k_sq += pars.x.k2 * (Nx - i) * (Nx - i);
                    tmp_psi.cx[id] /= pars.x.k * ((int)i - (int)Nx);
                }
                else if (2 * i == Nx || i == 0)
                {
                    k_sq += pars.x.k2 * i * i;
                    tmp_psi.cx[id] = 0.0;
                }
                else
                {
                    k_sq += pars.x.k2 * i * i;
                    tmp_psi.cx[id] /= pars.x.k * i;
                }
                if (j > Ny / 2)
                {
                    k_sq += pars.y.k2 * (Ny - j) * (Ny - j);
                }
                else
                {
                    k_sq += pars.y.k2 * j * j;
                }
                if (fabs(k_sq) > 1.0e-12)
                {
                    tmp_psi.c[id] = 0.5 * a2 *
                        (tmp_psi.cy[id] + 3.0 * hubble * tmp_psi.cx[id]) /
                        (k_sq * N);
                }
                else
                {
                    tmp_psi.c[id] = 0.0;
                }
                tmp_psi.cz[id] = 0.5 * tmp_psi.cx[id] / N -
                        hubble * tmp_psi.c[id];
            }
        }
    }

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw, tmp_psi.c, psi);
    fftw_execute_dft_c2r(p_bw, tmp_psi.cz, dpsi);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    poisson_time += get_wall_time();
    #endif
}

// computes a crude estimation of the power spectrum, more info in main.h
// and stores it in global pow_spec
void mk_power_spectrum(const fftw_complex *in) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t bins = pars.file.bins_powspec;
    size_t Mx = pars.x.M;
    size_t My = pars.y.M;
    size_t Mz = pars.z.M;

    double k_max2 = pars.x.k2 * (Nx/2) * (Nx/2) +
                    pars.y.k2 * (Ny/2) * (Ny/2) +
                    pars.z.k2 * (Nz/2) * (Nz/2);
    double dk2 = k_max2 / bins;
    double k2_tmp = 0.0;
    double pow2_tmp = 0.0;

    #pragma omp parallel for
    for (size_t i = 0; i < bins; ++i)
    {
        pow_spec[i] = 0.0;
    }

    size_t osx, osy, idx;
    // todo parallelize?
    for (size_t i = 0; i < Mx; ++i)
    {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j)
        {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k)
            {
                if (k == 0 || 2 * k == Nz)
                {
                    pow2_tmp = in[osy + k] * conj(in[osy + k]);
                }
                else
                {
                    pow2_tmp = 2.0 * in[osy + k] * conj(in[osy + k]);
                }
                if (pow2_tmp > 0.0)
                {
                    k2_tmp = pars.z.k2 * k * k;
                    if (i > Nx / 2)
                    {
                        k2_tmp += pars.x.k2 * (Nx - i) * (Nx - i);
                    }
                    else
                    {
                        k2_tmp += pars.x.k2 * i * i;
                    }

                    if (j > Ny / 2)
                    {
                        k2_tmp += pars.y.k2 * (Ny - j) * (Ny - j);
                    }
                    else
                    {
                        k2_tmp += pars.y.k2 * j * j;
                    }
                    idx = (int)(k2_tmp / dk2 - 1e-10);
                    pow_spec[idx] += pow2_tmp / N;
                }
            }
        }
    }
}

// filter the real input field, input gets overwritten with filtered data
void apply_filter_real(double *inout) {
    size_t N = pars.N;

    #ifdef SHOW_TIMING_INFO
    filter_time -= get_wall_time();
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw, inout, tmp_phi.c);
    fftw_execute_dft_r2c(p_fw, inout + N, tmp_phi.cx);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    #endif

    apply_filter_fourier(tmp_phi.c, tmp_phi.cx);

    #ifdef SHOW_TIMING_INFO
    fftw_time_exe -= get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw, tmp_phi.c, inout);
    fftw_execute_dft_c2r(p_bw, tmp_phi.cx, inout + N);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time();
    filter_time += get_wall_time();
    #endif
}

// filtering in fourier domain for two fields (phi and dphi) simultaneously
void apply_filter_fourier(fftw_complex *inout, fftw_complex *dinout) {
    size_t N = pars.N;
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t Mx = pars.x.M;
    size_t My = pars.y.M;
    size_t Mz = pars.z.M;

    double x, y, z;
    size_t osx, osy;
    #pragma omp parallel for private(osx, osy)
    for (size_t i = 0; i < Mx; ++i)
    {
        osx = i * My * Mz;
        for (size_t j = 0; j < My; ++j)
        {
            osy = osx + j * Mz;
            for (size_t k = 0; k < Mz; ++k)
            {
                x = filter_window_function(2.0 *
                    (i > Nx / 2 ? (int)Nx - (int)i : i) / (double) Nx);
                y = filter_window_function(2.0 *
                    (j > Ny / 2 ? (int)Ny - (int)j : j) / (double) Ny);
                z = filter_window_function(2.0 * k / (double) Nz);
                inout[osy + k] *= x * y * z / (double) N;
                dinout[osy + k] *= x * y * z / (double) N;
            }
        }
    }
}

// the cutoff function for filtering, use either two thirds or fourier smoothing
inline double filter_window_function(const double x) {
    // fourier smoothing
    return exp(-36.0 * pow(x, 36));

    // two thirds rule
    // return x < 2.0/3.0 ? x : 0.0;

    // miscellaneous, have been used in testing
    /* return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. ); */
    /* return exp(1. + 1. / ( pow(x, 8) - 1. )); */
    /* return 0.5 * ( 1. + cos( pow(x, 8) * PI ) ); */
    /* return 0.0; */
}

// recompute current power spectrum and rho and save current timeslice to buffer
// we might want to save this extra call to to mk_rho and write out data at a
// point where everything is available anyway
void prepare_and_save_timeslice() {
    evo_flags.compute_pow_spec = 1;
    mk_rho(field);
    evo_flags.compute_pow_spec = 0;
    save();
}

//TODO: do i need that function? delete?
double mean(const double *in, size_t N) {
    double mean = 0.0;
    #pragma omp parallel for reduction(+: mean)
    for (size_t i = 0; i < N; ++i)
    {
        mean += in[i];
    }
    return mean / N;
}

void contains_nan(double *f, size_t N) {
    size_t count = 0;
    for (size_t i = 0; i < N; ++i)
    {
        if (isnan(f[i]))
        {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}

void contains_nanc(complex *f, size_t N) {
    size_t count = 0;
    for (size_t i = 0; i < N; ++i)
    {
        if (isnan(creal(f[i])) || isnan(cimag(f[i])))
        {
            ++count;
        }
    }
    printf("found %zu nans\n", count);
}
