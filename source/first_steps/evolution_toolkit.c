#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"
#include "filehandling.h"

// see header file for more info
evolution_flags_t evo_flags = {.filter = 0, .compute_pow_spec = 0};

// compute the right hand side of the pde, i.e. the first temporal derivatives
// of all fields (scalar field, its first temporal derivative and a)
void mk_rhs(const double t, double *f, double *result) {
    size_t N = pars.N;
    size_t N2 = 2 * N;
    double a = f[N2];

    rho_avg = mk_rho(f);
    double hubble = sqrt(rho_avg / 3.0);

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        result[i] = f[N + i];
    }

    #pragma omp parallel for
    for (size_t i = N; i < N2; ++i)
    {
        result[i] = dtmp_lap[i - N] / (a * a);
        result[i] -= ( 3.0 * hubble * f[i]
                        + potential_prime(f[i - N]) );
    }
    result[N2] = a * hubble;
}

// compute energy density rho, i.e. 00 of stress energy, return average value
double mk_rho(double *f) {
    size_t N = pars.N;
    double a = f[2 * N];
    double T00 = 0.0;

    mk_gradient_squared_and_laplacian(f);

    double df, grad2_a;
    #pragma omp parallel for default(shared) private(df, grad2_a) reduction(+:T00)
    for (size_t i = 0; i < N; ++i)
    {
        df = f[N + i];
        grad2_a = dtmp_grad2[i] / (a * a);
        rho[i] = (df * df + grad2_a) / 2. + potential(f[i]);
        T00 += rho[i];
    }
    return T00 / N;
}

// compute the laplacian and the squared gradient of the input and store them in
// dtmp_lap and dtmp_grad2
void mk_gradient_squared_and_laplacian(double *in) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t ncz = Nz / 2 + 1;

    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw_3d, in, cfftw_tmp);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start;
    #endif

    if (evo_flags.compute_pow_spec == 1)
    {
        mk_power_spectrum(cfftw_tmp);
    }

    // TODO[performance]: precompute these factors only once and reuse them
    double Lx = pars.x.L;
    double Ly = pars.y.L;
    double Lz = pars.z.L;

    complex prefac = 2. * PI * I;
    complex factor_x = prefac / (Lx * N);
    complex factor_y = prefac / (Ly * N);
    complex factor_z = prefac / (Lz * N);

    double prefac2 = -4. * PI * PI;
    double factor_x2 = prefac2 / (Lx * Lx * N);
    double factor_y2 = prefac2 / (Ly * Ly * N);
    double factor_z2 = prefac2 / (Lz * Lz * N);
    double k_sq;

    size_t osx, osy, id;
    #pragma omp parallel for private(osx, osy, id, k_sq)
    for (size_t i = 0; i < Nx; ++i)
    {
        osx = i * Ny * ncz;
        for (size_t j = 0; j < Ny; ++j)
        {
            osy = osx + j * ncz;
            for (size_t k = 0; k < ncz; ++k)
            {
                id = osy + k;

                // for lagrangian
                k_sq = factor_z2 * k * k;

                // x derivative
                if (i > Nx / 2)
                {
                    cfftw_tmp_x[id] = cfftw_tmp[id] * factor_x
                                                        * ((int)i - (int)Nx);
                    k_sq += factor_x2 * (Nx - i) * (Nx - i);
                }
                else if (2 * i == Nx)
                {
                    cfftw_tmp_x[id] = 0.0;
                    k_sq += factor_x2 * i * i;
                }
                else
                {
                    cfftw_tmp_x[id] = cfftw_tmp[id] * factor_x * i;
                    k_sq += factor_x2 * i * i;
                }

                // y derivative
                if (j > Ny / 2)
                {
                    cfftw_tmp_y[id] = cfftw_tmp[id] * factor_y
                                                        * ((int)j - (int)Ny);
                    k_sq += factor_y2 * (Ny - j) * (Ny - j);
                }
                else if (2 * j == Ny)
                {
                    cfftw_tmp_y[id] = 0.0;
                    k_sq += factor_y2 * j * j;
                }
                else
                {
                    cfftw_tmp_y[id] = cfftw_tmp[id] * factor_y * j;
                    k_sq += factor_y2 * j * j;
                }

                // z derivative
                if (k == ncz - 1)
                {
                    cfftw_tmp_z[id] = 0.0;
                }
                else
                {
                    cfftw_tmp_z[id] = cfftw_tmp[id] * factor_z * k;
                }

                // lagrangian
                cfftw_tmp[id] *= k_sq;
            }
        }
    }

    #ifdef SHOW_TIMING_INFO
    start = get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_x, dtmp_x);
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_y, dtmp_y);
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_z, dtmp_z);
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp, dtmp_lap);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start;
    #endif

    // gradient squared
    double gx, gy, gz;
    #pragma omp parallel for private(gx, gy, gz)
    for (size_t i = 0; i < N; ++i)
    {
        gx = dtmp_x[i];
        gy = dtmp_y[i];
        gz = dtmp_z[i];
        dtmp_grad2[i] = gx * gx + gy * gy + gz * gz;
    }
}

// A selection of potentials one can try, make sure to set the corresponding
// potential_prime, the derivative is not computed automatically
// TODO: change that?
inline double potential(const double f) {
    // higgs metastability potential
    double l = LAMBDA / 1.0e-10;
    double a = 0.01 * l, b = 1.0 * l;
    return f == 0.0 ? LAMBDA :
        (LAMBDA + a * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) * pow(f, 4) +
        b * pow(f, 6));

    // notch or step potential (LAMBDA = 1.876e-4)
    // double lambda = 100.0;
    // return LAMBDA / (1.0 + exp(-lambda * f));

    // standard f squared potential
    // return MASS * MASS * f * f / 2.0;

    // standard f to the fourth (with f squared) potential
    // return MASS * MASS * f * f / 2.0 + COUPLING * f * f * f * f / 24.0;

    // return 0.0;
}

inline double potential_prime(const double f) {
    // higgs metastability potential
    double l = LAMBDA / 1.0e-10;
    double a = 0.01 * l, b = 1.0 * l;
    return f == 0.0 ? 0 :
        (4.0 * a * pow(f, 3) * (1.0 - 0.1 * log10(fabs(f) * 1.0e19)) -
        (0.1 * a * pow(f, 4) * ((f > 0.0) - (f < 0.0))) / (fabs(f) * log(10.0))
        + 6.0 * b * pow(f, 5));

    // notch or step potential (LAMBDA = 1.876e-4)
    // double lambda = 100.0;
    // double tmp = exp(lambda * f);
    // return LAMBDA * lambda * tmp / ((1.0 + tmp) * (1.0 + tmp));

    // standard f squared potential
    // return MASS * MASS * f;

    // standard f to the fourth (with f squared) potential
    // return MASS * MASS * f + COUPLING * f * f * f / 6.0;

    // return 0.0;
}

void solve_poisson_eq() {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t ncz = Nz / 2 + 1;

    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw_3d, rho, cfftw_tmp);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start;
    #endif
    
    // TODO[performance]: precompute these factors only once and reuse them
    double Lx = pars.x.L;
    double Ly = pars.y.L;
    double Lz = pars.z.L;

    double prefac2 = -4. * PI * PI;
    double factor_x2 = prefac2 * N / (Lx * Lx);
    double factor_y2 = prefac2 * N / (Ly * Ly);
    double factor_z2 = prefac2 * N / (Lz * Lz);
    double k_sq;

    size_t osx, osy;
    #pragma omp parallel for private(osx, osy, k_sq)
    for (size_t i = 0; i < Nx; ++i)
    {
        osx = i * Ny * ncz;
        for (size_t j = 0; j < Ny; ++j)
        {
            osy = osx + j * ncz;
            for (size_t k = 0; k < ncz; ++k)
            {
                k_sq = factor_z2 * k * k;
                if (i > Nx / 2)
                {
                    k_sq += factor_x2 * (Nx - i) * (Nx - i);
                }
                else
                {
                    k_sq += factor_x2 * i * i;
                }
                if (j > Ny / 2)
                {
                    k_sq += factor_y2 * (Ny - j) * (Ny - j);
                }
                else
                {
                    k_sq += factor_y2 * j * j;
                }
                if (k_sq > 1.0e-16)
                {
                    // factor two comes from equation: grad^2 psi = rho / 2
                    // using 8 pi G = 1
                    cfftw_tmp[osy + k] /= 2.0 * k_sq;
                }
                else
                {
                    cfftw_tmp[osy + k] = 0.0;
                }
            }
        }
    }

    #ifdef SHOW_TIMING_INFO
    start = get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp, psi);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start;
    #endif
}

// computes a crude estimation of the power spectrum, more info in main.h
// and stores it in global pow_spec
void mk_power_spectrum(const fftw_complex *in) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t N = pars.N;
    size_t ncz = Nz / 2 + 1;
    size_t bins = pars.file.bins_powspec;

    // todo[performance]: precompute bins only once and reuse
    double Lx = pars.x.L;
    double Ly = pars.y.L;
    double Lz = pars.z.L;

    double k_x2 = 1.0 / (Lx * Lx);
    double k_y2 = 1.0 / (Ly * Ly);
    double k_z2 = 1.0 / (Lz * Lz);

    double k_max2 = k_x2 * (Nx/2) * (Nx/2) +
                    k_y2 * (Ny/2) * (Ny/2) +
                    k_z2 * (Nz/2) * (Nz/2);

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
    for (size_t i = 0; i < Nx; ++i)
    {
        osx = i * Ny * ncz;
        for (size_t j = 0; j < Ny; ++j)
        {
            osy = osx + j * ncz;
            for (size_t k = 0; k < ncz; ++k)
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
                    k2_tmp = k_z2 * k * k;
                    if (i > Nx / 2)
                    {
                        k2_tmp += k_x2 * (Nx - i) * (Nx - i);
                    }
                    else
                    {
                        k2_tmp += k_x2 * i * i;
                    }

                    if (j > Ny / 2)
                    {
                        k2_tmp += k_y2 * (Ny - j) * (Ny - j);
                    }
                    else
                    {
                        k2_tmp += k_y2 * j * j;
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
    double start_filter = get_wall_time();
    #endif

    #ifdef SHOW_TIMING_INFO
    double start_fft = get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw_3d, inout, cfftw_tmp);
    fftw_execute_dft_r2c(p_fw_3d, inout + N, cfftw_tmp_x);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start_fft;
    #endif
    
    apply_filter_fourier(cfftw_tmp, cfftw_tmp_x);

    #ifdef SHOW_TIMING_INFO
    start_fft = get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp, inout);
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_x, inout + N);
    #ifdef SHOW_TIMING_INFO
    fftw_time_exe += get_wall_time() - start_fft;
    #endif

    #ifdef SHOW_TIMING_INFO
    filter_time += get_wall_time() - start_filter;
    #endif
}

// filtering in fourier domain for two fields (phi and dphi) simultaneously
void apply_filter_fourier(fftw_complex *inout, fftw_complex *dinout) {
    size_t N = pars.N;
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t ncz = Nz / 2 + 1;

    double x, y, z;
    size_t osx, osy;
    #pragma omp parallel for private(osx, osy)
    for (size_t i = 0; i < Nx; ++i)
    {
        osx = i * Ny * ncz;
        for (size_t j = 0; j < Ny; ++j)
        {
            osy = osx + j * ncz;
            for (size_t k = 0; k < ncz; ++k)
            {
                x = filter_window_function(2.0 *
                    (i > Nx / 2 ? Nx - i : i) / (double) Nx);
                y = filter_window_function(2.0 *
                    (j > Ny / 2 ? Ny - j : j) / (double) Ny);
                z = filter_window_function(2.0 * k / (double) Nz) / (double) N;
                inout[osy + k] *= x * y * z; 
                dinout[osy + k] *= x * y * z; 
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

    // miscellaneous
    // return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. );
    // return exp(1. + 1. / ( pow(x, 8) - 1. ));
    // return 0.5 * ( 1. + cos( pow(x, 8) * PI ) );
    // return 0.0;
}

// recompute current power spectrum and rho and save current timeslice to buffer
// we might want to save this extra call to to mk_rho and write out data at a
// point where everything is available anyway
void prepare_and_save_timeslice() {
    evo_flags.compute_pow_spec = 1;
    rho_avg = mk_rho(field);
    evo_flags.compute_pow_spec = 0;
    save();
}
