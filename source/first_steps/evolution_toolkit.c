#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <fftw3.h>
#include "evolution_toolkit.h"
#include "main.h"
#include "filehandling.h"

evolution_flags_t evo_flags = {.filter = 0, .compute_pow_spec = 0};

void mk_gradient_squared_and_laplacian(double *in, double *grad2,
                                            double *laplacian) {
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
    double end = get_wall_time();
    fftw_time_exe += end - start;
    #endif

    if (evo_flags.compute_pow_spec == 1)
    {
        mk_power_spectrum(cfftw_tmp);
    }

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
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp,  laplacian);
    #ifdef SHOW_TIMING_INFO
    end = get_wall_time();
    fftw_time_exe += end - start;
    #endif

    double gx, gy, gz;
    #pragma omp parallel for private(gx, gy, gz)
    for (size_t i = 0; i < N; ++i)
    {
        gx = dtmp_x[i];
        gy = dtmp_y[i];
        gz = dtmp_z[i];
        grad2[i] = gx * gx + gy * gy + gz * gz;
    }
}

void mk_power_spectrum(fftw_complex *in) {
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

void apply_filter_real(double *inout) {
    size_t N = pars.N;

    #ifdef SHOW_TIMING_INFO
    double start = get_wall_time();
    #endif
    fftw_execute_dft_r2c(p_fw_3d, inout, cfftw_tmp);
    fftw_execute_dft_r2c(p_fw_3d, inout + N, cfftw_tmp_x);
    #ifdef SHOW_TIMING_INFO
    double end = get_wall_time();
    fftw_time_exe += end - start;
    #endif
    
    apply_filter_fourier(cfftw_tmp);
    apply_filter_fourier(cfftw_tmp_x);

    #ifdef SHOW_TIMING_INFO
    start = get_wall_time();
    #endif
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp, inout);
    fftw_execute_dft_c2r(p_bw_3d, cfftw_tmp_x, inout + N);
    #ifdef SHOW_TIMING_INFO
    end = get_wall_time();
    fftw_time_exe += end - start;
    #endif
}

void apply_filter_fourier(fftw_complex *inout) {
    double cutoff_fraction = pars.cutoff_fraction;
    if (cutoff_fraction < 0.0 || cutoff_fraction > 1.0)
    {
        fputs("cutoff_fraction must be between 0 and 1.", stderr);
        exit(EXIT_FAILURE);
    }

    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
    size_t ncz = Nz / 2 + 1;

    size_t i_a = (size_t) floor((Nx / 2) * (1.0 - cutoff_fraction));
    size_t i_b = Nx - i_a;
    size_t j_a = (size_t) floor((Ny / 2) * (1.0 - cutoff_fraction));
    size_t j_b = Ny - j_a;
    size_t k_a = (size_t) floor((Nz / 2) * (1.0 - cutoff_fraction));

    // TODO[performance]: run loop only over _a to _b in each dimension for two
    // thirds rule, leave it like that for possible fourier smoothing
    size_t osx, osy, nok = 0, nbad = 0;
    #pragma omp parallel for private(osx, osy) reduction(+:nok, nbad)
    for (size_t i = 0; i < Nx; ++i)
    {
        osx = i * Ny * ncz;
        for (size_t j = 0; j < Ny; ++j)
        {
            osy = osx + j * ncz;
            for (size_t k = 0; k < ncz; ++k)
            {
                if ( k >= k_a ||
                    (j >= j_a && j <= j_b) ||
                    (i >= i_a && i <= i_b) )
                {
                    inout[osy + k] = 0.0;
                    ++nbad;
                }
                else
                {
                    ++nok;
                }
            }
        }
    }
    RUNTIME_INFO(printf("filtered: %f\n", (double)nbad / (double)(nok + nbad)));
}

/*
compute the right hand side of the pde, ie the first order temporal derivatives
*/
void mk_velocities(double t, double *f, double *result) {
    size_t N = pars.N;
    size_t N2 = 2 * N;
    double a = f[N2];

    rho = mk_rho(f);
    double hubble = sqrt(rho / 3.0);

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

/*
A selection of potentials one can try, make sure to set the corresponding
potential_prime, the derivative is not computed automatically yet
TODO: change that?
*/
inline double potential(double f){
    // return 0.0;
    double lambda = 100.0;
    return LAMBDA / (1.0 + exp(-lambda * f));

    // return MASS * MASS * f * f / 2.0;

    // return MASS * MASS * f * f / 2.0 + COUPLING * f * f * f * f / 24.0;

    // return 0.0;
}

inline double potential_prime(double f) {
    // return 0.0;
    double lambda = 100.0;
    double tmp = exp(lambda * f);
    return LAMBDA * lambda * tmp / ((1.0 + tmp) * (1.0 + tmp));

    // return MASS * MASS * f;

    // return MASS * MASS * f + COUPLING * f * f * f / 6.0;

    // return 20.0 * tanh(pow(f, 50));

    // return 0.0;
}

/*
compute average 00 component of stress energy
*/
double mk_rho(double *f) {
    size_t N = pars.N;
    double a = f[2 * N];
    double T00 = 0.0;

    mk_gradient_squared_and_laplacian(f, dtmp_grad2, dtmp_lap);

    double df, grad2_a;
    #pragma omp parallel for default(shared) private(df, grad2_a) reduction(+:T00)
    for (size_t i = 0; i < N; ++i)
    {
        df = f[N + i];
        grad2_a = dtmp_grad2[i] / (a * a);
        T00 += (df * df + grad2_a) / 2. + potential(f[i]);
    }
    return T00 / N;
}

void mk_filter_window(double *inout, size_t cutoffindex, size_t windowlength) {

    #pragma omp parallel for
    for (size_t i = 0; i < cutoffindex; ++i)
    {
        inout[i] = filter_window_function((double)i / cutoffindex);
    }
    #pragma omp parallel for
    for (size_t i = cutoffindex; i < windowlength; ++i)
    {
        inout[i] = 0.0;
    }
}

inline double filter_window_function(double x) {
    return exp(-36.0 * pow(x, 36));
    // return 1. - tanh( 1. / ( 1. - pow(x, 8) ) - 1. );
    // return exp(1. + 1. / ( pow(x, 8) - 1. ));
    // return 0.5 * ( 1. + cos( pow(x, 8) * PI ) );
    // return 0.0;
}
