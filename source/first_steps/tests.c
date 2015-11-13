#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "main.h"
#include "setup.h"
#include "evolution_toolkit.h"

void run_all_tests(parameters_t *pars) {
    test_fft_first_derivative(pars);
    test_fft_second_derivative(pars);
    test_mk_laplacian(pars);
    test_mk_gradient(pars);
    test_fft_apply_filter(pars);
}

void test_fft_first_derivative(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;

    fft_D(field, dmisc_tmp_x, 1, 1, pars);
    fill_field(dmisc_tmp_x + Ntot, test_func_Dx, pars);
    puts("test first derivative in x direction:");
    if (are_fields_equal(dmisc_tmp_x, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testDx");
        print_vector(dmisc_tmp_x, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_x + Ntot, Ntot);
        puts("\n");
#endif

    // fft_D(field, dmisc_tmp_y, 2, 1, pars);
    // fill_field(dmisc_tmp_y + Ntot, test_func_Dy, pars);
    // puts("test first derivative in y direction:");
    // if (are_fields_equal(dmisc_tmp_y, pars) == 0)
    // {
    //     puts("passed\n");
    // }
    // else
    // {
    //     puts("failed\n");
    // }
// #ifdef DEBUG
//         puts("testDy");
//         print_vector(dmisc_tmp_y, Ntot);
//         puts("exact");
//         print_vector(dmisc_tmp_y + Ntot, Ntot);
//         puts("\n");
// #endif

    fft_D(field, dmisc_tmp_z, 3, 1, pars);
    fill_field(dmisc_tmp_z + Ntot, test_func_Dz, pars);
    puts("test first derivative in z direction:");
    if (are_fields_equal(dmisc_tmp_z, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testDz");
        print_vector(dmisc_tmp_z, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_z + Ntot, Ntot);
        puts("\n");
#endif
}

void test_fft_second_derivative(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;

    fft_D(field, dmisc_tmp_x, 1, 2, pars);
    fill_field(dmisc_tmp_x + Ntot, test_func_D2, pars);
    puts("test second derivative in x direction:");
    if (are_fields_equal(dmisc_tmp_x, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testDxx");
        print_vector(dmisc_tmp_x, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_x + Ntot, Ntot);
        puts("\n");
#endif

    // fft_D(field, dmisc_tmp_y, 2, 2, pars);
    // fill_field(dmisc_tmp_y + Ntot, test_func_D2, pars);
    // puts("test second derivative in y direction:");
    // if (are_fields_equal(dmisc_tmp_y, pars) == 0)
    // {
    //     puts("passed\n");
    // }
    // else
    // {
    //     puts("failed\n");
    // }
#ifdef DEBUG
        puts("testDyy");
        print_vector(dmisc_tmp_y, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_y + Ntot, Ntot);
        puts("\n");
#endif

    fft_D(field, dmisc_tmp_z, 3, 2, pars);
    fill_field(dmisc_tmp_z + Ntot, test_func_D2, pars);
    puts("test second derivative in z direction:");
    if (are_fields_equal(dmisc_tmp_z, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testDzz");
        print_vector(dmisc_tmp_z, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_z + Ntot, Ntot);
        puts("\n");
#endif
}

void test_mk_gradient(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;

    mk_gradient_squared(field, dmisc_tmp_x, pars);
    fill_field(dmisc_tmp_x + Ntot, test_func_gradsq, pars);
    puts("test mk gradient squared:");
    if (are_fields_equal(dmisc_tmp_x, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testgradsq");
        print_vector(dmisc_tmp_x, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_x + Ntot, Ntot);
        puts("\n");
#endif
}

void test_mk_laplacian(parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;

    mk_laplacian(field, dmisc_tmp_x, pars);
    fill_field(dmisc_tmp_x + Ntot, test_func_lap, pars);
    puts("test mk laplacian:");
    if (are_fields_equal(dmisc_tmp_x, pars) == 0)
    {
        puts("passed\n");
    }
    else
    {
        puts("failed\n");
    }
#ifdef DEBUG
        puts("testlap");
        print_vector(dmisc_tmp_x, Ntot);
        puts("exact");
        print_vector(dmisc_tmp_x + Ntot, Ntot);
        puts("\n");
#endif
}

void test_fft_apply_filter(parameters_t *pars) {
    //todo
}

double test_func(double x, double y, double z) {
    return sin(x) * sin(y) * sin(z);
}

double test_func_gradsq(double x, double y, double z) {
    double grad_x = test_func_Dx(x,y,z);
    // double grad_y = test_func_Dy(x,y,z);
    double grad_z = test_func_Dz(x,y,z);
    return  grad_x * grad_x + grad_z * grad_z; // + grad_y * grad_y
}

double test_func_lap(double x, double y, double z) {
    return -3.0 * sin(x) * sin(y) * sin(z);
}

double test_func_Dx(double x, double y, double z) {
    return cos(x) * sin(y) * sin(z);
}

double test_func_Dy(double x, double y, double z) {
    return sin(x) * cos(y) * sin(z);
}

double test_func_Dz(double x, double y, double z) {
    return sin(x) * sin(y) * cos(z);
}

double test_func_D2(double x, double y, double z) {
    return -sin(x) * sin(y) * sin(z);
}

void fill_field(double *f, double (*func)(double, double, double),
                parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t osx, osy;
    double x, y, z;

    for (size_t i = 0; i < Nx; ++i)
    {
        x = grid[i];
        osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j)
        {
            y = grid[Nx + j];
            osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k)
            {
                z = grid[Nx + Ny + k];
                f[osy + k] = func(x, y, z);
            }
        }
    }
}

int are_fields_equal(double *f, parameters_t *pars) {
    size_t Nx = pars->x.N;
    size_t Ny = pars->y.N;
    size_t Nz = pars->z.N;
    size_t Ntot = Nx * Ny * Nz;

    for (size_t i = 0; i < Ntot; ++i)
    {
        if (equal(f[i], f[Ntot + i]) != 0)
        {
            return -1;
        }
    }
    return 0;
}

int equal(double a, double b) {
    return fabs(a - b) < 1e-8 ? 0 : -1;
}