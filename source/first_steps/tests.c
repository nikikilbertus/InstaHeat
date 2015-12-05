#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "main.h"
#include "setup.h"
#include "evolution_toolkit.h"

void run_all_tests() {
    test_mk_gradient_squared_and_laplacian();
    // test_fft_apply_filter();
}

void test_mk_gradient_squared_and_laplacian() {
    size_t N = pars.N;

    mk_gradient_squared_and_laplacian(field, dtmp_grad2, dtmp_lap);
    puts("test mk gradient squared and laplacian:");
    fill_field(dtmp_x + N, test_func_Dx);
    if (are_fields_equal(dtmp_x, dtmp_x + N) == 0)
    {
        puts("Dx passed\n");
    }
    else
    {
        puts("Dx failed\n");
    }
    fill_field(dtmp_x + N, test_func_Dy);
    if (are_fields_equal(dtmp_y, dtmp_x + N) == 0)
    {
        puts("Dy passed\n");
    }
    else
    {
        puts("Dy failed\n");
    }
    fill_field(dtmp_x + N, test_func_Dz);
    if (are_fields_equal(dtmp_z, dtmp_x + N) == 0)
    {
        puts("Dz passed\n");
    }
    else
    {
        puts("Dz failed\n");
    }

    fill_field(dtmp_x, test_func_gradsq);
    fill_field(dtmp_y, test_func_lap);
    if (are_fields_equal(dtmp_x, dtmp_grad2) == 0)
    {
        puts("gradient squared passed\n");
    }
    else
    {
        puts("gradient squared failed\n");
    }
    if (are_fields_equal(dtmp_y, dtmp_lap) == 0)
    {
        puts("laplace passed\n");
    }
    else
    {
        puts("laplace failed\n");
    }
#ifdef DEBUG
        puts("testgradsq");
        print_vector(dtmp_grad2, N);
        puts("exact");
        print_vector(dtmp_x, N);
        puts("\n");
        puts("testlap");
        print_vector(dtmp_lap, N);
        puts("exact");
        print_vector(dtmp_y, N);
        puts("\n");
#endif
}

void test_fft_apply_filter() {
    //todo
}

double test_func(double x, double y, double z) {
    return sin(x) * sin(y) * sin(z);
}

double test_func_gradsq(double x, double y, double z) {
    double grad_x = test_func_Dx(x,y,z);
    double grad_y = test_func_Dy(x,y,z);
    double grad_z = test_func_Dz(x,y,z);
    return grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
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

void fill_field(double *f, double (*func)(double, double, double)) {
    size_t Nx = pars.x.N;
    size_t Ny = pars.y.N;
    size_t Nz = pars.z.N;
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

int are_fields_equal(double *f, double *g) {
    for (size_t i = 0; i < pars.N; ++i)
    {
        if (equal(f[i], g[i]) != 0)
        {
            return -1;
        }
    }
    return 0;
}

int equal(double a, double b) {
    return fabs(a - b) < 1e-8 ? 0 : -1;
}