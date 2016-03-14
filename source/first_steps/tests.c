#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "main.h"
#include "setup.h"
#include "evolution_toolkit.h"

void run_all_tests()
{
    test_mk_gradient_squared_and_laplacian();
}

void test_mk_gradient_squared_and_laplacian()
{
    const size_t N = pars.N;

    fill_field(field, test_func);
    mk_gradient_squared_and_laplacian(field);

    puts("test mk gradient squared and laplacian:");
    fill_field(tmp.xphi + N, test_func_Dx);
    if (are_fields_equal(tmp.xphi, tmp.xphi + N) == 0) {
        puts("Dx passed\n");
    }
    else {
        puts("Dx failed\n");
    }
    fill_field(tmp.xphi + N, test_func_Dy);
    if (are_fields_equal(tmp.yphi, tmp.xphi + N) == 0) {
        puts("Dy passed\n");
    }
    else {
        puts("Dy failed\n");
    }
    fill_field(tmp.xphi + N, test_func_Dz);
    if (are_fields_equal(tmp.zphi, tmp.xphi + N) == 0) {
        puts("Dz passed\n");
    }
    else {
        puts("Dz failed\n");
    }

    fill_field(tmp.xphi, test_func_gradsq);
    fill_field(tmp.yphi, test_func_lap);
    if (are_fields_equal(tmp.xphi, tmp.grad) == 0) {
        puts("gradient squared passed\n");
    }
    else {
        puts("gradient squared failed\n");
    }
    if (are_fields_equal(tmp.yphi, tmp.lap) == 0) {
        puts("laplace passed\n");
    }
    else {
        puts("laplace failed\n");
    }
    #ifdef DEBUG
    puts("testgradsq");
    print_vector(tmp.grad, N);
    puts("exact");
    print_vector(tmp.xphi, N);
    puts("\n");
    puts("testlap");
    print_vector(tmp.lap, N);
    puts("exact");
    print_vector(tmp.yphi, N);
    puts("\n");
    #endif
}

double test_func_gradsq(const double x, const double y, const double z)
{
    const double grad_x = test_func_Dx(x,y,z);
    const double grad_y = test_func_Dy(x,y,z);
    const double grad_z = test_func_Dz(x,y,z);
    return grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
}

double test_func_lap(const double x, const double y, const double z)
{
    return test_func_D2x(x, y, z) + test_func_D2y(x, y, z) +
        test_func_D2z(x, y, z);
}

double test_func(const double x, const double y, const double z)
{
    /* return exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(2.0 * x) * cos(y) * cos(4.0 * z); */
    return sin(x) * sin(y) * sin(z);
}

double test_func_Dx(const double x, const double y, const double z)
{
    /* return -2.0 * exp(-2. * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(y) * cos(4.0 * z) * (2.0 * x * cos(2.0 * x) + sin(2.0 * x)); */
    return cos(x) * sin(y) * sin(z);
}

double test_func_Dy(const double x, const double y, const double z)
{
    /* return -exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(2.0 * x) * cos(4.0 * z) * (8.0 * y * cos(y) + sin(y)); */
    return sin(x) * cos(y) * sin(z);
}

double test_func_Dz(const double x, const double y, const double z)
{
    /* return exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(2.0 * x) * cos(y) * (-3.0 * z * cos(4.0 * z) - 4.0 * sin(4.0 * z)); */
    return sin(x) * sin(y) * cos(z);
}

double test_func_D2x(const double x, const double y, const double z)
{
    /* return 8.0 * exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(y) * cos(4.0 * z) * ((-1.0 + 2.0 * pow(x, 2)) * cos(2.0 * x) + */
    /*     2.0 * x * sin(2.0 * x)); */
    return -sin(x) * sin(y) * sin(z);
}

double test_func_D2y(const double x, const double y, const double z)
{
    /* return exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(2.0 * x) * cos(4.0 * z) * ((-9.0 + 64.0 * pow(y, 2)) * cos(y) + */
    /*     16.0 * y * sin(y)); */
    return -sin(x) * sin(y) * sin(z);
}

double test_func_D2z(const double x, const double y, const double z)
{
    /* return exp(-2.0 * pow(x, 2) - 4.0 * pow(y, 2) - 1.5 * pow(z, 2)) * */
    /*     cos(2.0 * x) * cos(y) * ((-19.0 + 9.0 * pow(z, 2)) * cos(4.0 * z) + */
    /*     24.0 * z * sin(4.0 * z)); */
    return -sin(x) * sin(y) * sin(z);
}

void fill_field(double *f, double (*func)(const double, const double,
                                                        const double))
{
    const size_t Nx = pars.x.N;
    const size_t Ny = pars.y.N;
    const size_t Nz = pars.z.N;
    size_t osx, osy;
    double x, y, z;

    for (size_t i = 0; i < Nx; ++i) {
        x = grid[i];
        osx = i * Ny * Nz;
        for (size_t j = 0; j < Ny; ++j) {
            y = grid[Nx + j];
            osy = osx + j * Nz;
            for (size_t k = 0; k < Nz; ++k) {
                z = grid[Nx + Ny + k];
                f[osy + k] = func(x, y, z);
            }
        }
    }
}

int are_fields_equal(const double *f, const double *g)
{
    for (size_t i = 0; i < pars.N; ++i) {
        if (equal(f[i], g[i]) != 0) {
            return -1;
        }
    }
    return 0;
}

int equal(const double a, const double b)
{
    return fabs(a - b) < 1e-8 ? 0 : -1;
}
