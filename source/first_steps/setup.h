#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>
#include "main.h"

void allocate_and_initialize_all(parameters_t *pars);
void initialize_parameters(parameters_t *pars);
void allocate_external(size_t Nx, size_t Nt);
void mk_fourier_spectral_operators(size_t N, double a, double b);
void mk_initial_conditions(size_t N, double (*f_init)(double),
								double (*df_init)(double));
double phi_init(double x);
double dphi_init(double x);
void free_all_external();
void matrix_matrix(double *matrixA, double *matrixB, double *result, size_t N);
void sq_matrix(double *matrix, double *result, size_t N);
void matrix_vector(double *matrix, double *vector, double *result, size_t N);
void print_matrix(double *matrix, size_t N);
void print_vector(double *vector, size_t N);
size_t idx(size_t N, size_t row, size_t col);

#endif