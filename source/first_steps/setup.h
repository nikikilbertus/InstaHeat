#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>

void allocate_and_initialize_all(void);
void allocate_external(size_t N);
void mk_fourier_spectral_operators(size_t N, double low_bnd, double up_bnd);
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