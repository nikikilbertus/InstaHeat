#ifndef __EVOLUTION_TOOLKIT__
#define __EVOLUTION_TOOLKIT__

#include <stddef.h>

void fft_D1(double *f, double *result, size_t N);
void fft_D2(double *f, double *result, size_t N);
void fft_apply_filter(double *f, size_t N);
void mk_filter_window(double *window, size_t nmax, size_t nc);
extern double filter_window_function(double x);
#ifdef USE_ACCELERATE_FRAMEWORK
#ifdef SPECTRAL_OPERATOR_DERIVATIVE
void spectral_op_D2(double *f, double *result, size_t N);
#endif
void matrix_matrix(double *matrixA, double *matrixB, double *result, size_t N);
void sq_matrix(double *matrix, double *result, size_t N);
void matrix_vector(double *matrix, double *vector, double *result, size_t N);
#endif
void print_matrix(double *matrix, size_t N);
void print_vector(double *vector, size_t N);
size_t idx(size_t N, size_t row, size_t col);

#endif