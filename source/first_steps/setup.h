#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>
#include "main.h"

void allocate_and_initialize_all();
void initialize_threading();
void initialize_parameters();
void allocate_external();
void mk_grid();
void mk_fftw_plans();
void mk_initial_conditions();
double phi_init(double x, double y, double z, double *phases);
double dphi_init(double x, double y, double z);
void free_and_destroy_all();
void destroy_fftw_plans();
void free_all_external();
void print_matrix(double *matrix, size_t N);
void print_vector(double *vector, size_t N);

#endif