#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>
#include "main.h"

void allocate_and_initialize_all(parameters_t *pars);
void initialize_parameters(parameters_t *pars);
void allocate_external(parameters_t *pars);
void mk_grid(parameters_t *pars);
void mk_fftw_plans(parameters_t *pars);
void mk_initial_conditions(parameters_t *pars);
double phi_init(double x, double y, double z);
double dphi_init(double x, double y, double z);
void free_and_destroy_all(parameters_t *pars);
void destroy_fftw_plans();
void free_all_external(parameters_t *pars);
void print_matrix(double *matrix, size_t N);
void print_vector(double *vector, size_t N);

#endif