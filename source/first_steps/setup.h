#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>
#include "main.h"

// everything here is only called once in the beginning or the end, thus there
// is no need for any optimization
void allocate_and_initialize_all();
void initialize_threading();
void initialize_parameters();
void allocate_external();
void mk_grid();
void mk_fftw_plans();
void mk_initial_conditions();
double phi_init(const double x, const double y, const double z,
        const double *phases);
double dphi_init(const double x, const double y, const double z);
void free_and_destroy_all();
void destroy_and_cleanup_fftw();
void free_external();
void print_vector(const double *vector, const size_t N);

#endif
