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
void mk_x_grid();
void mk_fftw_plans();
void mk_initial_conditions();
void mk_k_grid();
void mk_filter_mask();
extern double filter_window(const double x);
void read_initial_data();
double phi_init(const double x, const double y, const double z,
        const double *ph);
double dphi_init(const double x, const double y, const double z,
        const double *ph);
double wrapped_gaussian(const double x, const double y, const double z);
void mk_initial_psi();
void free_and_destroy_all();
void destroy_and_cleanup_fftw();
void free_external();
void mk_bunch_davies(double *f, const double H, const double homo,
        const double gamma);
extern complex box_muller();
void print_vector(const double *vector, const size_t N);

#endif
