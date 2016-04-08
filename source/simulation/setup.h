#ifndef __SETUP__
#define __SETUP__

#include <stddef.h>
#include "main.h"

/**
 * @file setup.h
 * @brief Function declarations for setup.c
 * @note The only two functions called from the outside of this file are
 * allocate_and_initialize_all() and free_and_destroy_all() before and after the
 * simulation respectively. Moreover they are only called once. Setup and
 * cleanup is usually fast, so we do not care about performance here.
 */

void allocate_and_initialize_all();
void initialize_threading();
void initialize_parameters();
void allocate_external();
void mk_fftw_plans();
void mk_k_grid();
void mk_filter_mask();
extern double filter_window(const double x);
void mk_initial_conditions();
void initialize_from_dat();
void initialize_from_bunch_davies();
void mk_bunch_davies(double *f, const double H, const double homo,
        const double gamma);
extern complex box_muller();
void initialize_from_internal_function();
void mk_x_grid(double *grid);
double phi_init(const double x, const double y, const double z,
        const double *ph);
double dphi_init(const double x, const double y, const double z,
        const double *ph);
double wrapped_gaussian(const double x, const double y, const double z);
void mk_initial_psi();
void free_and_destroy_all();
void destroy_and_cleanup_fftw();
void free_external();
void print_vector(const double *vector, const size_t N);

#endif
