#ifndef __RKF45__
#define __RKF45__

#include "main.h"

/**
 * @file RKF45.h
 * @brief Function declarations for RKF45.c
 */

void run_rkf45();
int mk_rhs_wrapper(double t, const double f[], double result[], void * params);

#endif
