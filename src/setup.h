#ifndef __SETUP__
#define __SETUP__

#include "main.h"

/**
 * @file setup.h
 *
 * @brief Function declarations for `setup.c`
 *
 * @note Both of the two publicly visible function
 * `allocate_and_initialize_all()` and `free_and_destroy_all()` are each called
 * exactly once before and after the simulation respectively. Setup and cleanup
 * is usually fast compared to the actual integration, so we do not have to
 * worry about performance here.
 */

void allocate_and_initialize_all();
void free_and_destroy_all();

#endif
