#ifndef __SETUP__
#define __SETUP__

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
void free_and_destroy_all();

#endif
