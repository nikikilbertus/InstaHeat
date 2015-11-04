#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>
#include "main.h"
#include "setup.h"
#include "RK4_stepper.h"
#include "filehandling.h"
#include "tests.h"

// simulation parameters
parameters_t pars;

//grid and spectral operators
double *x;
double *D1;
double *D2;

// evolution of the field and spatial derivatives (2*Nx * Nt space required)
double *field;

// evolution of the scale parameter a for the FRW equations
double *a;

int main(int argc, const char * argv[]) {

#ifdef RUN_TESTS_ONLY
    run_all_tests();
    return 0;
#endif

    int count = 0;
    for (double dt = 0.1; dt > 1e-2; dt /= 2, count++)
    {
    	pars.dt = dt;
    	allocate_and_initialize_all(&pars);
    	run_RK4_stepper(&pars);
		print_vector_to_file(field, 2*pars.Nx*pars.Nt, 1, count);
    	free_all_external();
    }

    return 0;
}
