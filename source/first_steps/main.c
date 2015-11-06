#include <stdio.h>
#include <stdlib.h>
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
double *frw_a;

// evolution of the total energy of the field
double *e_tot;

int main(int argc, const char * argv[]) {

#ifdef RUN_TESTS_ONLY
    pars.dt = 0.1;
    allocate_and_initialize_all(&pars);
    run_all_tests();
    free_all_external();
    return 0;
#endif

    int count = 0;
    for (double dt = 0.1; dt > 1e-2; dt /= 2, count++)
    {
    	pars.dt = dt;
    	allocate_and_initialize_all(&pars);
    	run_RK4_stepper(&pars);
        char *prefix_field = "field";
		print_vector_to_file(field, 2*pars.Nx*pars.Nt, 1, prefix_field, count);
        char *prefix_frw_a = "a";
        print_vector_to_file(frw_a, pars.Nt, 1, prefix_frw_a, count);
        char *prefix_e_tot = "energy";
        print_vector_to_file(e_tot, pars.Nt, 1, prefix_e_tot, count);
    	free_all_external();
    }

    return 0;
}
