#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>
#include "main.h"
#include "setup.h"
#include "RK4_stepper.h"
#include "tests.h"

//grid and spectral operators
double *x;
double *D1;
double *D2;

// timestep counter and number of timesteps
size_t nt;
size_t TS;

//evolution of phi and spatial derivatives (NOGP * TSMAX space)
double *phi;
double *phiD1;
double *phiD2;

//solution for dsigma
double *dphi;
double *dphiD1;
double *dphiD2;

int main(int argc, const char * argv[]) {

#ifdef RUN_TESTS_ONLY
    run_all_tests();
    return 0;
#endif

    allocate_and_initialize_all();

    run_RK4_stepper(DT, NOGP);

    free_all_external();
    return 0;
}
