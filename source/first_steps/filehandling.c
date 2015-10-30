#include <stdio.h>
#include <stdlib.h>
#include "main.h"

void print_vector_to_file(double *v, size_t N, int filenumber) {
	if (filenumber < 0 || filenumber > 999)
	{
		fputs("Need an odd number of gridpoints for a fourier grid", stderr);
		return;
	}

	char filename[sizeof "file000.txt"];
    sprintf(filename, "file%03d.txt", filenumber);

    FILE *phi_f;
    phi_f = fopen(filename, "w");

    for (size_t i = 0; i < TS; i+=10)
    {
    	for (size_t j = 0; j < N; j+=2)
    	{
    		fprintf(phi_f, "%.10f, ", v[i*N+j]);
    	}
    	fputs("\n", phi_f);
    }
    fclose(phi_f);
}