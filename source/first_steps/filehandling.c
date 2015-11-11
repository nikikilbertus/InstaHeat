#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"

/*
print a N \times 1 vector to a .txt file with the name "<prefix>_<filenumber>".
row_skip: only print every row_skip-th entry
*/
void print_vector_to_file(double *v, size_t N, size_t row_skip,
							char *prefix, int filenumber) {
	if (filenumber < 0 || filenumber > 999)
	{
		fputs("Filenumber should not exceed three digits.", stderr);
		exit(EXIT_FAILURE);
	}
	if (row_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		exit(EXIT_FAILURE);
	}
    if (prefix == NULL)
    {
        prefix = "vector";
    }

    RUNTIME_INFO(printf("Start writing %s to file...\n", prefix));

    size_t prefix_length = strlen(prefix);
	char filename[prefix_length + (sizeof "../../../data/_000.txt")];
    sprintf(filename, "../../../data/%s_%03d.txt", prefix, filenumber);

    FILE *vector_f;
    vector_f = fopen(filename, "w");

    if (!vector_f)
    {
    	fputs("Could not open file.", stderr);
		exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N; i+=row_skip)
    {
    	fprintf(vector_f, "%.15f\n", v[i]);
    }
    fclose(vector_f);
    RUNTIME_INFO(puts("Finished writing to file.\n"));
}

/*
print a N \times M matrix to a .txt file with the name "<prefix>_<filenumber>".
row_skip and col_skip: only every <>-th row/column is printed
*/
void print_matrix_to_file(double *A, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *prefix, int filenumber) {
	if (filenumber < 0 || filenumber > 999)
	{
		fputs("Filenumber should not exceed three digits.", stderr);
		return;
	}
	if (row_skip < 1 || col_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		return;
	}
    if (prefix == NULL)
    {
        prefix = "matrix";
    }

    RUNTIME_INFO(printf("Start writing %s to file...\n", prefix));

    size_t prefix_length = strlen(prefix);
    char filename[prefix_length + (sizeof "../../../data/_000.txt")];
    sprintf(filename, "../../../data/%s_%03d.txt", prefix, filenumber);

    FILE *matrix_f;
    matrix_f = fopen(filename, "w");

    if (!matrix_f)
    {
    	fputs("Could not open file.", stderr);
		exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N; i+=row_skip)
    {
    	for (size_t j = 0; j < M; j+=col_skip)
    	{
    		fprintf(matrix_f, "%.15f, ", A[i*N+j]);
    	}
    	fputs("\n", matrix_f);
    }
    fclose(matrix_f);
    RUNTIME_INFO(puts("Finished writing to file.\n"));
}