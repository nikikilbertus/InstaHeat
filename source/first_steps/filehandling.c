#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"

void file_create_empty(char *filename) {
    if (filename == NULL)
    {
        fputs("Filename is NULL.", stderr);
        exit(EXIT_FAILURE);
    }
    FILE *empty;
    empty = fopen(filename, "w");
    if (!empty)
    {
        fputs("Could not open file.", stderr);
        exit(EXIT_FAILURE);
    }
    fclose(empty);
}

/*
print a N \times 1 vector to a .txt file with the name "<prefix>_<filenumber>".
row_skip: only print every row_skip-th entry
*/
void file_append_1d(double *data, size_t N, size_t row_skip, char *filename) {
	if (row_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		exit(EXIT_FAILURE);
	}
    if (filename == NULL)
    {
        fputs("Filename is NULL.", stderr);
        exit(EXIT_FAILURE);
    }

    RUNTIME_INFO(puts("Start writing to file..."));

    FILE *vector_f;
    vector_f = fopen(filename, "a");
    if (!vector_f)
    {
    	fputs("Could not open file.", stderr);
		exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N; i+=row_skip)
    {
    	fprintf(vector_f, "%f\n", data[i]);
    }
    fclose(vector_f);
    RUNTIME_INFO(puts("Finished writing to file.\n"));
}

/*
print a N \times M matrix to a .txt file with the name "<prefix>_<filenumber>".
row_skip and col_skip: only every <>-th row/column is printed
*/
void file_append_2d(double *data, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *filename) {
	if (row_skip < 1 || col_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		return;
	}
    if (filename == NULL)
    {
        fputs("Filename is NULL.", stderr);
        exit(EXIT_FAILURE);
    }

    RUNTIME_INFO(puts("Start writing to file..."));

    FILE *matrix_f;
    matrix_f = fopen(filename, "a");
    if (!matrix_f)
    {
    	fputs("Could not open file.", stderr);
		exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N; i+=row_skip)
    {
    	for (size_t j = 0; j < M; j+=col_skip)
    	{
    		fprintf(matrix_f, "%f, ", data[i*N+j]);
    	}
    	fputs("\n", matrix_f);
    }
    fclose(matrix_f);
    RUNTIME_INFO(puts("Finished writing to file.\n"));
}