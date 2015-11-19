#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "filehandling.h"
#include "main.h"

/*
shortcut for writing a single vector to a new file.
the filepath is composed, the file is created/emptied and the data is written.
row_skip: only print every row_skip-th entry
*/
void file_single_write_1d(double *data, size_t N, size_t row_skip,
                            char *name, int number) {
    RUNTIME_INFO(puts("Start writing to file..."));
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[16];
    sprintf(suffix, "_%03d.txt", number);
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_create_empty(filepath);
    file_append_1d(data, N, row_skip, filepath);
    RUNTIME_INFO(puts("Finished writing to file."));
}

/*
create empty file <filepath>
*/
void file_create_empty(char *filepath) {
    if (filepath == NULL)
    {
        fputs("Filepath is NULL.", stderr);
        exit(EXIT_FAILURE);
    }
    RUNTIME_INFO(puts("Create empty file..."));
    FILE *empty;
    empty = fopen(filepath, "w");
    if (!empty)
    {
        fputs("Could not open file.", stderr);
        exit(EXIT_FAILURE);
    }
    fclose(empty);
    RUNTIME_INFO(puts("File created."));
}

/*
create empty file by name, path is determined
*/
void file_create_empty_by_name(char *name, int number) {
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[16];
    sprintf(suffix, "_%03d.txt", number);
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_create_empty(filepath);
}

/*
append a N \times 1 vector to a file with the given name, the path is determined
row_skip: only print every row_skip-th entry
*/
void file_append_by_name_1d(double *data, size_t N, size_t row_skip,
                                char *name, int number) {
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[16];
    sprintf(suffix, "_%03d.txt", number);
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_append_1d(data, N, row_skip, filepath);
}

/*
append a N \times 1 vector to a file with the path <filepath>.
row_skip: only print every row_skip-th entry
*/
void file_append_1d(double *data, size_t N, size_t row_skip, char *filepath) {
	if (row_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		exit(EXIT_FAILURE);
	}
    if (filepath == NULL)
    {
        fputs("Filepath is NULL.", stderr);
        exit(EXIT_FAILURE);
    }

    FILE *vector_f;
    vector_f = fopen(filepath, "a");
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
}

/*
append a N \times M matrix to a file with the path <filepath>.
row_skip and col_skip: only every <>-th row/column is printed
*/
void file_append_2d(double *data, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *filepath) {
	if (row_skip < 1 || col_skip < 1)
	{
		fputs("Skip amount has to be at least 1.", stderr);
		return;
	}
    if (filepath == NULL)
    {
        fputs("Filepath is NULL.", stderr);
        exit(EXIT_FAILURE);
    }

    FILE *matrix_f;
    matrix_f = fopen(filepath, "a");
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
}