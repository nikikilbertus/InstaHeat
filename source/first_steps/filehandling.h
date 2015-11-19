#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include <stddef.h>

void file_single_write_1d(double *data, size_t N, size_t row_skip,
                            char *name, int number);
void file_create_empty(char *filename);
void file_create_empty_by_name(char *name, int number);
void file_append_by_name_1d(double *data, size_t N, size_t row_skip,
                                char *name, int number);
void file_append_1d(double *v, size_t N, size_t row_skip, char *filename);
void file_append_2d(double *A, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *filename);

#endif