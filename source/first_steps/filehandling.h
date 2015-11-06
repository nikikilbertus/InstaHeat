#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include <stddef.h>

void print_vector_to_file(double *v, size_t N, size_t row_skip,
							char *prefix, int filenumber);
void print_matrix_to_file(double *A, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *prefix, int filenumber);

#endif