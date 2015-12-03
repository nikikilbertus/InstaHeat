#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include <stddef.h>
#include "hdf5.h"

void h5_create_empty_by_path(char *name);
void h5_write_buffers_to_disk(hsize_t Nt);
void h5_close();
void save();
void file_single_write_by_name_1d(double *data, size_t N, size_t row_skip,
                            char *name);
void file_create_empty_by_path(char *filepath);
void file_create_empty_by_name(char *name);
void file_append_by_name_1d(double *data, size_t N, size_t row_skip,
                                char *name);
void file_append_by_path_1d(double *v, size_t N, size_t row_skip,
							char *filename);
void file_append_2d(double *A, size_t N, size_t M, size_t row_skip,
                            size_t col_skip, char *filename);

#endif