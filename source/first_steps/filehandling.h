#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "hdf5.h"

void h5_create_empty_by_path(const char *name);
void h5_write_parameter(const hid_t file, const char *name, const double *val,
        size_t N);
void h5_write_buffers_to_disk(const hsize_t Nt);
void h5_close();
void save();

#endif
