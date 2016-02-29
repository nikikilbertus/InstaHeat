#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "hdf5.h"

void h5_create_empty_by_path(const char *name);
void h5_create_dset(const hsize_t rank, const hsize_t *dim,
        const hsize_t *max_dim, const hsize_t *chunk, hsize_t *dset_id,
        const char *name);
void h5_write_parameter(const char *name, const double *val, size_t N);
void h5_write_all_buffers(const hsize_t Nt);
void h5_get_extent(hsize_t *max, hsize_t *cur);
void h5_close();
void save();

#endif
