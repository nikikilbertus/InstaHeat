#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "hdf5.h"

void h5_create_empty_by_path(char *name);
void h5_write_buffers_to_disk(hsize_t Nt);
void h5_close();
void save();

#endif
