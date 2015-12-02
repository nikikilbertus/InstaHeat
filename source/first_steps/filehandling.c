#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#include "filehandling.h"
#include "main.h"

void h5_create_empty_by_path(char *name) {
    hsize_t rank = 2;
    // TODO DELETE AFTER TESTING
    hsize_t N = 10;//pars.Ntot;
    herr_t status;

    // create file
    hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    pars.file.id = file;

    // --------------------------phi------------------------------------------
    // create dataspace for the phi
    hsize_t dims[2] = {WRITE_OUT_BUFFER_NUMBER, N};
    hsize_t max_dims[2] = {H5S_UNLIMITED, N};
    hid_t dspace_phi = H5Screate_simple(rank, dims, max_dims);

    // create property list for the phi
    hid_t plist_phi = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_phi, H5D_CHUNKED);
    hsize_t chunk_dims[2] = {WRITE_OUT_BUFFER_NUMBER, N};
    H5Pset_chunk(plist_phi, rank, chunk_dims);

    // create dataset for the phi
    hid_t dset_phi = H5Dcreate(file, "phi", H5T_NATIVE_DOUBLE,
                            dspace_phi, H5P_DEFAULT, plist_phi, H5P_DEFAULT);
    pars.file.dset_phi = dset_phi;

    // close property list and dataspace
    status = H5Pclose(plist_phi);
    status = H5Sclose(dspace_phi);

    // --------------------------time-------------------------------------------
    // create dataspace for the time
    rank = 1;
    dims[0] = WRITE_OUT_BUFFER_NUMBER;
    max_dims[0] = H5S_UNLIMITED;
    hid_t dspace_time = H5Screate_simple(rank, dims, max_dims);

    // create property list for the time
    hid_t plist_time = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_time, H5D_CHUNKED);
    chunk_dims[0] = WRITE_OUT_BUFFER_NUMBER;
    H5Pset_chunk(plist_time, rank, chunk_dims);

    // create dataset for the time
    hid_t dset_time = H5Dcreate(file, "time", H5T_NATIVE_DOUBLE,
                            dspace_time, H5P_DEFAULT, plist_time, H5P_DEFAULT);
    pars.file.dset_time = dset_time;

    // close property list and dspace_time
    status = H5Pclose(plist_time);
    status = H5Sclose(dspace_time);

    // ---------------------------a---------------------------------------------
    // create dataspace for the a
    hid_t dspace_a = H5Screate_simple(rank, dims, max_dims);

    // create property list for the a
    hid_t plist_a = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_a, H5D_CHUNKED);
    H5Pset_chunk(plist_a, rank, chunk_dims);

    // create dataset for the a
    hid_t dset_a = H5Dcreate(file, "a", H5T_NATIVE_DOUBLE,
                            dspace_a, H5P_DEFAULT, plist_a, H5P_DEFAULT);
    pars.file.dset_a = dset_a;

    // close property list and dspace_a
    status = H5Pclose(plist_a);
    status = H5Sclose(dspace_a);

    // ---------------------------rho-------------------------------------------
    // create dataspace for the rho
    hid_t dspace_rho = H5Screate_simple(rank, dims, max_dims);

    // create property list for the rho
    hid_t plist_rho = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_rho, H5D_CHUNKED);
    H5Pset_chunk(plist_rho, rank, chunk_dims);

    // create dataset for the rho
    hid_t dset_rho = H5Dcreate(file, "rho", H5T_NATIVE_DOUBLE,
                            dspace_rho, H5P_DEFAULT, plist_rho, H5P_DEFAULT);
    pars.file.dset_rho = dset_rho;

    // close property list and dspace_rho
    status = H5Pclose(plist_rho);
    status = H5Sclose(dspace_rho);
}

void h5_append_dataset_by_id_2d(hid_t dset, double *field,
                                    hsize_t N, hsize_t Nt) {
    const hsize_t rank = 2;
    herr_t status;

    hsize_t add_dims[2] = {Nt, N};
    hid_t mem_space = H5Screate_simple(rank, add_dims, NULL);

    hid_t dataspace = H5Dget_space(dset);
    hsize_t curr_dims[rank];
    hsize_t max_dims[rank];
    H5Sget_simple_extent_dims(dataspace, curr_dims, max_dims);

    hsize_t new_dims[2] = {curr_dims[0] + Nt, N};
    H5Dset_extent(dset, new_dims);

    hsize_t start_dims[2] = {curr_dims[0] + 1, 0};
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dataspace, H5P_DEFAULT, field);

    H5Sclose(mem_space);
    H5Sclose(dataspace);
}

void h5_close(hid_t file) {
    H5Fflush(file, H5F_SCOPE_GLOBAL);
    hid_t obj_ids[10];
    hsize_t obj_count = H5Fget_obj_ids(file, H5F_OBJ_DATASET, -1, obj_ids);
    for (hsize_t i = 0; i < obj_count; ++i)
    {
        H5Dclose(obj_ids[i]);
    }
    H5Fclose(file);
}

/*
shortcut for writing a single vector to a new file.
the filepath is composed, the file is created/emptied and the data is written.
row_skip: only print every row_skip-th entry
*/
void file_single_write_by_name_1d(double *data, size_t N, size_t row_skip,
                            char *name) {
    RUNTIME_INFO(puts("Start writing to file..."));
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[] = ".txt";
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_create_empty_by_path(filepath);
    file_append_by_path_1d(data, N, row_skip, filepath);
    RUNTIME_INFO(puts("Finished writing to file."));
}

/*
create empty file <filepath>
*/
void file_create_empty_by_path(char *filepath) {
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
void file_create_empty_by_name(char *name) {
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[] = ".txt";
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_create_empty_by_path(filepath);
}

/*
append a N \times 1 vector to a file with the given name, the path is determined
row_skip: only print every row_skip-th entry
*/
void file_append_by_name_1d(double *data, size_t N, size_t row_skip,
                                char *name) {
    int length = strlen(DATAPATH) + 32;
    char filepath[length];
    char suffix[] = ".txt";
    strcpy(filepath, DATAPATH);
    strcat(filepath, name);
    strcat(filepath, suffix);
    file_append_by_path_1d(data, N, row_skip, filepath);
}

/*
append a N \times 1 vector to a file with the path <filepath>.
row_skip: only print every row_skip-th entry
*/
void file_append_by_path_1d(double *data, size_t N, size_t row_skip,
                                char *filepath) {
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
deprecated
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