#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "hdf5.h"
#include "filehandling.h"
#include "main.h"

void h5_create_empty_by_path(const char *name) {
    hsize_t rank = 2;
    hsize_t N = pars.outN;
    hsize_t Nt = pars.file.buf_size;
    hsize_t bins = pars.file.bins_powspec;

    // create file
    hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    pars.file.id = file;

    // --------------------------phi--------------------------------------------
    // create dataspace for phi
    hsize_t dims[2] = {0, N};
    hsize_t max_dims[2] = {H5S_UNLIMITED, N};
    hid_t dspace_phi = H5Screate_simple(rank, dims, max_dims);

    // create property list for phi
    hid_t plist_phi = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_phi, H5D_CHUNKED);
    hsize_t chunk_dims[2] = {Nt, N};
    H5Pset_chunk(plist_phi, rank, chunk_dims);

    // create dataset for phi
    hid_t dset_phi = H5Dcreate(file, "phi", H5T_NATIVE_DOUBLE,
                            dspace_phi, H5P_DEFAULT, plist_phi, H5P_DEFAULT);
    pars.file.dset_phi = dset_phi;

    // close property list and dataspace
    H5Pclose(plist_phi);
    H5Sclose(dspace_phi);

    // --------------------------psi--------------------------------------------
    // create dataspace for psi
    hid_t dspace_psi = H5Screate_simple(rank, dims, max_dims);

    // create property list for psi
    hid_t plist_psi = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_psi, H5D_CHUNKED);
    H5Pset_chunk(plist_psi, rank, chunk_dims);

    // create dataset for psi
    hid_t dset_psi = H5Dcreate(file, "psi", H5T_NATIVE_DOUBLE,
                            dspace_psi, H5P_DEFAULT, plist_psi, H5P_DEFAULT);
    pars.file.dset_psi = dset_psi;

    // close property list and dataspace
    H5Pclose(plist_psi);
    H5Sclose(dspace_psi);

    // --------------------------rho--------------------------------------------
    // create dataspace for rho
    hid_t dspace_rho = H5Screate_simple(rank, dims, max_dims);

    // create property list for rho
    hid_t plist_rho = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_rho, H5D_CHUNKED);
    H5Pset_chunk(plist_rho, rank, chunk_dims);

    // create dataset for rho
    hid_t dset_rho = H5Dcreate(file, "rho", H5T_NATIVE_DOUBLE,
                            dspace_rho, H5P_DEFAULT, plist_rho, H5P_DEFAULT);
    pars.file.dset_rho = dset_rho;

    // close property list and dataspace
    H5Pclose(plist_rho);
    H5Sclose(dspace_rho);

    // --------------------------power spectrum---------------------------------
    // create dataspace for the power_spectrum
    dims[1] = bins;
    max_dims[1] = bins;
    hid_t dspace_powspec = H5Screate_simple(rank, dims, max_dims);

    // create property list for the power_spectrum
    hid_t plist_powspec = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_powspec, H5D_CHUNKED);
    chunk_dims[1] = bins;
    H5Pset_chunk(plist_powspec, rank, chunk_dims);

    // create dataset for the power_spectrum
    hid_t dset_powspec = H5Dcreate(file, "power_spectrum", H5T_NATIVE_DOUBLE,
                    dspace_powspec, H5P_DEFAULT, plist_powspec, H5P_DEFAULT);
    pars.file.dset_powspec = dset_powspec;

    // close property list and dataspace
    H5Pclose(plist_powspec);
    H5Sclose(dspace_powspec);

    // --------------------------time-------------------------------------------
    // create dataspace for the time
    rank = 1;
    hid_t dspace_time = H5Screate_simple(rank, dims, max_dims);

    // create property list for the time
    hid_t plist_time = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_time, H5D_CHUNKED);
    H5Pset_chunk(plist_time, rank, chunk_dims);

    // create dataset for the time
    hid_t dset_time = H5Dcreate(file, "time", H5T_NATIVE_DOUBLE,
                            dspace_time, H5P_DEFAULT, plist_time, H5P_DEFAULT);
    pars.file.dset_time = dset_time;

    // close property list and dspace_time
    H5Pclose(plist_time);
    H5Sclose(dspace_time);

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
    H5Pclose(plist_a);
    H5Sclose(dspace_a);

    RUNTIME_INFO(puts("Created hdf5 file with datasets for "
                "phi, psi, t, a, rho.\n"));
}

void h5_write_buffers_to_disk(const hsize_t Nt) {
    hsize_t N = pars.outN;
    hsize_t bins = pars.file.bins_powspec;
    hsize_t rank;
    // TODO[performance] maybe use static variable to count dataset size instead
    // of reading it from the file each time
    // static hsize_t counter;

    #ifdef SHOW_TIMING_INFO
    h5_time_write -= get_wall_time();
    #endif

    // --------------------------phi--------------------------------------------
    rank = 2;
    hid_t dset = pars.file.dset_phi;

    hsize_t add_dims[2] = {Nt, N};
    hid_t mem_space = H5Screate_simple(rank, add_dims, NULL);
    hid_t dspace = H5Dget_space(dset);
    hsize_t curr_dims[rank];
    hsize_t max_dims[rank];
    H5Sget_simple_extent_dims(dspace, curr_dims, max_dims);
    hsize_t new_dims[2] = {curr_dims[0] + Nt, N};
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    hsize_t start_dims[2] = {curr_dims[0], 0};
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT,
                    field_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    // --------------------------psi--------------------------------------------
    dset = pars.file.dset_psi;

    mem_space = H5Screate_simple(rank, add_dims, NULL);
    dspace = H5Dget_space(dset);
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT,
                    psi_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    // --------------------------rho--------------------------------------------
    dset = pars.file.dset_rho;

    mem_space = H5Screate_simple(rank, add_dims, NULL);
    dspace = H5Dget_space(dset);
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT,
                    rho_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    // --------------------------power spectrum---------------------------------
    dset = pars.file.dset_powspec;

    add_dims[1] = bins;
    mem_space = H5Screate_simple(rank, add_dims, NULL);
    dspace = H5Dget_space(dset);
    new_dims[1] = bins;
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT,
                    pow_spec_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    // --------------------------time-------------------------------------------
    rank = 1;
    dset = pars.file.dset_time;

    mem_space = H5Screate_simple(rank, add_dims, NULL);
    dspace = H5Dget_space(dset);
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT, time_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    // --------------------------a----------------------------------------------
    dset = pars.file.dset_a;

    mem_space = H5Screate_simple(rank, add_dims, NULL);
    dspace = H5Dget_space(dset);
    H5Dset_extent(dset, new_dims);
    dspace = H5Dget_space(dset);

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start_dims, NULL,
                            add_dims, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT, f_a_buf);

    H5Sclose(mem_space);
    H5Sclose(dspace);

    #ifdef SHOW_TIMING_INFO
    h5_time_write += get_wall_time();
    #endif
}

void h5_close() {
    hid_t file = pars.file.id;
    if (pars.file.index != 0)
    {
        h5_write_buffers_to_disk(pars.file.index);
    }
    H5Fflush(file, H5F_SCOPE_GLOBAL);
    hid_t obj_ids[10];
    hsize_t obj_count = H5Fget_obj_ids(file, H5F_OBJ_DATASET, -1, obj_ids);
    for (size_t i = 0; i < obj_count; ++i)
    {
        H5Dclose(obj_ids[i]);
    }
    H5Fclose(file);
}

void save() {
    hsize_t index = pars.file.index;
    hsize_t Nt = pars.file.buf_size;
    hsize_t N = pars.N;
    hsize_t bins = pars.file.bins_powspec;

    hsize_t os = index * N;
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        field_buf[os + i] = field[i];
        psi_buf[os + i] = psi[i];
        rho_buf[os + i] = rho[i];
        #ifdef CHECK_FOR_NAN
        if (isnan(field[i]) || isnan(psi[i]) || isnan(rho[i]))
        {
            fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
        }
        #endif
    }

    os = index * bins;
    #pragma omp parallel for
    for (size_t i = 0; i < bins; ++i)
    {
        pow_spec_buf[os + i] = pow_spec[i];
        #ifdef CHECK_FOR_NAN
        if (isnan(pow_spec[i]))
        {
            fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
        }
        #endif
    }

    time_buf[index] = pars.t.t;
    f_a_buf[index] = field[2 * N];

    #ifdef CHECK_FOR_NAN
    if (isnan(pars.t.t) || isnan(field[2 * N]))
    {
        fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
    }
    #endif

    if (index == Nt - 1)
    {
        h5_write_buffers_to_disk(Nt);
        pars.file.index = 0;
    }
    else
    {
        pars.file.index += 1;
    }
    RUNTIME_INFO(printf("Writing to file at t = %f\n", pars.t.t));
}
