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

    // initial setup of dimensions
    hsize_t dim[2] = {0, N};
    hsize_t max[2] = {H5S_UNLIMITED, N};
    hsize_t chunk[2] = {Nt, N};

    // ---------------------------full fields: phi, psi, rho--------------------
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.field), "phi");
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.dfield), "dphi");
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.field), "psi");
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.dfield), "dpsi");
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_rho.field), "rho");

    // ---------------------------power spectrum--------------------------------
    dim[1] = bins;
    max[1] = bins;
    chunk[1] = bins;
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_powspec),
            "power_spectrum");

    // ---------------------------time, a, means, variances---------------------
    rank = 1;
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_time), "time");
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_a), "a");

    // ---------------------------parameters------------------------------------
    double val[3] = {MASS, 0.0, 0.0};
    h5_write_parameter(file, "mass", val, 1);

    val[0] = pars.dim;
    h5_write_parameter(file, "dimension", val, 1);

    val[0] = SEED;
    h5_write_parameter(file, "seed", val, 1);

    val[0] = pars.file.skip;
    h5_write_parameter(file, "strides_time", val, 1);

    val[0] = RELATIVE_TOLERANCE;
    val[1] = ABSOLUTE_TOLERANCE;
    h5_write_parameter(file, "tolerances", val, 2);

    val[0] = pars.x.N;
    val[1] = pars.y.N;
    val[2] = pars.z.N;
    h5_write_parameter(file, "gridpoints_internal", val, 3);

    val[0] = pars.x.outN;
    val[1] = pars.y.outN;
    val[2] = pars.z.outN;
    h5_write_parameter(file, "gridpoints_output", val, 3);

    val[0] = pars.x.stride;
    val[1] = pars.y.stride;
    val[2] = pars.z.stride;
    h5_write_parameter(file, "strides_space", val, 3);

    // ---------------------------commit hash-----------------------------------
    #if VERSION_CONTROL != VERSION_CONTROL_NONE
    #if VERSION_CONTROL == VERSION_CONTROL_HG
    char *cmd = "hg id -i";
    #elif VERSION_CONTROL == VERSION_CONTROL_GIT
    char *cmd = "git rev-parse --short HEAD";
    #endif
    size_t len = 16;
    char hash[len];
    FILE *output;

    if ((output = popen(cmd, "r")) == NULL)
    {
        fputs("Could not get hash of current commit.\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (fgets(hash, len, output) != NULL)
    {
        hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, len - 1);
        hid_t memtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(memtype, len);
        hid_t dspace_str = H5Screate_simple (1, dim, NULL);

        hid_t dset_str = H5Dcreate(file, "commit-hash", filetype, dspace_str,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_str, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, hash);
        H5Dclose(dset_str);
        H5Sclose(dspace_str);
        H5Tclose(filetype);
        H5Tclose(memtype);
    }
    else
    {
        fputs("Could not parse hash of current commit.\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (pclose(output))
    {
        fputs("Could not close file of commit hash.\n", stderr);
        exit(EXIT_FAILURE);
    }
    #endif

    RUNTIME_INFO(puts("Created hdf5 file with parameters and datasets for "
                "phi, dphi, psi, dpsi, t, a, rho, powerspec.\n"));
}

void h5_create_dset(const hsize_t rank, const hsize_t *dim,
        const hsize_t *max, const hsize_t *chunk, hsize_t *dset,
        const char *name) {
    // create dataspace
    hid_t dspace = H5Screate_simple(rank, dim, max);

    // create property list
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    H5Pset_chunk(plist, rank, chunk);

    // create and save dataset
    hid_t dset_id = H5Dcreate(pars.file.id, name, H5T_NATIVE_DOUBLE,
                            dspace, H5P_DEFAULT, plist, H5P_DEFAULT);
    (*dset) = dset_id;

    // close property list and dataspace
    H5Pclose(plist);
    H5Sclose(dspace);
}

void h5_write_parameter(const hid_t file, const char *name, const double *val,
        size_t N) {
    hsize_t rank = 1;
    hsize_t dim[1] = {N};
    hsize_t max[1] = {N};

    hid_t dspace_par = H5Screate_simple(rank, dim, max);
    hid_t plist_par = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist_par, H5D_COMPACT);
    hid_t dset_par = H5Dcreate(file, name, H5T_NATIVE_DOUBLE,
                            dspace_par, H5P_DEFAULT, plist_par, H5P_DEFAULT);
    H5Dwrite(dset_par, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
    H5Pclose(plist_par);
    H5Dclose(dset_par);
    H5Sclose(dspace_par);
}

void h5_get_extent(const hid_t dset, hsize_t *max, hsize_t *cur) {
    hid_t dspace = H5Dget_space(dset);
    H5Sget_simple_extent_dims(dspace, cur, max);
}

void h5_write_buffer(const hsize_t rank, const hsize_t *start,
        const hsize_t *add, const hsize_t *new_dim, hsize_t dset,
        const double *buf) {
    hid_t mem_space = H5Screate_simple(rank, add, NULL);
    hid_t dspace = H5Dget_space(dset);
    H5Dset_extent(dset, new_dim);
    //TODO: necessary to call this again?
    dspace = H5Dget_space(dset);
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, add, NULL);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, dspace, H5P_DEFAULT, buf);
    H5Sclose(mem_space);
    H5Sclose(dspace);
}

void h5_write_all_buffers(const hsize_t Nt) {
    hsize_t N = pars.outN;
    hsize_t bins = pars.file.bins_powspec;
    hsize_t rank;
    file_parameters_t f = pars.file;
    // TODO[performance] maybe use static variable to count dataset size instead
    // of reading it from the file each time
    // static hsize_t counter;

    #ifdef SHOW_TIMING_INFO
    h5_time_write -= get_wall_time();
    #endif

    rank = 2;
    hsize_t add[2] = {Nt, N};
    hsize_t curr_dim[rank];
    hsize_t max[rank];
    hsize_t new_dim[rank];
    hsize_t start[rank];

    h5_get_extent(pars.file.dset_phi.field, max, curr_dim);
    new_dim[0] = curr_dim[0] + Nt;
    new_dim[1] = N;
    start[0] = curr_dim[0];
    start[1] = 0;

    // ---------------------------full fields: phi, psi, rho--------------------
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.field, phi_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.dfield, dphi_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.field, psi_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.dfield, dpsi_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_rho.field, rho_buf);

    // --------------------------power spectrum---------------------------------
    add[1] = bins;
    new_dim[1] = bins;
    h5_write_buffer(rank, start, add, new_dim, f.dset_powspec, pow_spec_buf);

    // ---------------------------time, a, means, variances---------------------
    rank = 1;
    h5_write_buffer(rank, start, add, new_dim, f.dset_time, time_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_a, f_a_buf);

    #ifdef SHOW_TIMING_INFO
    h5_time_write += get_wall_time();
    #endif
}

void h5_close() {
    hid_t file = pars.file.id;
    if (pars.file.index != 0)
    {
        h5_write_all_buffers(pars.file.index);
    }
    H5Fflush(file, H5F_SCOPE_GLOBAL);
    hid_t obj_ids[30];
    hsize_t obj_count = H5Fget_obj_ids(file, H5F_OBJ_DATASET, -1, obj_ids);
    for (size_t i = 0; i < obj_count; ++i)
    {
        H5Dclose(obj_ids[i]);
    }
    H5Fclose(file);
}

void save() {
    hsize_t index = pars.file.index;
    hsize_t Nt    = pars.file.buf_size;
    hsize_t Nx    = pars.x.N;
    hsize_t Ny    = pars.y.N;
    hsize_t Nz    = pars.z.N;
    hsize_t outy  = pars.y.outN;
    hsize_t outz  = pars.z.outN;
    hsize_t N     = pars.N;
    hsize_t outN  = pars.outN;
    hsize_t bins  = pars.file.bins_powspec;

    hsize_t os = index * outN;
    size_t osx, osy, id;
    size_t osxb, osyb, idb;
    #pragma omp parallel for private(osx, osxb, osy, osyb, id, idb)
    for (size_t i = 0; i < Nx; i += pars.x.stride)
    {
        osx = i * Ny * Nz;
        osxb = i * outy * outz / pars.x.stride;
        for (size_t j = 0; j < Ny; j += pars.y.stride)
        {
            osy = osx + j * Nz;
            osyb = osxb + j * outz / pars.y.stride;
            for (size_t k = 0; k < Nz; k += pars.z.stride)
            {
                id = osy + k;
                idb = osyb + k / pars.z.stride;
                phi_buf[os + idb]  = field[id];
                dphi_buf[os + idb] = field[N + id];
                psi_buf[os + idb]  = psi[id];
                dpsi_buf[os + idb] = dpsi[id];
                rho_buf[os + idb]  = rho[id];
                #ifdef CHECK_FOR_NAN
                if (isnan(field[id]) || isnan(psi[id]) || isnan(rho[id]))
                {
                    fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
                    exit(EXIT_FAILURE);
                }
                #endif
            }
        }
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
        h5_write_all_buffers(Nt);
        pars.file.index = 0;
    }
    else
    {
        pars.file.index += 1;
    }
    RUNTIME_INFO(printf("Writing to file at t = %f\n", pars.t.t));
}
