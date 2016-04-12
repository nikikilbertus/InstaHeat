#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "hdf5.h"
#include "filehandling.h"
#include "main.h"

void h5_create_empty_by_path(const char *name)
{
    const hsize_t Nt = pars.file.buf_size;
    hsize_t rank;

    // create file
    const hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    pars.file.id = file;

    rank = 2;
    // ---------------------------full fields: phi, dphi, psi, dpsi, rho--------
    #ifdef OUTPUT_PHI
    h5_create_dset(rank, Nt, phi.dim, &(phi.id), H5_PHI_NAME);
    #endif
    #ifdef OUTPUT_DPHI
    h5_create_dset(rank, Nt, dphi.dim, &(dphi.id), H5_DPHI_NAME);
    #endif
    #ifdef OUTPUT_PSI
    h5_create_dset(rank, Nt, psi.dim, &(psi.id), H5_PSI_NAME);
    #endif
    #ifdef OUTPUT_DPSI
    h5_create_dset(rank, Nt, dpsi.dim, &(dpsi.id), H5_DPSI_NAME);
    #endif
    #ifdef OUTPUT_RHO
    h5_create_dset(rank, Nt, rho_out.dim, &(rho_out.id), H5_RHO_NAME);
    #endif

    // ---------------------------summaries-------------------------------------
    #ifdef OUTPUT_PHI_SMRY
    h5_create_dset(rank, Nt, phi_smry.dim, &(phi_smry.id), H5_PHI_SMRY_NAME);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    h5_create_dset(rank, Nt, dphi_smry.dim, &(dphi_smry.id), H5_DPHI_SMRY_NAME);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    h5_create_dset(rank, Nt, psi_smry.dim, &(psi_smry.id), H5_PSI_SMRY_NAME);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    h5_create_dset(rank, Nt, dpsi_smry.dim, &(dpsi_smry.id), H5_DPSI_SMRY_NAME);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    h5_create_dset(rank, Nt, rho_smry.dim, &(rho_smry.id), H5_RHO_SMRY_NAME);
    #endif

    // ---------------------------power spectra---------------------------------
    #ifdef OUTPUT_PHI_PS
    h5_create_dset(rank, Nt, phi_ps.dim, &(phi_ps.id), H5_PHI_PS_NAME);
    #endif

    rank = 1;
    // ---------------------------time, a---------------------------------------
    h5_create_dset(rank, Nt, time.dim, &(time.id), H5_TIME_NAME);
    h5_create_dset(rank, Nt, a_out.dim, &(a_out.id), H5_A_NAME);

    // ---------------------------parameters------------------------------------
    double val[3] = {MASS, 0.0, 0.0};
    h5_write_parameter(H5_MASS_NAME, val, 1);

    val[0] = pars.dim;
    h5_write_parameter(H5_DIMENSION_NAME, val, 1);

    val[0] = SEED;
    h5_write_parameter(H5_SEED_NAME, val, 1);

    val[0] = pars.file.skip;
    h5_write_parameter(H5_STRIDES_TIME_NAME, val, 1);

    val[0] = RELATIVE_TOLERANCE;
    val[1] = ABSOLUTE_TOLERANCE;
    h5_write_parameter(H5_TOLERANCES_NAME, val, 2);

    val[0] = pars.x.N;
    val[1] = pars.y.N;
    val[2] = pars.z.N;
    h5_write_parameter(H5_GRIDPOINTS_INTERNAL_NAME, val, 3);

    val[0] = SPATIAL_LOWER_BOUND_X;
    val[1] = SPATIAL_UPPER_BOUND_X;
    h5_write_parameter(H5_SPATIAL_BOUNDS_X_NAME, val, 2);

    val[0] = SPATIAL_LOWER_BOUND_Y;
    val[1] = SPATIAL_UPPER_BOUND_Y;
    h5_write_parameter(H5_SPATIAL_BOUNDS_Y_NAME, val, 2);

    val[0] = SPATIAL_LOWER_BOUND_Z;
    val[1] = SPATIAL_UPPER_BOUND_Z;
    h5_write_parameter(H5_SPATIAL_BOUNDS_Z_NAME, val, 2);

    val[0] = pars.x.outN;
    val[1] = pars.y.outN;
    val[2] = pars.z.outN;
    h5_write_parameter(H5_GRIDPOINTS_OUTPUT_NAME, val, 3);

    val[0] = pars.x.stride;
    val[1] = pars.y.stride;
    val[2] = pars.z.stride;
    h5_write_parameter(H5_STRIDES_SPACE_NAME, val, 3);

    #ifdef MAX_DT_HUBBLE_FRACTION
    val[0] = MAX_DT_HUBBLE_FRACTION;
    #else
    val[0] = -1.0;
    #endif
    h5_write_parameter(H5_MAX_DT_HUBBLE_FRACTION_NAME, val, 1);

    // ---------------------------commit hash-----------------------------------
    hid_t filetype, memtype, dspace_str, dset_str;
    size_t len;
    #if VERSION_CONTROL != VERSION_CONTROL_NONE
    #if VERSION_CONTROL == VERSION_CONTROL_HG
    const char *cmd = "hg id -i";
    #elif VERSION_CONTROL == VERSION_CONTROL_GIT
    const char *cmd = "git rev-parse --short HEAD";
    #endif
    hsize_t dim[1] = {1};
    len = 16;
    char hash[len];
    FILE *output;

    if ((output = popen(cmd, "r")) == NULL) {
        fputs("Could not get hash of current commit.\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (fgets(hash, len, output) != NULL) {
        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, len - 1);
        memtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(memtype, len);
        dspace_str = H5Screate_simple (1, dim, NULL);

        dset_str = H5Dcreate(file, H5_COMMIT_HASH_NAME, filetype,
                dspace_str, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_str, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, hash);
        H5Dclose(dset_str);
        H5Sclose(dspace_str);
        H5Tclose(filetype);
        H5Tclose(memtype);
    } else {
        fputs("Could not parse hash of current commit.\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (pclose(output)) {
        fputs("Could not close file of commit hash.\n", stderr);
        exit(EXIT_FAILURE);
    }
    #endif

    //-----------------write out psi method-------------------------------------
    #if PSI_METHOD == PSI_ELLIPTIC
    len = 9;
    const char *method = "elliptic";
    #elif PSI_METHOD == PSI_PARABOLIC
    len = 10;
    const char *method = "parabolic";
    #elif PSI_METHOD == PSI_HYPERBOLIC
    len = 11;
    const char *method = "hyperbolic";
    #endif
    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, len - 1);
    memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, len);
    dspace_str = H5Screate_simple (1, dim, NULL);

    dset_str = H5Dcreate(file, H5_METHOD_NAME, filetype,
            dspace_str, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_str, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, method);
    H5Dclose(dset_str);
    H5Sclose(dspace_str);
    H5Tclose(filetype);
    H5Tclose(memtype);

    INFO(puts("Created hdf5 file with datasets for output.\n"));
}

void h5_create_dset(const hsize_t rank, const hsize_t Nt, const hsize_t N,
        hsize_t *dset, const char *name)
{
    hsize_t dim[2] = {0, N};
    hsize_t max[2] = {H5S_UNLIMITED, N};
    hsize_t chunk[2] = {Nt, N};
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

void h5_write_parameter(const char *name, const double *val, const size_t N)
{
    hsize_t rank = 1;
    hsize_t dim[1] = {N};
    hsize_t max[1] = {N};
    hid_t file = pars.file.id;

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

void h5_get_extent(hsize_t *max, hsize_t *cur)
{
    hid_t dspace = H5Dget_space(time.id);
    H5Sget_simple_extent_dims(dspace, cur, max);
}

void h5_write_all_buffers(const hsize_t Nt)
{
    // TODO[performance] maybe use static variable to count dataset size instead
    // of reading it from the file each time
    // static hsize_t counter;
    #ifdef SHOW_TIMING_INFO
    h5_time_write -= get_wall_time();
    #endif

    hsize_t rank;
    hsize_t curr_dim[1];
    hsize_t max[1];
    h5_get_extent(max, curr_dim);
    const hsize_t os = curr_dim[0];

    rank = 2;
    // ---------------------------full fields: phi, psi, rho--------------------
    #ifdef OUTPUT_PHI
    h5_write_buffer(rank, Nt, phi.dim, os, phi.id, phi.buf);
    #endif
    #ifdef OUTPUT_DPHI
    h5_write_buffer(rank, Nt, dphi.dim, os, dphi.id, dphi.buf);
    #endif
    #ifdef OUTPUT_PSI
    h5_write_buffer(rank, Nt, psi.dim, os, psi.id, psi.buf);
    #endif
    #ifdef OUTPUT_DPSI
    h5_write_buffer(rank, Nt, dpsi.dim, os, dpsi.id, dpsi.buf);
    #endif
    #ifdef OUTPUT_RHO
    h5_write_buffer(rank, Nt, rho.dim, os, rho_out.id, rho_out.buf);
    #endif

    // --------------------------power spectra----------------------------------
    #ifdef OUTPUT_PHI_PS
    h5_write_buffer(rank, Nt, phi_ps.dim, os, phi_ps.id, phi_ps.buf);
    #endif

    // --------------------------summaries--------------------------------------
    #ifdef OUTPUT_PHI_SMRY
    h5_write_buffer(rank, Nt, phi_smry.dim, os, phi_smry.id, phi_smry.buf);
    #endif
    #ifdef OUTPUT_DPHI_SMRY
    h5_write_buffer(rank, Nt, dphi_smry.dim, os, dphi_smry.id, dphi_smry.buf);
    #endif
    #ifdef OUTPUT_PSI_SMRY
    h5_write_buffer(rank, Nt, psi_smry.dim, os, psi_smry.id, psi_smry.buf);
    #endif
    #ifdef OUTPUT_DPSI_SMRY
    h5_write_buffer(rank, Nt, dpsi_smry.dim, os, dpsi_smry.id, dpsi_smry.buf);
    #endif
    #ifdef OUTPUT_RHO_SMRY
    h5_write_buffer(rank, Nt, rho_smry.dim, os, rho_smry.id, rho_smry.buf);
    #endif

    rank = 1;
    // ---------------------------time and  a-----------------------------------
    h5_write_buffer(rank, Nt, time.dim, os, time.id, time.buf);
    h5_write_buffer(rank, Nt, a_out.dim, os, a_out.id, a_out.buf);

    #ifdef SHOW_TIMING_INFO
    h5_time_write += get_wall_time();
    #endif
    INFO(printf("Dumping to disk at t = %f\n", pars.t.t));
}

void h5_write_buffer(const hsize_t rank, const hsize_t Nt,
        const hsize_t N, const hsize_t os, const hsize_t dset,
        const double *buf)
{
    hsize_t add[2] = {Nt, N};
    hsize_t new_dim[2] = {os + Nt, N};
    hsize_t start[2] = {os, 0};
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

void h5_close()
{
    hid_t file = pars.file.id;
    if (pars.file.index != 0) {
        h5_write_all_buffers(pars.file.index);
    }
    H5Fflush(file, H5F_SCOPE_GLOBAL);
    hid_t obj_ids[30];
    hsize_t obj_count = H5Fget_obj_ids(file, H5F_OBJ_DATASET, -1, obj_ids);
    for (size_t i = 0; i < obj_count; ++i) {
        H5Dclose(obj_ids[i]);
    }
    H5Fclose(file);
}

void save()
{
    hsize_t index = pars.file.index;
    hsize_t Nt = pars.file.buf_size;
    hsize_t N = pars.N;
    hsize_t N2 = 2 * N;
    hsize_t bins = pars.file.bins_powspec;

    time.buf[index] = pars.t.t;
    a_buf[index] = field[N2];

    #ifdef CHECK_FOR_NAN
    if (isnan(pars.t.t) || isnan(field[N2])) {
        fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
    }
    #endif

    #ifdef OUTPUT_PHI_MEAN
    phi_mean_buf[index] = phi_mean;
    #endif
    #ifdef OUTPUT_PHI_VARIANCE
    phi_var_buf[index] = phi_var;
    #endif
    #ifdef OUTPUT_DPHI_MEAN
    dphi_mean_buf[index] = dphi_mean;
    #endif
    #ifdef OUTPUT_DPHI_VARIANCE
    dphi_var_buf[index] = dphi_var;
    #endif
    #ifdef OUTPUT_PSI_MEAN
    psi_mean_buf[index] = psi_mean;
    #endif
    #ifdef OUTPUT_PSI_VARIANCE
    psi_var_buf[index] = psi_var;
    #endif
    #ifdef OUTPUT_DPSI_MEAN
    dpsi_mean_buf[index] = dpsi_mean;
    #endif
    #ifdef OUTPUT_DPSI_VARIANCE
    dpsi_var_buf[index] = dpsi_var;
    #endif
    #ifdef OUTPUT_RHO_MEAN
    rho_mean_buf[index] = rho_mean;
    #endif
    #ifdef OUTPUT_RHO_VARIANCE
    rho_var_buf[index] = rho_var;
    #endif

    #ifdef LARGE_OUTPUT
    hsize_t N2p = N2 + 2;
    hsize_t N3p = 3 * N + 2;
    hsize_t Nx = pars.x.N;
    hsize_t Ny = pars.y.N;
    hsize_t Nz = pars.z.N;
    hsize_t outy = pars.y.outN;
    hsize_t outz = pars.z.outN;
    hsize_t outN = pars.outN;
    hsize_t os = index * outN;
    size_t osx, osy, id;
    size_t osxb, osyb, idb;
    #pragma omp parallel for private(osx, osxb, osy, osyb, id, idb)
    for (size_t i = 0; i < Nx; i += pars.x.stride) {
        osx = i * Ny * Nz;
        osxb = i * outy * outz / pars.x.stride;
        for (size_t j = 0; j < Ny; j += pars.y.stride) {
            osy = osx + j * Nz;
            osyb = osxb + j * outz / pars.y.stride;
            for (size_t k = 0; k < Nz; k += pars.z.stride) {
                id = osy + k;
                idb = osyb + k / pars.z.stride;
                #ifdef OUTPUT_PHI
                phi_buf[os + idb] = field[id];
                #endif
                #ifdef OUTPUT_DPHI
                dphi_buf[os + idb] = field[N + id];
                #endif
                #ifdef OUTPUT_PSI
                psi_buf[os + idb] = field[N2p + id];
                #endif
                #ifdef OUTPUT_DPSI
                    #if PSI_METHOD == PSI_PARABOLIC
                    dpsi_buf[os + idb] = dfield[N2p + id];
                    #else
                    dpsi_buf[os + idb] = field[N3p + id];
                    #endif
                #endif
                #ifdef OUTPUT_RHO
                rho_buf[os + idb] = rho[id];
                #endif
                #ifdef CHECK_FOR_NAN
                if (isnan(field[id]) || isnan(rho[id])) {
                    fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
                    exit(EXIT_FAILURE);
                }
                #endif
            }
        }
    }
    #endif

    #ifdef OUTPUT_POWER_SPECTRUM
    hsize_t os1 = index * bins;
    #pragma omp parallel for
    for (size_t i = 0; i < bins; ++i) {
        pow_spec_buf[os1 + i] = pow_spec[i];
        #ifdef CHECK_FOR_NAN
        if (isnan(pow_spec[i])) {
            fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
        }
        #endif
    }
    #endif

    if (index == Nt - 1) {
        h5_write_all_buffers(Nt);
        pars.file.index = 0;
    } else {
        pars.file.index += 1;
    }
    #ifdef DEBUG
    printf("Writing to file at t = %f\n", pars.t.t);
    #endif
}

void h5_read_timeslice()
{
    hid_t file, dset, dspace;
    size_t N = pars.N;
    size_t N2p = 2 * N + 2;
    size_t N3p = 3 * N + 2;
    double t = pars.t.ti;

    file = H5Fopen(INITIAL_DATAPATH, H5F_ACC_RDONLY, H5P_DEFAULT);

    // ---------------------------get time and find index-----------------------
    dset = H5Dopen(file, H5_TIME_NAME, H5P_DEFAULT);
    dspace = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    if (ndims != 1) {
        INFO(fputs("Could not read time properly.\n", stderr));
        exit(EXIT_FAILURE);
    }
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    size_t Nt = dims[0];

    double *time_tmp = calloc(Nt, sizeof *time_tmp);
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_tmp);

    hsize_t index = Nt;
    for (size_t i = 0; i < Nt; ++i) {
        if (time_tmp[i] + DBL_EPSILON >= t) {
            index = i;
            break;
        }
    }
    if (index == Nt) {
        INFO(puts("The initial time is larger than the maximal time in"
                    " the h5 file. Starting at last existing timeslice."));
        index = Nt - 1;
    }
    pars.t.ti = time_tmp[index];
    pars.t.t = time_tmp[index];

    free(time_tmp);
    H5Dclose(dset);
    H5Sclose(dspace);

    // ---------------------------read fields at index--------------------------
    h5_read_and_fill(file, index, H5_PHI_NAME, field);
    h5_read_and_fill(file, index, H5_DPHI_NAME, field + N);
    h5_read_and_fill(file, index, H5_PSI_NAME, field + N2p);
    h5_read_and_fill(file, index, H5_DPSI_NAME, field + N3p);

    // ---------------------------read a at index-------------------------------
    dset = H5Dopen(file, H5_A_NAME, H5P_DEFAULT);
    dspace = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_ndims(dspace);
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    if (dims[0] != Nt) {
        INFO(fputs("Dimensions of dataset does not agree with specified"
                    " values.\n", stderr));
        exit(EXIT_FAILURE);
    }
    hsize_t dimsm[1] = {1};
    hid_t mspace = H5Screate_simple(1, dimsm, NULL);
    hsize_t start_m[1] = {0};
    hsize_t count_m[1] = {1};
    H5Sselect_hyperslab(mspace, H5S_SELECT_SET, start_m, NULL, count_m, NULL);
    hsize_t start[1] = {index};
    hsize_t count[1] = {1};
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, field + 2 * N);
    H5Dclose(dset);
    H5Sclose(dspace);
    H5Fclose(file);
}

void h5_read_and_fill(const hid_t file, const hsize_t index, const char *name,
        double *out)
{
    size_t N = pars.N;
    hid_t dset = H5Dopen(file, name, H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    if (dims[1] != N || ndims != 2) {
        INFO(fputs("Dimensions of dataset does not agree with specified"
                    " values.\n", stderr));
        exit(EXIT_FAILURE);
    }
    hsize_t dimsm[1] = {N};
    hid_t mspace = H5Screate_simple(1, dimsm, NULL);
    hsize_t start_m[1] = {0};
    hsize_t count_m[1] = {N};
    H5Sselect_hyperslab(mspace, H5S_SELECT_SET, start_m, NULL, count_m, NULL);
    hsize_t start[2] = {index, 0};
    hsize_t count[2] = {1, N};
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, out);
    H5Dclose(dset);
    H5Sclose(dspace);
    H5Sclose(mspace);
}

void read_initial_data()
{
    size_t N = pars.N;
    size_t N2p = 2 * N + 2;
    size_t N3p = 3 * N + 2;

    FILE *file = fopen(INITIAL_DATAPATH, "r");
    if (!file) {
        fputs("Could not read initial data file.\n", stderr);
        exit(EXIT_FAILURE);
    }

    int ii, jj, kk;
    for (size_t i = 0; i < N; ++i) {
        if(!fscanf(file, " %d %d %d %lf %lf %lf %lf\n",
                    &ii, &jj, &kk, &field[i], &field[i + N], &field[i + N2p],
                    &field[i + N3p])) {
            fputs("Could not read initial data file.\n", stderr);
            exit(EXIT_FAILURE);
        }
    }
    fclose(file);
}
