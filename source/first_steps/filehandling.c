#include <stdio.h>
#include <stdlib.h>
#include <float.h>
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
    #ifdef OUTPUT_PHI
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.field),
            H5_PHI_NAME);
    #endif
    #ifdef OUTPUT_DPHI
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.dfield),
            H5_DPHI_NAME);
    #endif
    #ifdef OUTPUT_PSI
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.field),
            H5_PSI_NAME);
    #endif
    #ifdef OUTPUT_DPSI
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.dfield),
            H5_DPSI_NAME);
    #endif
    #ifdef OUTPUT_RHO
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_rho.field),
            H5_RHO_NAME);
    #endif

    // ---------------------------power spectrum--------------------------------
    #ifdef OUTPUT_POWER_SPECTRUM
    dim[1] = bins;
    max[1] = bins;
    chunk[1] = bins;
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_powspec),
            H5_POWER_SPECTRUM_NAME);
    #endif

    // ---------------------------time, a, means, variances---------------------
    rank = 1;
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_time), H5_TIME_NAME);
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_a), H5_A_NAME);
    #ifdef OUTPUT_PHI_MEAN
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.mean),
            H5_PHI_MEAN_NAME);
    #endif
    #ifdef OUTPUT_PHI_VARIANCE
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.var),
            H5_PHI_VARIANCE_NAME);
    #endif
    #ifdef OUTPUT_DPHI_MEAN
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.dmean),
            H5_DPHI_MEAN_NAME);
    #endif
    #ifdef OUTPUT_DPHI_VARIANCE
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_phi.dvar),
            H5_DPHI_VARIANCE_NAME);
    #endif
    #ifdef OUTPUT_PSI_MEAN
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.mean),
            H5_PSI_MEAN_NAME);
    #endif
    #ifdef OUTPUT_PSI_VARIANCE
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.var),
            H5_PSI_VARIANCE_NAME);
    #endif
    #ifdef OUTPUT_DPSI_MEAN
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.dmean),
            H5_DPSI_MEAN_NAME);
    #endif
    #ifdef OUTPUT_DPSI_VARIANCE
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_psi.dvar),
            H5_DPSI_VARIANCE_NAME);
    #endif
    #ifdef OUTPUT_RHO_MEAN
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_rho.mean),
            H5_RHO_MEAN_NAME);
    #endif
    #ifdef OUTPUT_RHO_VARIANCE
    h5_create_dset(rank, dim, max, chunk, &(pars.file.dset_rho.var),
            H5_RHO_VARIANCE_NAME);
    #endif

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

    val[0] = pars.x.outN;
    val[1] = pars.y.outN;
    val[2] = pars.z.outN;
    h5_write_parameter(H5_GRIDPOINTS_OUTPUT_NAME, val, 3);

    val[0] = pars.x.stride;
    val[1] = pars.y.stride;
    val[2] = pars.z.stride;
    h5_write_parameter(H5_STRIDES_SPACE_NAME, val, 3);

    // ---------------------------commit hash-----------------------------------
    #if VERSION_CONTROL != VERSION_CONTROL_NONE
    #if VERSION_CONTROL == VERSION_CONTROL_HG
    char *cmd = "hg id -i";
    #elif VERSION_CONTROL == VERSION_CONTROL_GIT
    char *cmd = "git rev-parse --short HEAD";
    #endif
    dim[0] = 1;
    size_t len = 16;
    char hash[len];
    FILE *output;

    hid_t filetype, memtype, dspace_str, dset_str;
    if ((output = popen(cmd, "r")) == NULL)
    {
        fputs("Could not get hash of current commit.\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (fgets(hash, len, output) != NULL)
    {
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

    //-----------------write out psi method-------------------------------------
    #if PSI_METHOD == PSI_ELLIPTIC
    len = 9;
    char *method = "elliptic";
    #elif PSI_METHOD == PSI_PARABOLIC
    len = 10;
    char *method = "parabolic";
    #elif PSI_METHOD == PSI_HYPERBOLIC
    len = 11;
    char *method = "hyperbolic";
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

    RUNTIME_INFO(puts("Created hdf5 file with datasets for output.\n"));
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

void h5_write_parameter(const char *name, const double *val, size_t N) {
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

void h5_get_extent(hsize_t *max, hsize_t *cur) {
    hid_t dspace = H5Dget_space(pars.file.dset_time);
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
    hsize_t curr_dim[1];
    hsize_t max[2];
    h5_get_extent(max, curr_dim);

    hsize_t add[2] = {Nt, N};
    hsize_t new_dim[2] = {curr_dim[0] + Nt, N};
    hsize_t start[2] = {curr_dim[0], 0};
    max[1] = N;

    // ---------------------------full fields: phi, psi, rho--------------------
    #ifdef OUTPUT_PHI
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.field, phi_buf);
    #endif
    #ifdef OUTPUT_DPHI
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.dfield, dphi_buf);
    #endif
    #ifdef OUTPUT_PSI
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.field, psi_buf);
    #endif
    #ifdef OUTPUT_DPSI
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.dfield, dpsi_buf);
    #endif
    #ifdef OUTPUT_RHO
    h5_write_buffer(rank, start, add, new_dim, f.dset_rho.field, rho_buf);
    #endif

    // --------------------------power spectrum---------------------------------
    #ifdef OUTPUT_POWER_SPECTRUM
    add[1] = bins;
    new_dim[1] = bins;
    h5_write_buffer(rank, start, add, new_dim, f.dset_powspec, pow_spec_buf);
    #endif

    // ---------------------------time, a, means, variances---------------------
    rank = 1;
    h5_write_buffer(rank, start, add, new_dim, f.dset_time, time_buf);
    h5_write_buffer(rank, start, add, new_dim, f.dset_a, a_buf);
    #ifdef OUTPUT_PHI_MEAN
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.mean, phi_mean_buf);
    #endif
    #ifdef OUTPUT_PHI_VARIANCE
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.var, phi_var_buf);
    #endif
    #ifdef OUTPUT_DPHI_MEAN
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.dmean, dphi_mean_buf);
    #endif
    #ifdef OUTPUT_DPHI_VARIANCE
    h5_write_buffer(rank, start, add, new_dim, f.dset_phi.dvar, dphi_var_buf);
    #endif
    #ifdef OUTPUT_PSI_MEAN
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.mean, psi_mean_buf);
    #endif
    #ifdef OUTPUT_PSI_VARIANCE
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.var, psi_var_buf);
    #endif
    #ifdef OUTPUT_DPSI_MEAN
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.dmean, dpsi_mean_buf);
    #endif
    #ifdef OUTPUT_DPSI_VARIANCE
    h5_write_buffer(rank, start, add, new_dim, f.dset_psi.dvar, dpsi_var_buf);
    #endif
    #ifdef OUTPUT_RHO_MEAN
    h5_write_buffer(rank, start, add, new_dim, f.dset_rho.mean, rho_mean_buf);
    #endif
    #ifdef OUTPUT_RHO_VARIANCE
    h5_write_buffer(rank, start, add, new_dim, f.dset_rho.var, rho_var_buf);
    #endif

    #ifdef SHOW_TIMING_INFO
    h5_time_write += get_wall_time();
    #endif
    RUNTIME_INFO(printf("Dumping to disk at t = %f\n", pars.t.t));
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
    hsize_t N     = pars.N;
    hsize_t N2    = 2 * N;
    hsize_t N2p   = N2 + 2;
    hsize_t N3p   = 3 * N + 2;
    hsize_t bins  = pars.file.bins_powspec;

    time_buf[index] = pars.t.t;
    a_buf[index] = field[N2];

    #ifdef CHECK_FOR_NAN
    if (isnan(pars.t.t) || isnan(field[N2]))
    {
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

    #if defined(OUTPUT_PHI) || defined(OUTPUT_DPHI) || defined(OUTPUT_PSI) || \
        defined(OUTPUT_DPSI) || defined(OUTPUT_RHO)
    hsize_t Nx    = pars.x.N;
    hsize_t Ny    = pars.y.N;
    hsize_t Nz    = pars.z.N;
    hsize_t outy  = pars.y.outN;
    hsize_t outz  = pars.z.outN;
    hsize_t outN  = pars.outN;
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
                #ifdef OUTPUT_PHI
                phi_buf[os + idb]  = field[id];
                #endif
                #ifdef OUTPUT_DPHI
                dphi_buf[os + idb] = field[N + id];
                #endif
                #ifdef OUTPUT_PSI
                psi_buf[os + idb]  = field[N2p + id];
                #endif
                #ifdef OUTPUT_DPSI
                    #if PSI_METHOD == PSI_PARABOLIC
                    dpsi_buf[os + idb] = dfield[N2p + id];
                    #else
                    dpsi_buf[os + idb] = field[N3p + id];
                    #endif
                #endif
                #ifdef OUTPUT_RHO
                rho_buf[os + idb]  = rho[id];
                #endif
                #ifdef CHECK_FOR_NAN
                if (isnan(field[id]) || isnan(rho[id]))
                {
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
    for (size_t i = 0; i < bins; ++i)
    {
        pow_spec_buf[os1 + i] = pow_spec[i];
        #ifdef CHECK_FOR_NAN
        if (isnan(pow_spec[i]))
        {
            fprintf(stderr, "Discovered nan at time: %f \n", pars.t.t);
            exit(EXIT_FAILURE);
        }
        #endif
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
    #ifdef DEBUG
    printf("Writing to file at t = %f\n", pars.t.t);
    #endif
}

void h5_read_timeslice() {
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
    if (ndims != 1)
    {
        RUNTIME_INFO(fputs("Could not read time properly.\n", stderr));
        exit(EXIT_FAILURE);
    }
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    size_t Nt = dims[0];

    double *time_tmp = calloc(Nt, sizeof *time_tmp);
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_tmp);

    hsize_t index = Nt;
    for (size_t i = 0; i < Nt; ++i)
    {
        if (time_tmp[i] + DBL_EPSILON >= t)
        {
            index = i;
            break;
        }
    }
    if (index == Nt)
    {
        RUNTIME_INFO(puts("The initial time is larger than the maximal time in"
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
    if (dims[0] != Nt)
    {
        RUNTIME_INFO(fputs("Dimensions of dataset does not agree with specified"
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
        double *out) {
    size_t N = pars.N;
    hid_t dset = H5Dopen(file, name, H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    if (dims[1] != N || ndims != 2)
    {
        RUNTIME_INFO(fputs("Dimensions of dataset does not agree with specified"
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

void read_initial_data() {
    size_t N = pars.N;

    FILE *file = fopen(INITIAL_DATAPATH, "r");
    if (!file)
    {
        fputs("Could not read initial data file.\n", stderr);
        exit(EXIT_FAILURE);
    }

    //TODO: adjust to actual file format, this is just a dummy
    for (size_t i = 0; i < N; ++i)
    {
        fscanf(file, "%lf", &field[i]);
    }
    fclose(file);
}
