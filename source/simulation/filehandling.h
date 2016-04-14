#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "hdf5.h"
#include "main.h"

/**
 * @file filehandling.h
 * @brief Contains preprocessor definitions of the names of the datasets in the
 * output `.h5` file and function declarations for `filehandling.c`.
 */

#define     H5_PHI_NAME                     ("phi")
#define     H5_DPHI_NAME                    ("dphi")
#define     H5_PSI_NAME                     ("psi")
#define     H5_DPSI_NAME                    ("dpsi")
#define     H5_RHO_NAME                     ("rho")
#define     H5_PHI_PS_NAME                  ("phi_power_spectrum")
#define     H5_DPHI_PS_NAME                 ("dphi_power_spectrum")
#define     H5_PSI_PS_NAME                  ("psi_power_spectrum")
#define     H5_DPSI_PS_NAME                 ("dpsi_power_spectrum")
#define     H5_RHO_PS_NAME                  ("rho_power_spectrum")
#define     H5_TIME_NAME                    ("time")
#define     H5_A_NAME                       ("a")
#define     H5_PHI_SMRY_NAME                ("phi_summary")
#define     H5_DPHI_SMRY_NAME               ("dphi_summary")
#define     H5_PSI_SMRY_NAME                ("psi_summary")
#define     H5_DPSI_SMRY_NAME               ("dpsi_summary")
#define     H5_RHO_SMRY_NAME                ("rho_summary")
#define     H5_CONSTRAINTS_NAME             ("constraints")
#define     H5_MASS_NAME                    ("mass")
#define     H5_DIMENSION_NAME               ("dimension")
#define     H5_SEED_NAME                    ("seed")
#define     H5_STRIDES_TIME_NAME            ("strides_time")
#define     H5_STRIDES_SPACE_NAME           ("strides_space")
#define     H5_TOLERANCES_NAME              ("tolerances")
#define     H5_GRIDPOINTS_INTERNAL_NAME     ("gridpoints_internal")
#define     H5_GRIDPOINTS_OUTPUT_NAME       ("gridpoints_output")
#define     H5_SPATIAL_BOUNDS_X_NAME        ("spatial_bounds_x")
#define     H5_SPATIAL_BOUNDS_Y_NAME        ("spatial_bounds_y")
#define     H5_SPATIAL_BOUNDS_Z_NAME        ("spatial_bounds_z")
#define     H5_COMMIT_HASH_NAME             ("commit_hash")
#define     H5_METHOD_NAME                  ("method")
#define     H5_RUNTIME_TOTAL_NAME           ("runtime_total")
#define     H5_RUNTIME_FFTW_NAME            ("runtime_fftw")
#define     H5_RUNTIME_FFTWPLAN_NAME        ("runtime_fftwplan")
#define     H5_RUNTIME_FILTER_NAME          ("runtime_filter")
#define     H5_RUNTIME_ELLIPTIC_NAME        ("runtime_elliptic")
#define     H5_RUNTIME_WRITEOUT_NAME        ("runtime_writeout")
#define     H5_RUNTIME_STEPPER_NAME         ("runtime_stepper")
#define     H5_STEPS_TOTAL_NAME             ("steps_total")
#define     H5_STEPS_OK_NAME                ("steps_ok")
#define     H5_STEPS_BAD_NAME               ("steps_bad")
#define     H5_MAX_DT_HUBBLE_FRACTION_NAME  ("max_dt_hubble_fraction")

void h5_create_empty_by_path();
void h5_create_dset(const hsize_t rank, const hsize_t N, hsize_t *dset,
        const char *name);
void h5_write_parameter(const char *name, const double *val, const size_t N);
void h5_get_extent(hsize_t *cur);
void h5_write_buffer(const hsize_t rank, const hsize_t Nt,
        const hsize_t N, const hsize_t os, const hsize_t dset,
        const double *buf);
void h5_write_all_buffers(const hsize_t Nt);
void h5_close();
void save();
void append_to_buffer(struct output f);
void h5_read_timeslice();
void h5_read_and_fill(const hid_t file, const hsize_t index, const char *name,
        double *out);
void read_initial_data();

#endif
