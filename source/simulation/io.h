#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "main.h"

/**
 * @file io.h
 * @brief Contains preprocessor definitions of the names of the datasets in the
 * output `.h5` file and function declarations for `io.c`.
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
#define     H5_INFLATON_MASS_NAME           ("inflaton_mass")
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
#define     H5_RUNTIME_COPY_BUFFER_NAME     ("runtime_copy_buffer")
#define     H5_RUNTIME_STEPPER_NAME         ("runtime_stepper")
#define     H5_RUNTIME_CSTR_NAME            ("runtime_constraints")
#define     H5_RUNTIME_SMRY_NAME            ("runtime_summaries")
#define     H5_RUNTIME_STIFF_NAME           ("runtime_stiffness_check")
#define     H5_COUNTER_RHS                  ("calls_rhs")
#define     H5_STEPS_TOTAL_NAME             ("steps_total")
#define     H5_STEPS_OK_NAME                ("steps_ok")
#define     H5_STEPS_BAD_NAME               ("steps_bad")
#define     H5_MAX_DT_HUBBLE_FRACTION_NAME  ("max_dt_hubble_fraction")

void h5_create_empty_by_path();
void h5_write_parameter(const char *name, const double *val, const size_t N);
void save();
void h5_close();
void h5_read_timeslice();
void read_initial_data();

#endif
