#ifndef __FILEHANDLING__
#define __FILEHANDLING__

#include "hdf5.h"

#define     H5_PHI_NAME                 ("phi")
#define     H5_DPHI_NAME                ("dphi")
#define     H5_PSI_NAME                 ("psi")
#define     H5_DPSI_NAME                ("dpsi")
#define     H5_RHO_NAME                 ("rho")
#define     H5_POWER_SPECTRUM_NAME      ("power_spectrum")
#define     H5_TIME_NAME                ("time")
#define     H5_A_NAME                   ("a")
#define     H5_PHI_MEAN_NAME            ("phi_mean")
#define     H5_PHI_VARIANCE_NAME        ("phi_variance")
#define     H5_DPHI_MEAN_NAME           ("dphi_mean")
#define     H5_DPHI_VARIANCE_NAME       ("dphi_variance")
#define     H5_PSI_MEAN_NAME            ("psi_mean")
#define     H5_PSI_VARIANCE_NAME        ("psi_variance")
#define     H5_DPSI_MEAN_NAME           ("dpsi_mean")
#define     H5_DPSI_VARIANCE_NAME       ("dpsi_variance")
#define     H5_RHO_MEAN_NAME            ("rho_mean")
#define     H5_RHO_VARIANCE_NAME        ("rho_variance")
#define     H5_MASS_NAME                ("mass")
#define     H5_DIMENSION_NAME           ("dimension")
#define     H5_SEED_NAME                ("seed")
#define     H5_STRIDES_TIME_NAME        ("strides_time")
#define     H5_STRIDES_SPACE_NAME       ("strides_space")
#define     H5_TOLERANCES_NAME          ("tolerances")
#define     H5_GRIDPOINTS_INTERNAL_NAME ("gridpoints_internal")
#define     H5_GRIDPOINTS_OUTPUT_NAME   ("gridpoints_output")
#define     H5_COMMIT_HASH_NAME         ("commit_hash")
#define     H5_RUNTIME_TOTAL_NAME       ("runtime_total")
#define     H5_RUNTIME_FFTW_NAME        ("runtime_fftw")
#define     H5_RUNTIME_FFTWPLAN_NAME    ("runtime_fftwplan")
#define     H5_RUNTIME_FILTER_NAME      ("runtime_filter")
#define     H5_RUNTIME_ELLIPTIC_NAME    ("runtime_elliptic")
#define     H5_RUNTIME_WRITEOUT_NAME    ("runtime_writeout")
#define     H5_RUNTIME_STEPPER_NAME     ("runtime_stepper")
#define     H5_STEPS_TOTAL_NAME         ("steps_total")
#define     H5_STEPS_OK_NAME            ("steps_ok")
#define     H5_STEPS_BAD_NAME           ("steps_bad")

void h5_create_empty_by_path(const char *name);
void h5_create_dset(const hsize_t rank, const hsize_t *dim,
        const hsize_t *max_dim, const hsize_t *chunk, hsize_t *dset_id,
        const char *name);
void h5_write_parameter(const char *name, const double *val, size_t N);
void h5_write_all_buffers(const hsize_t Nt);
void h5_get_extent(hsize_t *max, hsize_t *cur);
void h5_close();
void save();
void h5_read_timeslice(double t, double *f);
void h5_read_and_fill(const hid_t file, const hsize_t index, const char *name,
        double *out);

#endif
