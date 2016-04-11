# common
DATAPATH="../../../data/testbunch.h5"
INITIAL_DATAPATH="../../../data/karsten/data_64psi_2.dat"
INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"
PSI_METHOD="PSI_HYPERBOLIC"
MASS="1.0"
# MASS="0.002000003836216" # compare_2 compare_psi
A_INITIAL="1.0"
# A_INITIAL="6.1625e2" # compare_2, compare_psi, 5450
# A_INITIAL="1.05249e3" # compare_2, compare_psi, 5500
# A_INITIAL="4.09376e3" # compare_2, compare_psi, 6000
# A_INITIAL="1.370074629050061e5" # data_64_0
# A_INITIAL="6.227758258677358e4" # data_64psi_1
# A_INITIAL="1.868327477603207e4" # data_64psi_2
# A_INITIAL="6.227611966276128e3" # data_64psi_3
FINAL_TIME="1.0e-5"
GRIDPOINTS_X="128"
GRIDPOINTS_Y="128"
GRIDPOINTS_Z="128"
SPATIAL_LOWER_BOUND_X="0.0"
SPATIAL_UPPER_BOUND_X="10.0"
SPATIAL_LOWER_BOUND_Y="0.0"
SPATIAL_UPPER_BOUND_Y="10.0"
SPATIAL_LOWER_BOUND_Z="0.0"
SPATIAL_UPPER_BOUND_Z="10.0"

# uncommon
FFTW_DEFAULT_FLAG="FFTW_ESTIMATE"
VERSION_CONTROL="VERSION_CONTROL_HG"
ENABLE_FFT_FILTER="0"
WRITE_OUT_BUFFER_NUMBER="1"
POWER_SPECTRUM_BINS="50"
TIME_STEP_SKIPS="1"
STRIDE_X="1"
STRIDE_Y="1"
STRIDE_Z="1"
MASS_PLANCK="2.0e6"
MASS_KARSTEN="1.0e-2"
RELATIVE_TOLERANCE="1.0e-8"
ABSOLUTE_TOLERANCE="1.0e-12"
MAX_DT_HUBBLE_FRACTION="1.0e-2"
THREAD_NUMBER="0"

# rare
INTEGRATION_METHOD="DOPRI853"
EVOLVE_WITHOUT_PSI="0"
INITIAL_TIME="0.0"
DELTA_T="1.0e-5"
MINIMAL_DELTA_T="1.0e-9"
MAX_STEPS="1e15"
SMALLEST_SCALING="0.333"
LARGEST_SCALING="6.0"
BETA="0.0"
SAFE="0.9"
SEED="93"
COUPLING="1.0"
LAMBDA="1.876e-4"

# outputs (0=no output, 1=output)
PHI="1"
DPHI="1"
PSI="1"
DPSI="1"
RHO="1"

POWER_SPECTRUM="1"

PHI_MEAN="1"
DPHI_MEAN="1"
PSI_MEAN="1"
DPSI_MEAN="1"
RHO_MEAN="1"

PHI_VARIANCE="1"
DPHI_VARIANCE="1"
PSI_VARIANCE="1"
DPSI_VARIANCE="1"
RHO_VARIANCE="1"
