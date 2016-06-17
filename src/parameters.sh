# common
DATAPATH="../../data/test_mass.h5"
INITIAL_DATAPATH="../../data/karsten/data_64psi_3.dat"
INITIAL_CONDITIONS="IC_FROM_BUNCH_DAVIES"
MASS="4.0" # bunch davies
# MASS="9.990962364968714" # bunch davies s.t. kmax,kmin in resonance at some a
# MASS="0.002000003836216" # compare_2 compare_psi
# MASS="1.0e-2" # data_64psi_{1,2,3}
# MASS="1.0e-3" # data_64psi_4
# MASS="1.0e-4" # data_64psi_5
A_INITIAL="1.0" # bunch davies
# A_INITIAL="6.1625e2" # compare_2, compare_psi, 5450
# A_INITIAL="1.05249e3" # compare_2, compare_psi, 5500
# A_INITIAL="4.09376e3" # compare_2, compare_psi, 6000
# A_INITIAL="1.370074629050061e5" # data_64_0
# A_INITIAL="6.227611966276128e4" # data_64psi_{1,4,5}
# A_INITIAL="1.868327477603207e4" # data_64psi_2
# A_INITIAL="6.227611966276128e3" # data_64psi_3
FINAL_TIME="1.0e4"
GRIDPOINTS_X="96"
GRIDPOINTS_Y="96"
GRIDPOINTS_Z="96"
SPATIAL_LOWER_BOUND_X="0.0"
SPATIAL_UPPER_BOUND_X="10.0"
SPATIAL_LOWER_BOUND_Y="0.0"
SPATIAL_UPPER_BOUND_Y="10.0"
SPATIAL_LOWER_BOUND_Z="0.0"
SPATIAL_UPPER_BOUND_Z="10.0"

# uncommon
FFTW_DEFAULT_FLAG="FFTW_PATIENT"
VERSION_CONTROL="VERSION_CONTROL_HG"
ENABLE_FFT_FILTER="0"
ENABLE_GW="1" # enable gravitational wave extraction
WRITE_OUT_BUFFER_NUMBER="5"
POWER_SPECTRUM_BINS="50"
TIME_STEP_SKIPS="1"
STRIDE_X="1"
STRIDE_Y="1"
STRIDE_Z="1"
INFLATON_MASS="5.0e-3"
MASS_KARSTEN="1.0e-2"
RELATIVE_TOLERANCE="1.0e-10"
ABSOLUTE_TOLERANCE="1.0e-14"
MAX_DT_HUBBLE_FRACTION="1.0e-2"
THREAD_NUMBER="0"
BUNCH_DAVIES_CUTOFF="0"

# rare
INTEGRATION_METHOD="DOPRI853"
INITIAL_TIME="0.0"
DELTA_T="1.0e-7"
MINIMAL_DELTA_T="1.0e-9"
MAX_STEPS="1e15"
SMALLEST_SCALING="0.333"
LARGEST_SCALING="3.0"
BETA="0.0"
SAFE="0.9"
SEED="93"
COUPLING="1.0"
LAMBDA="1.876e-4"
MAX_RUNTIME="30" # slightly less than 48h

# outputs (0=no output, 1=output)
# full fields
OUTPUT_PHI="0"
OUTPUT_DPHI="0"
OUTPUT_PSI="0"
OUTPUT_DPSI="0"
OUTPUT_RHO="0"

# power spectra
OUTPUT_PHI_PS="1"
OUTPUT_PSI_PS="1"
OUTPUT_RHO_PS="1"

# summaries
OUTPUT_PHI_SMRY="1"
OUTPUT_DPHI_SMRY="1"
OUTPUT_PSI_SMRY="1"
OUTPUT_DPSI_SMRY="1"
OUTPUT_RHO_SMRY="1"
OUTPUT_PRESSURE_SMRY="1"
OUTPUT_H1_SMRY="0" # not fully functional yet
OUTPUT_H2_SMRY="0" # not fully functional yet

OUTPUT_CONSTRAINTS="0"
