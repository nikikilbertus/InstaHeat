#!/bin/bash

function setdef {
    if [ "$1" -eq "1" ]
    then
        sed -i -e "s/$2/$3/g" main.h
    else
        sed -i -e "s/$2/MISSING_INPUT/g" main.h
    fi
}

source parameters.sh

cp main_template.h main.h

sed -i -e "s;_PATH_;${DATAPATH};g" main.h
sed -i -e "s;_IPATH_;${INITIAL_DATAPATH};g" main.h
sed -i -e "s/_IC_/${INITIAL_CONDITIONS}/g" main.h
sed -i -e "s/_PM_/${PSI_METHOD}/g" main.h
sed -i -e "s/_M_/${MASS}/g" main.h
sed -i -e "s/_A_/${A_INITIAL}/g" main.h
sed -i -e "s/_TF_/${FINAL_TIME}/g" main.h
sed -i -e "s/_GPX_/${GRIDPOINTS_X}/g" main.h
sed -i -e "s/_GPY_/${GRIDPOINTS_Y}/g" main.h
sed -i -e "s/_GPZ_/${GRIDPOINTS_Z}/g" main.h
sed -i -e "s/_LX_/${SPATIAL_LOWER_BOUND_X}/g" main.h
sed -i -e "s/_LY_/${SPATIAL_LOWER_BOUND_Y}/g" main.h
sed -i -e "s/_LZ_/${SPATIAL_LOWER_BOUND_Z}/g" main.h
sed -i -e "s/_UX_/${SPATIAL_UPPER_BOUND_X}/g" main.h
sed -i -e "s/_UY_/${SPATIAL_UPPER_BOUND_Y}/g" main.h
sed -i -e "s/_UZ_/${SPATIAL_UPPER_BOUND_Z}/g" main.h

sed -i -e "s/_MPLANCK_/${MASS_PLANCK}/g" main.h
sed -i -e "s/_MKARSTEN_/${MASS_KARSTEN}/g" main.h
sed -i -e "s/_BUF_/${WRITE_OUT_BUFFER_NUMBER}/g" main.h
sed -i -e "s/_BINS_/${POWER_SPECTRUM_BINS}/g" main.h
sed -i -e "s/_TSKIP_/${TIME_STEP_SKIPS}/g" main.h
sed -i -e "s/_XSKIP_/${STRIDE_X}/g" main.h
sed -i -e "s/_YSKIP_/${STRIDE_Y}/g" main.h
sed -i -e "s/_ZSKIP_/${STRIDE_Z}/g" main.h
sed -i -e "s/_RELTOL_/${RELATIVE_TOLERANCE}/g" main.h
sed -i -e "s/_ABSTOL_/${ABSOLUTE_TOLERANCE}/g" main.h
sed -i -e "s/_HFRAC_/${MAX_DT_HUBBLE_FRACTION}/g" main.h
sed -i -e "s/_THREADS_/${THREAD_NUMBER}/g" main.h
sed -i -e "s/_FFTW_/${FFTW_DEFAULT_FLAG}/g" main.h
sed -i -e "s/_VC_/${VERSION_CONTROL}/g" main.h

sed -i -e "s/_IM_/${INTEGRATION_METHOD}/g" main.h
sed -i -e "s/_TI_/${INITIAL_TIME}/g" main.h
sed -i -e "s/_DELT_/${DELTA_T}/g" main.h
sed -i -e "s/_MINDELT_/${MINIMAL_DELTA_T}/g" main.h
sed -i -e "s/_MAXSTEPS_/${MAX_STEPS}/g" main.h
sed -i -e "s/_MINSCAL_/${SMALLEST_SCALING}/g" main.h
sed -i -e "s/_MAXSCAL_/${LARGEST_SCALING}/g" main.h
sed -i -e "s/_BETA_/${BETA}/g" main.h
sed -i -e "s/_SAFE_/${SAFE}/g" main.h
sed -i -e "s/_SEED_/${SEED}/g" main.h
sed -i -e "s/_COUPLING_/${COUPLING}/g" main.h
sed -i -e "s/_LAMBDA_/${LAMBDA}/g" main.h

setdef ${ENABLE_FFT_FILTER} "_FILTER_" "ENABLE_FFT_FILTER"
setdef ${EVOLVE_WITHOUT_PSI} "_WOPSI_" "EVOLVE_WITHOUT_PSI"

setdef ${PHI} "_OPHI_" "OUTPUT_PHI"
setdef ${DPHI} "_ODPHI_" "OUTPUT_DPHI"
setdef ${PSI} "_OPSI_" "OUTPUT_PSI"
setdef ${DPSI} "_ODPSI_" "OUTPUT_DPSI"
setdef ${RHO} "_ORHO_" "OUTPUT_RHO"

setdef ${PHI_SUMMARY} "_OPHIS_" "OUTPUT_PHI_SMRY"
setdef ${DPHI_SUMMARY} "_ODPHIS_" "OUTPUT_DPHI_SMRY"
setdef ${PSI_SUMMARY} "_OPSIS_" "OUTPUT_PSI_SMRY"
setdef ${DPSI_SUMMARY} "_ODPSIS_" "OUTPUT_DPSI_SMRY"
setdef ${RHO_SUMMARY} "_ORHOS_" "OUTPUT_RHO_SMRY"

setdef ${PHI_POWER_SPECTRUM} "_OPHIPS_" "OUTPUT_PHI_PS"

rm main.h-e
