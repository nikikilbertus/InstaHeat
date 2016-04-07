#!/bin/bash

source parameters.sh

cp main_template.h main.h

sed -i -e "s;_PATH_;${DATAPATH};g" main.h
sed -i -e "s;_IPATH_;${INITIAL_DATAPATH};g" main.h
sed -i -e "s/_IC_/${INITIAL_CONDITIONS}/g" main.h
sed -i -e "s/_METHOD_/${PSI_METHOD}/g" main.h
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

if [ "${ENABLE_FFT_FILTER}" -eq "1" ]
then
    sed -i -e "s/_FILTER_/ENABLE_FFT_FILTER/g" main.h
else
    sed -i -e "s/_FILTER_/DISABLE_FFT_FILTER/g" main.h
fi

rm main.h-e
#TODO: output stuff
