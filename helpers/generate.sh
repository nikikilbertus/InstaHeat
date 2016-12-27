#!/bin/bash

# This script compiles InstaHeat wit varying parameters. Recall that InstaHeat
# sets most parameter values during compile time. Therefore we cannot run the
# same binary with different parameters, but have to compile separate binaries
# each time.

module load hdf5
module load gsl/2.1

# declare -a arr=("RKF45" "RKCK" "DOPRI89" "DOPRI853")

for var1 in 32 48 64 96
    do
    # for var2 in "${arr[@]}"
    for var2 in 0.005 0.002 0.001 0.0005 0.0002 0.0001 0.00005 0.00002 0.00001
    do
        cp parameters.sh parameters_tmp.sh
        sed -i "s/XXX/${var1}/g" parameters_tmp.sh
        sed -i "s/YYY/${var2}/g" parameters_tmp.sh
        make parameters=parameters_tmp.sh
        mv run /scratch/users/username/run${var1}_${var2}
        make clean
    done
done
rm parameters_tmp.sh
