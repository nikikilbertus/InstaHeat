#!/bin/bash
# This script uses the template.sbatch file to create multiple sbatch files for
# different simulation parameters. This example varies two parameters.

for var1 in 32 48 64 96
do
    for var2 in 0.005 0.002 0.001 0.0005 0.0002 0.0001 0.00005 0.00002 0.00001
    do
        NAME=scr${var1}.sbatch
        cp template.sbatch $NAME
        sed -i "s/XXX/${var1}/g" $NAME
        sed -i "s/YYY/${var2}/g" $NAME
        if [ $# -gt 0 ]
        then
            sed -i "s/HOURS/$1/g" $NAME
        else
            sed -i "s/HOURS/48/g" $NAME
        fi
        if [ $# -gt 1 ]
        then
            sed -i "s/RUNFILE/$2${var1}_${var2}/g" $NAME
        else
            sed -i "s/RUNFILE/run${var1}_${var2}/g" $NAME
        fi
        sbatch "$NAME"
    done
done
