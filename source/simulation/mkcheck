#!/bin/bash
cp main.h main.h.keep
for num in 32 64 96 128 192 256
do
    cp main.h.keep main.h.tmp
    sed -i -e "s/=GRIDNUM=/${num}/g" main.h.tmp
    for tol in {3..13}
    do
        cp main.h.tmp main.h
        sed -i -e "s/=TOLNUM=/${tol}/g" main.h
        make
        ./first_steps
    done
done
rm main.h.tmp
rm main.h-e
rm main.h.tmp-e
mv main.h.keep main.h
make clean
