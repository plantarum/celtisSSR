#!/bin/bash

KMIN=1
KMAX=10
REPS=20

for REP in $(seq 1 $REPS)
do
    for i in $(seq $KMIN $KMAX)
    do
        structure -K $i -o output/output${i}_${REP} > output/screen${i}_${REP}.txt
    done
done

