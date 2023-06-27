#!/bin/bash

N_PROC=7

a1=26
a2=0

for CAP in Nh2 Nme; do

    i=0
    for a3 in $(seq 1 26); do

        if [ $a3 == 7 ]; then continue; fi

        python cap_library_charges.py $a1 $a2 $a3 $CAP &

        if [ $((++i)) == $N_PROC ]; then i=0; wait; fi

    done
done

wait

