#!/bin/bash

N_PROC=12

# Loop over terminal two residues in X-Y-Val-Nme and Ace-Val-Y-X
a2=0
a3=26

i=0
for a1 in $(seq 1 26); do

    for charge_state in 'charged' 'neutral'; do

        # N terminus
        python terminal_library_charges.py $a1 $a2 $a3 'N' $charge_state &

        # C terminus
        python terminal_library_charges.py $a3 $a2 $a1 'C' $charge_state &

        if [ $((++i)) == $N_PROC ]; then i=0; wait; fi

    done
done

wait

