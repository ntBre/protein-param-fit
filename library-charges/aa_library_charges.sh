#!/bin/bash

N_PROC=12

# Loop over central three residues in Ace-Val-X-Y-Z-Val-Nme
a1=26
a4=0
a5=26

i=0
for a3 in $(seq 1 26); do
    for a2 in $(seq 1 26); do

        python aa_library_charges.py $a1 $a2 $a3 $a4 $a5 &

        if [ $((++i)) == $N_PROC ]; then i=0; wait; fi

    done
done

wait

