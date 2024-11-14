#!/bin/bash

LEAP_DIR=$AMBERHOME/dat/leap
LEAP_LIB_DIR=$LEAP_DIR/lib

QAMBER_FILE=tmp.dat

rm -f $QAMBER_FILE

for LIB_FILE in amino12.lib aminont12.lib aminoct12.lib; do

    # Print residue name, atom name, and partial charge
    awk '
        /unit.atomspertinfo/ {do_print=0}
        do_print {printf "%4s %6s %9.6f\n", resname, $1, $8}
        /unit.atoms / {
            split($1, s, ".")
            resname=s[2]
            if (resname != "CYM" && resname != "HYP" && resname != "CHYP") do_print=1
        }
    ' $LEAP_LIB_DIR/$LIB_FILE >> $QAMBER_FILE

done

# Put cap charges at beginning of file, rename NHE to NH2
TMP_FILE=${QAMBER_FILE}.tmp
mv $QAMBER_FILE $TMP_FILE
awk '$1 == "ACE"' $TMP_FILE > $QAMBER_FILE
awk '$1 == "NHE"' $TMP_FILE | sed 's/NHE/NH2/' >> $QAMBER_FILE
awk '$1 == "NME"' $TMP_FILE >> $QAMBER_FILE
awk '$1 != "ACE" && $1 != "NHE" && $1 != "NME"' $TMP_FILE >> $QAMBER_FILE
rm $TMP_FILE
