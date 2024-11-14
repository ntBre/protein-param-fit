#!/bin/bash

LEAP_DIR=$AMBERHOME/dat/leap
LEAP_LIB_DIR=$LEAP_DIR/lib
LEAP_PARM_DIR=$LEAP_DIR/parm

LJAMBER_FILE=tmp-lj.dat
NBAMBER_FILE=tmp-nb.dat

# Print Amber atom type, Lennard-Jones radius, and Lennard-Jones epsilon
awk '
    do_print {
        if ($1 == "N" || $1 == "C*") {
            split(redundant_atom_types[$1], lj_types, " ")
            for (a in lj_types) printf "%2s %8.6f %8.6f\n", lj_types[a], $2, $3
        }
        else printf "%2s %8.6f %8.6f\n", $1, $2, $3
    }
    $1 == "I" {do_print=0}
    ! do_print && ($1 == "N" || $1 == "C*") {redundant_atom_types[$1] = $0}
    $1 == "MOD4" {do_print=1}
' $LEAP_PARM_DIR/parm10.dat > $LJAMBER_FILE
awk '
    do_print {printf "%2s %8.6f %8.6f\n", $1, $2, $3}
    $1 == "CO" {do_print=0}
    $1 == "NONB" {do_print=1}
' $LEAP_PARM_DIR/frcmod.ff14SB >> $LJAMBER_FILE

rm -f $NBAMBER_FILE

for LIB_FILE in amino12.lib aminont12.lib aminoct12.lib; do

    # Print residue name, atom name, partial charge. LJ radius, LJ epsilon
    awk '
        FNR==NR {lj_radius[$1] = $2; lj_epsilon[$1] = $3}
        /unit.atomspertinfo/ {do_print=0}
        do_print {
            atom_type = substr($2, 2, length($2)-2)
            atom_name = substr($1, 2, length($1)-2)
            printf "%4s %4s %9.6f %8.6f %8.6f\n", \
                resname, atom_name, $8, lj_radius[atom_type], lj_epsilon[atom_type]
        }
        /unit.atoms / {
            split($1, s, ".")
            resname=s[2]
            if (resname != "CYM" && resname != "HYP" && resname != "CHYP") do_print=1
        }
    ' $LJAMBER_FILE $LEAP_LIB_DIR/$LIB_FILE >> $NBAMBER_FILE

done

# Put cap parameters at beginning of file, rename NHE to NH2
TMP_FILE=${NBAMBER_FILE}.tmp
mv $NBAMBER_FILE $TMP_FILE
awk '$1 == "ACE"' $TMP_FILE > $NBAMBER_FILE
awk '$1 == "NHE"' $TMP_FILE | sed 's/NHE/NH2/' >> $NBAMBER_FILE
awk '$1 == "NME"' $TMP_FILE >> $NBAMBER_FILE
awk '$1 != "ACE" && $1 != "NHE" && $1 != "NME"' $TMP_FILE >> $NBAMBER_FILE
rm $TMP_FILE
