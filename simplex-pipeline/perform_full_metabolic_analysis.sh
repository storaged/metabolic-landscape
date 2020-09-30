#!/bin/bash

model="recon22_with_biomass.sfba"
actmat=$1
solutions_dir="solutions_dir_biomass"
arrIN=(${actmat///// })
actmatfile=${arrIN[1]}
arrINFILE=(${actmatfile//.// })
soldir="$solutions_dir/${arrINFILE[0]}"

actmatext=(${actmat//.tab/_ext.tab})

mkdir $soldir
echo "Dir created."
python2.7 extend_header.py $actmat
echo "Headers extended."
python2.7 define.py $model $actmatext
echo "Problems defined."
python2.7 optimize.py
echo "Problem optimized."
mv *.lp *.sol $soldir
rm *.tmp
python2.7 get_solutions.py $model $actmatext $soldir
mv matrix_* $soldir
