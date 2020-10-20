#!/bin/bash

model="example-datasets/recon22_with_biomass.sfba"
actmat="example-datasets/matrix_5_samples.tab" #$1
soldir="example-datasets/solutions" #"$solutions_dir/${arrINFILE[0]}"
actmatext=(${actmat//.tab/_ext.tab})

### deprecated parameters ###
#solutions_dir="example-datasets/solutions_dir_biomass"
#arrIN=(${actmat///// })
#actmatfile=${arrIN[1]}
#arrINFILE=(${actmatfile//.// })
#############################

mkdir $soldir

echo "Dir created."
python2.7 extend_header.py $actmat
echo "Headers extended."
python2.7 define.py $model $actmatext
mv *.lp $soldir
echo "Problems defined."
python2.7 optimize.py
echo "Problem optimized."
#mv *.sol $soldir
#rm *.tmp
python2.7 get_solutions.py $model $actmatext $soldir
mv matrix_* $soldir
