#!/bin/bash

# Authors - Bora & Tanmoy - May16, 2022
# This script was written with the intention to first identify and then collect corresponding VASP calculation folder to a given .res file from a "CONFIG-ENUM" calculation directory (where StepXX folders exist).
# To use this code run: run_collect_ConfigEnum.sh in the directory where StepX folders exist. Or one can separately run this code inside a given StepX folder.

for i in *.res
do
VAR1=`head -1 $i | awk '{print $5}' | awk -F. '{print $1"."substr($2,1,4)}'`
VAR2=$(grep -- "$VAR1" `find -name degeneracies` | head -n1  | awk '{print $1}' | awk -F. '{print $1}')
VAR3=$(basename $i .res)
mkdir CoEn-$VAR3
cp -rp input_0/$VAR2/* CoEn-$VAR3/
cp -rp $i CoEn-$VAR3/
done
