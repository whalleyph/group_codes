#!/bin/bash
# Author - Tanmoy
# This script assumes that all the associated scripts (along with itself) are placed in somewhere whose PATH is given in your ~/.bashrc file
# Run this script in your AIRSS calcuilation directory; NOTE: you also need "check_conv.sh" and "generate_xsf_aenet.sh" scripts
# The script first collects the converged structures into a folder called 'converged' and create a list. In the second step, the script goes through all the converged structures and convert them to .xsf file format and place them in collect_xsf folder

sh check_conv.sh

cd converged
printf '%s\n' * > list_converged.txt
cd ../

sed -i "s/.res//g" "list_converged.txt"

mkdir collect_xsf

while read -r line
do
    cd $line/
    sh ../generate_xsf_aenet.sh
    cd struc/
    mv *.xsf $line.xsf
    cp -rp *.xsf ../../collect_xsf/
    cd ..
    rm -rf struc/
    cd ../
done < list_converged.txt

rm list_converged.txt
