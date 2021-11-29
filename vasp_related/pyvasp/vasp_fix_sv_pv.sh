#! /bin/bash

sv='Ba Ca Cs K Rb Sr Y Zr'
pv='Nb'

for ele in $sv; do
    oldfile=POSCAR."$ele".vasp
    newfile=POSCAR."$ele"_sv.vasp
    > $newfile
    echo "$ele"_sv >> $newfile
    tail -n +2 POSCAR.$ele.vasp >> $newfile
    rm $oldfile
done

for ele in $pv; do
    oldfile=POSCAR."$ele".vasp
    newfile=POSCAR."$ele"_sv.vasp
    > $newfile
    echo "$ele"_pv >> $newfile
    tail -n +2 $oldfile >> $newfile
    rm $oldfile
done
