#! /bin/bash

# list of required files
req_files="POSCAR.rct.vasp POSCAR.prd.vasp INCAR.min INCAR.neb KPOINTS sub_min sub_neb POTCAR"

# check if all the required files are present
for f in $req_files; do
    if [ ! -f $f ]; then
       echo $f does not exist 
       exit 1
    fi
done

# make the directories for reactant and product
mkdir rct prd

# copy the POSCAR files
cp POSCAR.rct.vasp rct/POSCAR
cp POSCAR.prd.vasp prd/POSCAR

cp POTCAR rct/
cp POTCAR prd/
# exit if the POTCARS differ
diff -q rct/POTCAR prd/POTCAR
if [ $? -ne  0 ]; then
    exit 1
fi

# copy the INCAR files
cp INCAR.min rct/INCAR
cp INCAR.min prd/INCAR

# copy the KPOINTS file
cp KPOINTS rct/
cp KPOINTS prd/

# copy the submit file
cp sub_min rct/sub
cp sub_min prd/sub

cd rct
rct_id=`qsub sub`

cd ../prd
prd_id=`qsub sub`

cd ..
qsub -W depend=afterok:$rct_id:$prd_id sub_neb
