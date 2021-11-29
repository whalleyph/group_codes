#! /bin/bash
#Edited by BK

run=`ls -d RUN* |  tail -n 1`
#echo $run
if [ "$run" == "" ]; then
run=RUN1
else
id=`echo ${run:3:1}`
id=`expr $id + 1`
run="RUN$id"
#echo $id $run
fi
echo "Data is saved in $run directory."

#exit 0

mkdir $run
cp -r 0* $run
cp * $run

dirs=`ls -d ??`

for dir in $dirs; do
    cd $dir
    pwd
    cp  POSCAR POSCAR.old
    mv CONTCAR POSCAR
    cd ..
done
rm -f submit_vasp.*
