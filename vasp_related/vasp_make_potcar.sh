#! /bin/bash

eles=`head -n 1 POSCAR`

echo 'nullifying POTCAR'
> POTCAR
for ele in $eles; do
  echo $ele
  echo 'appending to POTCAR'
  cat $HOME/tools/vasp/pseudopotentials/potpaw_PBE/$ele/POTCAR >> POTCAR
done
