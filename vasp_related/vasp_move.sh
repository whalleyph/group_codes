#! /bin/bash

#if [ $# != seq ]; then
#  echo "Usage: vasp_move.sh < sequence number >"
#  exit seq
#fi

ls *.tgz 2> /dev/null | grep ^[0-9].*

if [ $? -eq 0 ]; then
  seq=$(ls *.tgz | grep ^[0-9].* | tail -1 | cut -c -1)
  let seq=seq+1
  echo $seq
else
  seq=0
  echo $seq
fi


cp -iv CONTCAR CONTCAR.$seq
mv -iv POSCAR POSCAR.$seq
mv -iv CONTCAR POSCAR
#cp -iv CHGCAR CHGCAR.$seq
cp -iv INCAR INCAR.$seq
mv -iv OSZICAR OSZICAR.$seq
mv -iv OUTCAR  OUTCAR.$seq
mv -iv vasp.out vasp.out.$seq
cp MODECAR MODECAR.$seq
cp NEWMODECAR NEWMODECAR.$seq
mv NEWMODECAR MODECAR

if [ -f $seq.tgz ]; then
  echo $seq'.tgz exists'
  exit 1
else
  tar czvf $seq.tgz *.$seq
fi
rm *.$seq
