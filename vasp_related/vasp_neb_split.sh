#! /bin/bash -x

if [ $# -ne 1 ]; then
  echo "Usage: split.sh <No. of images>"
  echo "Make sure that POSCARs in 0* folders have the same format (vasp4/vasp5)."
  exit 1
fi

mkdir splitneb

dirs=$(seq 0 "$1")
echo $dirs
let np1=$1+1
let np2=2*$1+2
cp 00/POSCAR 00/CONTCAR
cp 0$np1/POSCAR 0$np1/CONTCAR
mkdir tmpsplit
for i in $dirs; do
  let ip1=$i+1
  cp 0$i/CONTCAR tmpsplit/POS1
  cp 0$ip1/CONTCAR tmpsplit/POS2
  cd tmpsplit
  nebmake.pl POS1 POS2 1
  #nebmakecorrect.sh 2
  mv 01 ../splitneb/0$i.y
  rm -rf ./*
  cd ..
  pwd
done
rm -rf tmpsplit
#rm 00/CONTCAR 0$np1/CONTCAR
for i in $dirs; do
  cp -r 0$i splitneb/0$i.x
  mv splitneb/0$i.x/CONTCAR splitneb/0$i.x/POSCAR
done

cp -r 0$np1 splitneb/0$np2; mv splitneb/0$np2/CONTCAR splitneb/0$np2/POSCAR
cd splitneb
cnt=0
for i in $dirs; do
  mv 0$i.x 0$cnt
  let cnt=$cnt+1
  mv 0$i.y 0$cnt
  let cnt=$cnt+1
done
