#! /bin/bash

src_dir=$1
dst_dir=$2

files='INCAR KPOINTS POSCAR POTCAR CONTCAR OSZICAR vasp.out'

for f in $files; do
    mkdir -p $dst_dir
    cp -i $src_dir/$f $dst_dir/$f
done
