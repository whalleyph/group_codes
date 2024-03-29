#!/bin/sh
#lipai@mail.ustc.edu.cn
#generate xsf files using OUTCAR and POSCAR
# Modified by Tanmoy on March 9, 2022

if [ $# = 0 ]; then
    out="OUTCAR"
else
    out=$1
fi
echo $out
rm *.temp

# total number of ions in the system
num_atom=`grep -m 1 "NIONS =" $out|awk '{print $12}'`
echo "Total number of ions: $num_atom"

# create temp files for writing xyz files
typenum=`grep -m 1 'ions per type' $out |head -1 |awk '{print NF}' `
typenum=$(($typenum-4)) # how many types of ions
for i in `seq $typenum`
do
    elename=`grep -m $i POTCAR $out |tail -1 |awk '{print $3}'`
    j=$(($i+4))
    elenum=`grep -m 1 "ions per type" $out |awk -v j=$j '{print $j}'`
    echo $elename $elenum
    for j in `seq $elenum`
    do
        echo $elename  >> type.temp
    done
done
grep -A 3 -m 1 "direct lattice vectors" $out \
|tail -3 |awk '{printf("%f %f %f \n",$1,$2,$3)}' >primvec.temp
#grep "energy  without " $out |awk '{print $4}' >energy.temp
grep E0 vasp.out | awk '{print $5}' > energy.temp  # Tanmoy
awk '/POSITION/,/drift/{ if(NF==6) print $0 }' $out  > pos.temp

lines=`wc pos.temp|awk '{print $1}'`
num_str=`echo "$lines/$num_atom" |bc` # how many structures
echo "Number of structures: $num_str"
if [  -f all.xyz ]; then
    rm all.xyz
fi

for i in `seq $typenum`
do
    elename=`grep -m $i POTCAR $out |tail -1 |awk '{print $3}'`
    j=$(($i+4))
    elenum=`grep -m 1 "ions per type" $out |awk -v j=$j '{print $j}'`
    echo $elename $elenum
    for j in `seq $elenum`
    do
        echo $elename  >> type.temp
    done
done

echo "num of str: $num_str"

for i in `seq $num_str`
do
    energy=`head -n $i energy.temp|tail -1`
    echo "# total energy = $energy eV" >> str_$i.xsf
    echo " " >> str_$i.xsf
    echo "CRYSTAL" >> str_$i.xsf
    echo "PRIMVEC" >> str_$i.xsf
    cat primvec.temp >> str_$i.xsf
    echo "PRIMCOORD" >> str_$i.xsf
    echo "$num_atom 1" >> str_$i.xsf
    end=`echo "$i*$num_atom" |bc `
    head -n $end pos.temp|tail -n $num_atom >pos_i.temp
    paste type.temp pos_i.temp >> str_$i.xsf
    mv str_$i.xsf $out-$i.xsf
done

rm *.temp
if [ ! -d "struc" ]; then
    mkdir struc
fi
mv *xsf struc

# Tanmoy added these 5 lines below - needed for AENET
cd struc
ls -t | tail -n +2 | xargs rm --
atom=`grep PRIMCOORD -A 1 *.xsf | tail -1 | awk '{print $1}'`
sum=$(($atom+10))
sed -i "$sum,\$d" *.xsf
