#!/bin/bash
LINES=$(tput lines)
FG=`tput setaf 30` # 0=Black, 1=Red, 3=Green, 3=Yellow, 4=Blue, 5=Magenta, 6=Cyan, 7=White
RESET=`tput sgr0`
### Most of the comments are some experimental progress indicator shenanigans
### Please ignore, or uncomment if you are feeling brave.
### It's just all writing to the terminal and so won't affect the base
### functionality of the script.


# set_window ()
# {
#     # Create a virtual window that is two lines small at the bottom.
#     tput csr 0 $(($LINES-2))
# }

# print_status ()
# {
#     # Move cursor to last line in your screen
#     tput cup $LINES 0;

#     echo -n "--- $current of $conv_count ---"

#     # Move cursor to home position, back in virtual window
#     tput cup 0 0
# }


function make_xsf () {
    struct_name=$1
    out="OUTCAR"
    rm -f *.temp
    
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
        mv str_$i.xsf $struct_name-$i.xsf
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
    cd ../
    
}

conv=()
for i in */
do
    cd $i
    VAR=`grep "reached required accuracy" vasp.out | tail -1 | awk '{if(NF>0){print "T" }else{print "F"}}'`
    if [ "$VAR" = "T" ] ; then
        conv[${#conv[@]}]="$i"
    fi
    cd ../
done

conv_count=${#conv[@]}
echo "Found $conv_count converged structures."
sleep 2

#clear -x # clear -x maintain scrollback

#set_window

for j in "${!conv[@]}"
do
    k=${conv[$j]}
    cd $k
    let "j=j+1"
    current=$j
    echo
    echo "${FG}$j of $conv_count ${RESET}"
    echo "${k%/}"
    
    if [ -f "struc/*.xsf" ]; then
        echo "XSF found. Skipping..."
    else
        make_xsf ${k%/}
    fi
    cd ..
    #print_status
done

mkdir converged-xsf
cp -np */struc/*.xsf converged-xsf
