#!/bin/sh

if [ ! -d "converged" ]; then
      mkdir converged
fi

if [ ! -d "unconverged" ]; then
      mkdir unconverged
fi

for i in input-*/
do
        cd $i
        VAR=`grep "reached required accuracy" vasp.out | tail -1 | awk '{if(NF>0){print "T" }else{print "F"}}'`
        #echo $VAR 
        if [ "$VAR" = "T" ] ; then
                #E0=`grep E0 OSZICAR | tail -1 | awk '{print $5}'`
                #echo $i.vasp $E0
                #echo done!
                cp -rp *.res ../converged/
		cd ../

        else
        	#mv CONTCAT POSCAR
                cp -rp *.res ../unconverged/
                #qsub run.sh
                #echo $i
                cd ../
        fi
done
