#!/bin/bash

# Author - Tanmoy - May 16, 2022

if [ ! -d "collect_CoEn" ]; then
      mkdir collect_CoEn
fi

for i in Step*/
do 
cd $i/
bash ~/SCRIPTS/collect_ConfigEnum.sh
mv CoEn* ../collect_CoEn/
cd ../
done
