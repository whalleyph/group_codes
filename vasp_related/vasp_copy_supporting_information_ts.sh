#! /bin/bash

ts_dir=$1
rct_dir=$2
prd_dir=$3

dst_dir=$4


vasp_copy_supporting_information.sh $ts_dir $dst_dir/ts
vasp_copy_supporting_information.sh $rct_dir $dst_dir/rct
vasp_copy_supporting_information.sh $prd_dir $dst_dir/prd
