#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --time=48:00:00



module purge
module load GCC GCCcore python

python ~/scripts/group/group_codes/mcsqs_related/many_structs.py -c 48 -d */ -m 250