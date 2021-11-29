#!/bin/bash

echo -n "This will clean the current VASP calculation. Are you sure? (y: yes): "
read ans

if [ "$ans" != "y" ]; then
echo "Nothing is done..."
exit 1
fi

/bin/rm -f CHG \
CHGCAR \
CONTCAR \
DOSCAR PROCAR dos.png bandstr*.png \
DYNMAT \
EIGENVAL \
OSZICAR \
OUTCAR \
PCDAT \
WAVECAR \
XDATCAR \
IBZKPT \
vasp.out \
vasprun.xml \
OUTCAR.xyz STOPCAR temp.tmp mag.tmp cpu.tmp rms*.dat \
opted.xyz REPORT HILLSPOT \
NEWMODECAR CENTCAR DIMCAR dimer.dat dimer.eps exts.dat mep* movie* neb* spline* 0*/WAVECAR \
exts.dat  mep.eps movie* neb* spline.dat eigs.dat freq.* modes.dat \
ACF.dat  AECCAR0  AECCAR1  AECCAR2  AVF.dat  BCF.dat  CHGCAR_sum \
nodefile.out  vasp.nodes deneme.e* deneme.o*  scr_dir home_dir vasp_std* \
slurm-*.out  machine.file.* all-traj.xyz *.dat
