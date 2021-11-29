#!/usr/bin/env python3
from subprocess import Popen,PIPE # as popen # check_output
from sys import exit,stdout,argv,version_info
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
import numpy as np
from copy import deepcopy as dc
from ase.neighborlist import neighbor_list
import ase.io
from ase import Atoms,Atom
import argparse,math
from ase.visualize import view
from os import system


def Popen4(cmd):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate()
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err

parser = argparse.ArgumentParser(description='Script for creating a homogenous gas using the ASE interface.')

parser.add_argument('-n','--noAtoms', type=int,required=1, help='No of atoms in the system')

parser.add_argument('-at','--atype', type=str,required=True, help='Atom type to be added to the box')

parser.add_argument('-d','--dims', nargs='*', required=1, type=float,help="Box dimensions, e.g. -d 15 or -d 15 20 25")

parser.add_argument('-ai','--airss', default=False,action='store_true', help='Use AIRSS to generate a random cluster')

parser.add_argument('-aopt','--aiOpt', default=False,action='store_true', help='Optimise the geometry created by AIRSS using LAMMPS')

args = parser.parse_args()

noAts=args.noAtoms
atype=args.atype

if len(args.dims)==1:
    dims=[float(args.dims[0]) for i in range(3)]
elif len(args.dims)==3:
    dims=[float(args.dims[i]) for i in range(3)]
else: exit()

if args.airss:
    str1="""%%BLOCK LATTICE_CART
%.1f 0 0
0 %.1f 0
0 0 %.1f
#FIX
%%ENDBLOCK LATTICE_CART

%%BLOCK POSITIONS_FRAC
%s 0.0 0.0 0.0 # %s1 %% NUM=%d
%%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#MINSEP=1.4
#CLUSTER
#POSAMP=3.0
"""%(dims[0],dims[1],dims[2],atype,atype,noAts)
    open('airss.cell','w').writelines(str1)
    build='-build'
    par=''
    if args.aiOpt:
        print('Geometry will be opted using LAMMPS...')
        str2="""
#pair_style reax/c NULL #lmp_control
#pair_coeff * * ffield-JPCA2015-Ganesh C
#fix 1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

pair_style edip/multi
#pair_style edip
pair_coeff * * SiC.edip C
"""
        open('airss.pp','w').writelines(str2)
        build=''
        #par='-mpinp 16'

    cmd='airss.pl  -lammps -cluster -max 1 -seed airss %s %s'%(build,par)
    system(cmd) #no geom opt will be done if -build
    exit()


#Initate the Atoms object
atoms=Atoms('H',cell=dims,pbc=1);del atoms[0]

#xd=at*dims[0]/float(noAts)/3.
#yd=at*dims[1]/float(noAts)/3.
#zd=at*dims[2]/float(noAts)/3.

steps=int(np.ceil(math.pow(float(noAts),1./3.)))
#xd=dims[0]/steps
#yd=dims[1]/steps
#zd=dims[2]/steps
disp=[(dims[i]-0.75)/(math.pow(float(noAts),1./3.)) for i in range(3)]

flag=0
for st1 in range(steps):
    if flag:break
    for st2 in range(steps):
        if flag:break
        for st3 in range(steps):
            #for st in range(steps):
            atoms.append(Atom('%s'%atype,position=[st1*disp[0], st2*disp[1], st3*disp[2]]))

            if len(atoms)==noAts:flag=1;break

ase.io.write('geom-inp.dat',atoms,format='lammps-data',atom_style='charge')
        
"""
for at in range(1,noAts+1):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                atoms.append(Atom('%s'%atype,position=[ ]))
    print(atoms[-1])
"""

#print (math.pow(float(noAts),1./3.))
print(atoms)
view(atoms)
