#!/usr/bin/env python3
import ase.io
from ase.visualize import view
from os import system,popen,chdir,getcwd,getenv,putenv,listdir,environ
import os.path
from re import search
from sys import exit
#from ase.build import *
import ase.build
from ase import Atoms,Atom
import fractions
from copy import deepcopy
import time
#import ase.build.tools
import matplotlib.pyplot as plt
import ase.calculators.castep,ase.calculators.vasp
from ase.calculators.vasp import Vasp,Vasp2 #requires ase 3.17.0
from ase.calculators.espresso import Espresso

from spglib import find_primitive,standardize_cell,get_spacegroup #, niggli_reduce #niggli_reduce from ASE-tools conflicts with that from spglib.
from ase.build import niggli_reduce
import argparse
import os.path
from ase.spacegroup import Spacegroup,get_spacegroup

from subprocess import Popen,PIPE # as popen # check_output
from sys import stdout,version_info
def Popen4(cmd):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate()
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err

def callQE(fname,calc):

#    exe=getenv('VASP_COMMAND')
#    if not exe: exe='vasp_std'
#    print("VASP command: ", exe)

    environ['ASE_ESPRESSO_COMMAND']="/usr/local/Cluster-Apps/quantum-espresso/6.4/skl/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"

    cwd=os.getcwd()

    seed=fname.split('.')[0]
    try:chdir(seed)
    except:None

    print(cwd,fname)
    system('mv %s/%s .'%(cwd,fname))

    atoms=ase.io.read(fname)
    atoms.set_calculator(calc)
    Ep = atoms.get_potential_energy()
    
    #except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);#print Ep  
    #TODO:use popen4 instead

    chdir(cwd)

    return Ep,atoms


def call_vasp_v2(fname):

    exe=getenv('VASP_COMMAND')
    if not exe: exe='vasp_std'
    print("VASP command: ", exe)

    cwd=os.getcwd()

    seed=fname.split('.')[0]
    try:chdir(seed)
    except:None

    try:
        calc = Vasp(restart=True)
        atoms = calc.get_atoms()
        print ("VASP run was read succesfully from OUTCAR.")
        #else:
    except:
        print ("VASP run could not be read, starting a new run...")
        calc=Vasp2()
        #calc = Vasp2(atoms,restart=restart, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt=None)
        calc.read_incar()
        vasp_continue()
        atoms=ase.io.read("POSCAR",format="vasp")
        #try:atoms=ase.io.read("POSCAR",type="vasp",format="vasp4)
        calc.set(xc="pbe",ibrion=2)
        calc.directory="."#cdir
        atoms.set_calculator(calc)


    #try: 
    Ep = atoms.get_potential_energy()
    
    #except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);#print Ep  
    #TODO:use popen4 instead

    chdir(cwd)

    return Ep,atoms


def call_lammps(seed,relax_cmd='lammps_relax',lexec="mpirun -np 32  $HOME/APPS/lammps-Mar2020/lmp_openmpi", tmp_pp='C.pp',tmp_cell='C.cell',P=0.0 ): #uses the lammps_relax from airss.pl #TODO: use ASE-LAMMPS interface instead
    if not os.path.exists('%s.res.tmp'%seed):
        #system('ase-file-converter-bk.py -i %s.res -ot lammps-data -o %s.conf -ow ; mv %s.conf.data %s.conf'%(seed,seed,seed, seed))
        system('ase-file-converter-bk.py -i %s.res -ot castep-cell '%(seed))
        #system('echo "FIX_ALL_CELL : true" >> %s.cell'%(seed))
        with open( "%s.cell"%seed,'a') as f:
          f.write("FIX_ALL_CELL : true\n")
        system('mv %s.res %s.res.bkp'%(seed,seed))
        system('touch %s.res.tmp'%seed)
        system('cp %s %s.pp'%(tmp_pp,seed))
        system('cp %s.cell %s.cell'%(seed,seed.split('-')[0]))
        system('cp %s.pp %s.pp'%(seed,seed.split('-')[0]))
        #system('cp %s %s.cell'%(tmp_cell,seed)) #no need
        if not args.dryrun: 
          cmd='%s "%s" %s %s'%(relax_cmd,lexec,P,seed)
          system(cmd)
        system('castep2res 1 %s >%s-final.res'%(seed,seed))
        #exit()



parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                    help='an integer for the accumulator')


parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')
#parser.add_argument('-ot','--otype', type=str,required=True, help='Output file type')
#parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: input name is used as root.')
parser.add_argument('-it','--itype', type=str,required=False, help='Input file type. Def: determined automatically from the extension.')

#parser.add_argument('-t', '--tol',type=float, default=1e-4,help="The symmetry tolerance.")
parser.add_argument('-np', '--nprocs',type=int, default=32,help="No of processes to start for each LAMMPS calculation throuh mpirun. Def:32")

parser.add_argument('-dry','--dryrun', default=False,action='store_true', help='Do not perform the actual LAMMPS run but prepare the inputs and collect the results. Def: No')

#parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')

args = parser.parse_args()

initT=time.time()



cwd=os.getcwd()
cnt=0
for inpf in args.inpf:
    print (inpf)

    if not os.path.isfile(inpf):
        print('%s does not exists, skipping...'%inpf)
        continue

    #atoms=ase.io.read(inpf)
    seed=inpf.split('.')[0]

    exe="$HOME/APPS/lammps-Mar2020/lmp_openmpi"
    call_lammps(seed,relax_cmd='lammps_relax_sd',lexec="mpirun -np %d  %s"%(args.nprocs,exe))
