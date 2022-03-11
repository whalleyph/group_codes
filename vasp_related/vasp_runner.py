#!/usr/bin/env python3
import ase.io
from ase.visualize import view
from os import system,popen,chdir,getcwd,getenv,putenv,listdir,environ
import os.path,glob
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
from spglib import find_primitive,standardize_cell,get_spacegroup #, niggli_reduce #niggli_reduce from ASE-tools conflicts with that from spglib.
from ase.build import niggli_reduce
import argparse

import os.path
from ase.phonons import Phonons

from random import uniform



from subprocess import Popen,PIPE # as popen # check_output
from sys import stdout,version_info
def Popen4(cmd,timeout=None):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate(timeout=timeout)
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err


def vasp_continue(ifDel=0):
        fl = listdir(".")
        fl.sort()
        #print fl
        flag=0
        ind=0
        for fn in fl:
                if fn[0:3]=="RUN" :
                        ind = int(fn[3:])

        run="RUN%d"%(ind+1)
        print ("Previous data stored in %s."%run)
        system('mkdir %s'%run)
        system('cp * %s/ 2> /dev/null '%run)
        system('cp CONTCAR POSCAR') #this was missing in the previous version !!
        if ifDel: system('rm WAVECAR CHGCAR')


def grep(key,fname,n=-1): 
    #Uses grep to get the nth line with the keyword from a given file. By default the last occurence is returned.
        try:
                return popen4('grep -m %d "%s" %s '%(n,key,fname),"r")[1].readlines()[-1][0:-1] #don't take \n at the end.  #Fastest one!!!
        except: 
                return popen('grep -m %d "%s" %s '%(n,key,fname),"r").read()[0:-1]
        else:
                return ""



def call_vasp_v2(fname,exe=None,xc='pbe'):

    if not exe: exe=getenv('VASP_COMMAND') ; #print ('bk')
    if not exe: exe='vasp_std'
    
    os.environ['VASP_COMMAND']=exe
    print("VASP command: ", exe)#, os.environ['VASP_COMMAND'])

    cwd=os.getcwd()

    seed=fname.split('.')[0]
    try:chdir(seed)
    except:None
    
    system('vasp_make_potcar -xc potpaw_PBE.54 &> /dev/null')

    #try:
    #       atoms = ase.io.read('vasprun.xml',index=0)
    #       print "vasprun.xml was read successfully."
    #except:

    flag=0 #whether to start a new/continuation run 
    try:
        #try:
        calc = Vasp(restart=True)
        #except:calc=None
        atoms = calc.get_atoms()
        print ("VASP run was read succesfully from OUTCAR.")
        if Vasp.read_convergence(calc): print('Geom opt already converged...')
        else:
          print('Geom opt not converged running a continuation job...')
          flag=1
        
        #else:
    except: 
        print ("VASP run could not be read, starting a new run...")
        flag=1

    if flag:
        calc=Vasp()
        #calc = Vasp2(atoms,restart=restart, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt=None)
        calc.read_incar(filename='INCAR')
        if os.path.exists("./OUTCAR"):   vasp_continue()
        atoms=ase.io.read("POSCAR",format="vasp")
        #try:atoms=ase.io.read("POSCAR",type="vasp",format="vasp4)
        calc.set(xc=xc)
        #setups='recommended'
        setups='minimal'
        calc.set(setups=setups)
        calc.directory="."#cdir
        atoms.set_calculator(calc)
        #calc.set(ldau_luj={'Si': {'L': 1, 'U': 3, 'J': 0}})



    #try: 
    Ep = atoms.get_potential_energy()
    
    #except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);#print Ep  
    #TODO:use popen4 instead

    chdir(cwd)

    return Ep,atoms


parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                    help='an integer for the accumulator')


parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')
parser.add_argument('-it','--itype', type=str,required=False, help='Input file type. Def: determined automatically from the extension.')

parser.add_argument('-t', '--tol',type=float, default=1e-4,help="The symmetry tolerance.")


parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')

parser.add_argument('-dry','--dryrun', default=False,action='store_true', help='Dry run: do not call VASP, just gather information from VASP run folders. Def: No')

parser.add_argument('-np', '--nprocs',type=int, default=32,help="No of processes to start for each VASP calculation throuh mpirun. Def:32")

parser.add_argument('-exe','--exe', type=str,required=False, default='vasp_std',help='Vasp exectuable. Def: vasp_std')

parser.add_argument('-xc','--xc', type=str,required=False, default='pbe',help='Exchange correlation functional to use (pseudopots will be selected accordingly, e.g. LDA, PBE,PBEsol, etc. Def: PBE')

parser.add_argument('-ph','--phonons', default=False,action='store_true',help='Do a phonon calculation if the geometry is converged.')
parser.add_argument('-ph_sc','--phonons_supercell', type=int,nargs=3,default=(2,2,2),help='Size of the supercell to be uused in the phonon calculation.Def: No supercell used, i.e. 1 1 1')
parser.add_argument('-ph_del', '--phonons_delta',type=float, default=0.03,help="Stepsize to use in the finite difference method to compute force constants. Def:0.03 A")

parser.add_argument('-ph_np','--phonons_noplot', default=False,action='store_true',help='Do a phonon calculation if the geometry is converged.')

parser.add_argument('-ph_path','--phonons_path', default=None,type=str,help='The high-symmetry Q-point path for a phonon calculation.')


args = parser.parse_args()

initT=time.time()

from ase.spacegroup import Spacegroup,get_spacegroup


exe="mpirun -np %d %s"%(args.nprocs,args.exe)
#exe="srun -n %d %s"%(args.nprocs,args.exe)


#########
if 0:
  files=Popen4("""grep "reached required accuracy - stopping structural energy minimisation" C-*/vasp.out -L | awk -F/ '{print $1}' """) [0]

  print('List of calculations with non-converged geometry: ', files)

  cwd=os.getcwd()
  for f in files:
    #os.chdir(f)
    #if os.path.exists("%s/%s-final.res"%(f,f)): system('cp %s/%s-final.res %s.res'%(f,f,f))
    #else:
      atoms=ase.io.read('%s/XDATCAR'%f, index=-1)
      ase.io.write('%s.res'%f,atoms,format='res')

    #os.chdir(cwd)
  #exit()

#######

system ('mkdir -p ./inputs')
cwd=os.getcwd()
cnt=0
import random
#for inpf in random.choices(args.inpf):
for inpf in args.inpf:
    print (inpf)

    if not os.path.isfile(inpf):
        print('%s does not exists, skipping...'%inpf)
        continue

    atoms=ase.io.read(inpf)
    seed=inpf.split('.')[0]
    if not args.dryrun:
      system('mkdir -p %s'%seed) 
      atoms.write('%s/POSCAR'%seed,format='vasp',vasp5=1)
      system('cp -f INCAR  POTCAR KPOINTS %s/'%(seed))
      system('cp %s inputs/; mv %s %s/'%(inpf,inpf,seed))

    chdir(seed)

    #if Popen4(""" grep 'reached required accuracy - stopping structural energy minimisation' """)[0]: 
    if not args.dryrun:
        time.sleep(uniform(0.5, 2.0)) #this to prevent overwriting when two seperate instances of vasp_runner starts running  in parallel in the same foelder at the same exact time.
        try: Ep,atoms=call_vasp_v2(inpf,exe,xc=args.xc)
        except Exception as err:print('%s: Problem with %s, skipping...'%(err,inpf));chdir(cwd);continue
    else: 
      try: Ep=atoms.info['energy']; 
      except Exception as err: print(err); Ep=0.0; #atoms=ase.io.read(inpf)
      if Ep==None: Ep=0.0


    atoms.info['energy']=Ep #atoms.get_total_energy()
    atoms.info['name']=inpf

    #ASE does not get the Pressure right when restarting
    #P=atoms.info.get('pressure')
    #if P is None: 
    P=float(Popen4("""grep pressure OUTCAR | tail -1 | awk '{print ($4+$9)*0.1}' """)[0][0]) #kB to GPa
    atoms.info['pressure']=P

#grep "free  energy   TOTEN  = " $seed/OUTCAR  | tail -1 | awk '{print " VASP: Final Enthalpy     = "$5}' >> $seed.castep

    SG=atoms.info.get('spacegroup')
    if SG is None: 

        try:SG=Popen4('symm %s'%seed,timeout=30)[0][0]  #part of the airss package #now only for cluster, also support the crystal mode.  # symm -cl
        except Exception as err: print(err); SG=''
        print (inpf,SG)#,get_spacegroup(atoms, symprec=args.tol))

    if SG is None or SG == '':   SG=str(get_spacegroup(atoms, symprec=args.tol).symbol.replace(' ',''))
    atoms.info['spacegroup']=SG.replace('(','').replace(')','')
    atoms.info['times_found']=1
    atoms.info['name']=seed
    
    print('%s: Ep=%.5f P= %.2f GPa SG=%s'%(inpf,Ep,P,atoms.info['spacegroup']))


    myres=ase.io.res.Res(atoms)
    myres.energy=Ep
    myres.write_file('%s.res'%(seed), write_info=0,  significant_figures=6) #works!

    #other ways
    #ase.io.res.write_res('%s.res'%(seed), atoms, write_info=0, write_results=1, significant_figures=6) #worksbut energy is written as None even though the value is stored.    
    #ase.io.write('%s-final.res'%(seed), atoms,format='res')


    #Do the phonons calculations (if requested by the user and if the geometry has converged)
    if args.phonons: 
        print('Doing phonons calculations as requested')
        calc=atoms.get_calculator()
        if not args.dryrun and not Vasp.read_convergence(calc): print('Geometry has not converged, skipping the phonons calculation...')
        else:
            
            if os.path.exists('phonon'): 
                print('Previous phonons calc was found in the %s/phonon directory, will read existing force data and do the missing points.'%seed)
                for f in glob.glob("phonon/*.json",recursive=0):
                    if os.path.getsize(f)==0: system('rm %s'%f); print (f,' is empty, deleted')#delete the empty files to recalculate

            calc.directory=calc.directory+"/Phonon-tmp/"
            print(calc.directory)
            calc.set(nsw=0,ibrion=-1,algo='fast') #do a single point at each phonon step
            calc.set(lwave=1) #save wavecar for restarting electronic wavefunction

            ph = Phonons(atoms, calc, supercell=args.phonons_supercell, delta=args.phonons_delta)
            if not args.dryrun: ph.run()
            else: print('Dry run, no phonon calculations wll be done for %s'%inpf)

            # Read forces and assemble the dynamical matrix
            method='Frederiksen' 
            method='standard'
            try: ph.read(acoustic=True,method=method) 
            except: print ('Problem while reading phonon calc from %s, skipping...'%seed);continue

            path = atoms.cell.bandpath(args.phonons_path, npoints=100) #if special_path is set to None,then it will be automatically determined.

            bs = ph.get_band_structure(path)
            dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)
            print(path)

            if args.phonons_noplot: chdir(cwd); continue

            # Plot the band structure and DOS:
            fig = plt.figure(1, figsize=(7, 4))
            ax = fig.add_axes([.12, .07, .67, .85])

            #emax = 0.035 
            emax=0.04
            #emax=max(dos.get_energies()) #TODO: get this from user
            bs.plot(ax=ax, emin=0.0, emax=emax)

            dosax = fig.add_axes([.8, .07, .17, .85])
            dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color='grey',
                               edgecolor='k', lw=1)

            dosax.set_ylim(0, emax)
            dosax.set_yticks([])
            dosax.set_xticks([])
            dosax.set_xlabel("DOS", fontsize=14)

            fig.savefig('%s_phonon.png'%seed)
            #plt.show()

            # Write modes for specific q-vector to trajectory files:
            #L = path.special_points['L']
            #ph.write_modes([l / 2 for l in L], branches=[2], repeat=(8, 8, 8), kT=3e-4,
            #   center=True)


            # Generate gif animation:
            # XXX Temporarily disabled due to matplotlib writer compatibility issue.
            # with Trajectory('phonon.mode.2.traj', 'r') as traj:
            #     write('Al_mode.gif', traj, interval=50,
            #           rotation='-36x,26.5y,-25z')


    chdir(cwd)
