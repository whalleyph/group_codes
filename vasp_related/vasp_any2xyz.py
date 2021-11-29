#!/usr/bin/env python3
from __future__ import print_function    # (at top of module)

import numpy as np
import ase.io
import time , argparse
from copy import deepcopy as dc
from re import search,match
from sys import exit,argv
from os import system

def writeCoords(outf,coords): #Write coordinates in an xyz file (can be interfacewd with ASE for other file types.#coords is a list of Atoms objects for the whole trajectory, in coords[ts][atomid] format.
    noSteps=len(coords)
    noAtoms=len(coords[0])

    for st in range(noSteps):
        atoms=coords[st]
        lattice=atoms.get_cell()
        str2=""
        for lt in lattice: 
            for l in lt:
                str2+="%9.6f "%float(l)
        str1="%d\n"%noAtoms
        str1+='Lattice="%s" \n'%str2[0:-1]
        outf.write(str1)

        crd=atoms.get_positions(wrap=False)
        for at in range(noAtoms):
            outf.write("%3s %9.6f %9.6f %9.6f\n"%(atoms[at].symbol,crd[at][0],crd[at][1],crd[at][2]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inp", type=str, default=["XDATCAR"],nargs="*",
                        help="Input file to read the data from. Def: ./XDATCAR")
    parser.add_argument("-o","--out", type=str, #default="OUTCAR.xyz",
                        help="Output file name; a single file to concetanate trajectories from multiple input files. Def: ${inp}.xyz")
    parser.add_argument("-uw","--unwrap", action='store_true', default=False, help="To unwrap the atomic coordinates along the trajectory. Unwrapped coordinates are required for diffusivitiy analysis (e.g. MSD, VFAC, etc.). Def: Wrapped coordinates are saved.")
    parser.add_argument('-skip', '--skipSteps',type=int, default=1, help='To skip some time steps in the MD trajectory. Def: all frames are considered. ')

    #parser.add_argument("-l","--lattice", action='store_true', default=False, help="To include the lattice information in the output XYZ file. Def: No lattice info to speed-up the conversion.")
    #parser.add_argument("-e","--energy", action='store_true', default=False, help="To include the total (free) energy information in the output XYZ file. Def: No energy info.")

    args = parser.parse_args()

    initT=time.time()
    tsteps=[]
    if args.out: system("rm -f %s"%args.out); out=args.out;outf=open(out,"w")
    #init_pos=None # initial positions for referencing
    if args.skipSteps>1: print(('Skipping %d steps as requested...'%args.skipSteps))

    for inp in args.inp:
        #print "Working on %s. "%inp

        try:tsteps.extend([tr for tr in ase.io.iread(inp,index='::%d'%args.skipSteps)]) #Two takes almost same time.
        #tsteps=ase.io.read(inp,index=":")#read all steps.
        except: continue
            

    print("%d steps were read from %s in %.2f sec "%(len(tsteps),", ".join(args.inp),time.time()-initT))
    if not args.out:  out=inp.split(".")[0]+".xyz";  outf=open(out,"w")

    initT2=time.time()
    if args.unwrap:
        print("Unwrapping coordinates as requested.")
        noAtoms=len(tsteps[0])
        #if not init_pos: #done only once while the first file being read.
        init_pos=tsteps[0].get_scaled_positions()
        shifts=np.zeros((noAtoms,3))
        ppos=dc(init_pos)

        for st in range(1,len(tsteps)):#,args.skipSteps):
            atoms=tsteps[st]
            cpos=atoms.get_scaled_positions(wrap=False)
            #shifts=np.zeros((noAtoms,3))
            for at in range(noAtoms):
                for j in range(3):
                    diff=cpos[at][j]-ppos[at][j]
                    if np.fabs(diff) >0.5: #also covers the >1.0 differences.
                        shifts[at][j]+=(diff - np.sign(diff))
                    else: shifts[at][j]+=diff 

            tsteps[st].set_scaled_positions(init_pos+shifts)
            ppos=dc(cpos)
            
        print("Elapsed time: %.2f sec."%( time.time()-initT2))

    print("Writing/Appending coordinates to %s..."%args.out)
    if 1: writeCoords(outf,tsteps) #this works faster than built-in ase,io,write function for xyz format!!
    else: ase.io.write(outf,tsteps,format="xyz")#,append=True ,parallel=False) #can also work on file objects. Two options only work with ase >= 3.16.0

    if not args.out:outf.close()
    print("Elapsed time: %.2f sec."%( time.time()-initT2))

outf.close()
print("Total elapsed time: %.2f sec."%( time.time()-initT))
