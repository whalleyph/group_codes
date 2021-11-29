#!/usr/bin/env python2
import numpy as np
import ase.io
import time 
from copy import deepcopy as dc

initT=time.time()
inpf="all-traj.xyz"
tsteps=[tr for tr in ase.io.iread(inpf)]

print "%d steps were read from %s"%(len(tsteps),inpf)

init_pos=tsteps[0].get_scaled_positions()
ppos=dc(init_pos)
noAtoms=len(tsteps[0])
shifts=np.zeros((noAtoms,3))

for st in range(1,len(tsteps)):
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

print "Writing unwrapped.xyz..."
ase.io.write("unwrapped.xyz",tsteps,format="xyz")

print "Elapsed time: %.2f sec."%( time.time()-initT)
