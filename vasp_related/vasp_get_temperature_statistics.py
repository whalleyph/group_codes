#! /usr/bin/env python

import subprocess as sp
import numpy as np
import ase.io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='-t Target Temparature \n-i initial time step to start the statistics')
    
    parser.add_argument('-t', '--reqtemp', type=float, required=True)
    parser.add_argument('-i', '--initstep', type=int, required=False)
   
    args = parser.parse_args()
    print Exception

    req_temp = args.reqtemp
    try:init_step = args.initstep
    except: init_step=0

    temp = sp.check_output("F | awk '{print $3}'", shell=True)

    temp = [float(i) for i in temp.strip().split()]
    temp = temp[init_step:]
    temp = np.array(temp)

    natoms = len(ase.io.read('POSCAR'))

    print '%30s %i' % ('No. of atoms:', natoms)
    print '%30s %i' % ('No. of time steps considered:', len(temp))
    print '%30s %.2f' % ('Required temperature:', req_temp)
    print '%30s %.2f' % ('Required sqrt. variance:', np.sqrt(2*req_temp**2/(3*natoms)))

    print '%30s %.2f' % ('Mean temperature:', np.mean(temp))
    print '%30s %.2f' % ('Sqrt. variance:', np.sqrt(np.var(temp)))
