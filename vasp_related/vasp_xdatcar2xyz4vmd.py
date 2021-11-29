#! /usr/bin/env python

"""Functions for handling XDATCAR files

   read_xdatcar() returns an array of Atoms() objects with postions taken 
   from an XDATCAR file. 

   view_xdatcar() uses the ase gui to view the positions from an XDATCAR  
   file. 

   Author: Matthew Dyer    Date of last revision: 21 January 2009
"""

import os
import ase.io
import numpy as np
import argparse

def read_xdatcar(filename ,poscarfile):
    """Import XDATCAR type file.

    Reads unitcell, atom positions and constraints from the POSCAR/CONTCAR
    file and tries to read atom types from POSCAR/CONTCAR header, if this fails
    the atom types are read from OUTCAR or POTCAR file.

    Then reads in positions from an XDATCAR file and returns an array of 
    Atoms() objects with the positions from this file.
    """

    atoms=ase.io.read(poscarfile, format='vasp')
    natoms = len(atoms)
   
    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    data=f.readlines()

    nimages=(len(data)-5)/(natoms+1)

    images = []
    for i in range(nimages):
        images.append(atoms.copy())
        pos=np.zeros((natoms,3),np.float)
        for j in range(natoms):
            pos[j,0]=float(data[8+i*(natoms+1)+j].split()[0])
            pos[j,1]=float(data[8+i*(natoms+1)+j].split()[1])
            pos[j,2]=float(data[8+i*(natoms+1)+j].split()[2])
        images[i].set_scaled_positions(pos)

    return images

def xdatcar2xyz4vmd(filename='XDATCAR', poscarfile='CONTCAR'):
    configs = read_xdatcar(filename, poscarfile)
    configs = [config*(2,2,1) for config in configs]
    return configs

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description = """Script to convert XDATCAR to xyz format and
            multiply the suercell so that bonds form properly along the
            periodic boundaries""")

    parser.add_argument('-x', '--xdatcar', default='XDATCAR')
    parser.add_argument('-c', '--contcar', default='CONTCAR')
    parser.add_argument('-o', '--out_file', default='XDATCAR_for_vmd.xyz')

    args = parser.parse_args()
    
    configs = xdatcar2xyz4vmd(args.xdatcar, args.contcar)
    ase.io.write(args.out_file, configs, format='xyz')
