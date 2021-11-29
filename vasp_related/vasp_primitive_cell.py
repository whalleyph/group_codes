#! /usr/bin/env python

import ase.io
from spglib import *
import argparse

from os import system
from sys import exit
from copy import deepcopy as dc
from ase import Atoms

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='POSCAR')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-t', '--tol',type=float, default=1e-4,help="The symmetry tolerance.")
    
    args = parser.parse_args()

    infile = args.infile
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = infile + '_prim.vasp'
    #print atoms
    atoms = ase.io.read(infile, format='vasp')
    scaled_positions= atoms.get_scaled_positions()#(wrap=True) #atoms.scaled_positions # does not exist
    #print scaled_positions

    if 1: #using SPGlib directly.
        cell=(atoms.cell, scaled_positions, atoms.numbers)
        #print cell
        print "Space group of the given cell using tolerance=%f: %s"%(args.tol,get_spacegroup(cell,symprec=args.tol))
        lattice, scaled_positions, numbers = find_primitive(cell, symprec=args.tol)
        #lattice, scaled_positions, numbers = standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=args.tol)

        atoms2=Atoms(symbols=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
        ase.io.write(outfile, atoms2, format='vasp',direct=True)

    #exit()
    else: #Using cellsym by RJM.
        ase.io.write(infile + '_prim.cell', atoms, format='castep-cell')
    

        system("cellsym --cell_abs --primitive -e=%f %s_prim.cell out.cell"%(args.tol,infile))
        atoms = ase.io.read("out.cell", format='castep-cell')
        ase.io.write(outfile, atoms, format='vasp')
