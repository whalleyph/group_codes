#! /usr/bin/env python

import ase
import ase.io
from vasp_get_nearest_neighbors import get_nearest_neighbors
import argparse
import numpy as np

def get_undercoordinated_atoms(atoms, atom_type=None,
        coordination_number=None, rcut=None):
    req_atoms = [atom.index for atom in atoms if atom.symbol == atom_type]

    undercoordinated_atoms = []
    for i in req_atoms:
        nn = get_nearest_neighbors(i, atoms, rcut=rcut)
        #print i, len(nn)
        if len(nn) < coordination_number+2:
            undercoordinated_atoms.append(i)

    return undercoordinated_atoms


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='POSCAR')
    parser.add_argument('-at', '--atom_type', required=True)
    parser.add_argument('-rc', '--rcut', required=True, type=float)
    parser.add_argument('-nc', '--coordination_number', required=True, type=int)
    parser.add_argument('-m', '--mark', action='store_true')
    parser.add_argument('-o', '--outfile', default='marked.vasp')
    
    args = parser.parse_args()
    infile = args.infile
    atom_type = args.atom_type
    rc = args.rcut
    nc = args.coordination_number

    atoms = ase.io.read(infile)

    undercoordinated_atoms = get_undercoordinated_atoms(atoms, atom_type, nc, rc)
    print np.array(undercoordinated_atoms) + 1

    if args.mark:
        marked_atoms = atoms.copy()
        for iatom in undercoordinated_atoms:
            marked_atoms[iatom].symbol = 'X'
        ase.io.write(args.outfile, marked_atoms, format='vasp', vasp5=True, direct=True)


    #ase.io.write('undercoord.vasp', atoms, format='vasp', vasp5=True, direct=False)
