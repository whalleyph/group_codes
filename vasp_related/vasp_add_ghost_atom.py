#! /usr/bin/env python

import ase
import ase.io
import argparse
import numpy as np

def get_ghost_coordinates(atoms, indices):
    """returns the coordinates of the ghost atom depending on the number 
    of the indices"""

    r = np.array([atoms.get_positions()[i-1] for i in indices])

    ghost_coordinates = np.average(r, axis=0)

    return ghost_coordinates

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', '--molecule', help='molecule structure file', default='molecule.xyz')
    parser.add_argument('-n', '--indices', nargs='*', type=int,
            help='indices of atoms at whose center the ghost atom \
will be added; index count starts at 1')
    parser.add_argument('-o', '--outfile', help='output file', default='molecule_with_ghost.xyz')

    args = parser.parse_args()

    atoms = ase.io.read(args.molecule, format='xyz')

    if args.indices:
        indices = args.indices
        natoms = len(indices)

        gx, gy, gz = get_ghost_coordinates(atoms, indices)

        ghost_atom = ase.Atoms('X', [(gx, gy, gz)])
        atoms_with_ghost = ghost_atom + atoms

    else:
        ghost_atom = ase.Atoms('X', [(0.0, 0.0, 0.0)])
        atoms_with_ghost = ghost_atom + atoms

    ase.io.write(args.outfile, atoms_with_ghost, format='xyz')
