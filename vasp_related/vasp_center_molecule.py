#! /usr/bin/env python

import ase.io
import argparse

# Move atoms to center of cell
def put_atom_at_center(atoms, center_atom = 1, side=10.0):

    cell = 3*[side]
    atoms.set_cell(cell)

    # Center of cell
    C = 0.5*(atoms.get_cell().sum(axis=1))
    AC = C - atoms.get_positions()[center_atom-1]
    atoms.translate(AC)
    return atoms

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='script to put a molecule in the center of a cubic box')
    
    parser.add_argument('-m', '--molecule', default='molecule.xyz',
            help='molecule in the xyz format')
    parser.add_argument('-n', '--index', type=int, default=1,
            help='index of the atom that needs to be at the \
            center of the box (count starts at 1)')
    parser.add_argument('-s', '--side', type=float, default=10.0,
            help='side of the box')
    parser.add_argument('-o', '--outfile')

    args = parser.parse_args()

    if args.outfile:
        outfile = args.outfile
    else:
        outfile = args.molecule + '.vasp'

    atoms = ase.io.read(args.molecule, format='xyz')

    centered_atoms = put_atom_at_center(atoms, args.index, args.side)

    ase.io.write(outfile, centered_atoms, format='vasp', sort=True, vasp5=True)
