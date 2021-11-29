#! /usr/bin/env python

import ase
import ase.io
import numpy as np
import argparse

def set_up_orientations(molecule, bond_atoms=[1,2], center_atom=0,
        center=[0., 0., 0.]):

    orientations = 'x -x y -y z -z'.split()

    for o in orientations:
        mol = molecule.copy()
        pos = mol.get_positions()
        # bond to be aligned
        ab = pos[bond_atoms[1]-1] - pos[bond_atoms[0]-1]
        if center_atom:
            mol.rotate(ab, o, center=pos[center_atom])
        elif center:
            mol.rotate(ab, o, center=center)
        else:
            mol.rotate(ab, o)
        # remove x
        #del mol[[i.index for i in mol if i.symbol == 'X']]
        #als.add_adsorbate(surf, mol, 4.5, position=(x,y), mol_index=0)
        ase.io.write('molecule.' + o.replace('-', 'm') + '.xyz', mol, format='xyz')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="""- set up a molecule in 6 different orientations
            - align the given bond along x, y, z -x, -y, -z directions""")
    
    parser.add_argument('-m', '--molecule', default='molecule.aligned.xyz',
            help='molecule')
    parser.add_argument('-n', '--bond-atoms', type=int, nargs=2, default=[1,2],
            help='bond to be aligned')
    parser.add_argument('-c', '--center-atom', type=int, default=1,
            help='center atom to be kept fixed during rotation')
    parser.add_argument('-cr', '--center', nargs=3, type=float, default=[0,0,0],
            help='cartesian coordinates of the point to be kept fixed during rotation')
    
    args = parser.parse_args()
    molecule = ase.io.read(args.molecule)
    bond_atoms = args.bond_atoms
    center_atom = args.center_atom
    center = args.center

    set_up_orientations(molecule, bond_atoms, center_atom, center)
