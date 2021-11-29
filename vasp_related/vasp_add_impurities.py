#!/usr/bin/env python
# encoding: utf-8

import sys
import ase.io
import numpy as np
import argparse


def add_impurity_at_random_position(atoms, impurity_atom):
    random_position = np.random.rand(1,3)
    cell = atoms.get_cell()
    impurity_atom.set_cell(cell)# = ase.Atoms(['H'], random_position, cell=cell)
    impurity_atom.set_scaled_positions(random_position)
    new_atoms = atoms.copy()
    new_atoms.extend(impurity_atom)
    return new_atoms

        
def get_atoms_with_impurity(atoms, impurity_atom):
    accept_new_position = False
    cnt = 0
    while not accept_new_position:
        new_atoms = add_impurity_at_random_position(atoms, impurity_atom)
        natoms = len(new_atoms)
        nimpurity = natoms - 1
        #print new_atoms[nimpurity]
        far_enough_from_all_atoms = (natoms-1)*[True]
        for i in range(natoms-1):
            rij = new_atoms.get_distance(nimpurity, i, mic=True)
            if rij < tolerance:
                far_enough_from_all_atoms[i] = False
                cnt += 1
                break

        if all(far_enough_from_all_atoms):
            print cnt
            accept_new_position = True

    return new_atoms

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add impurity atoms to a POSCAR file')

    parser.add_argument('-i', '--input-file', default='POSCAR')
    parser.add_argument('-o', '--output-file',
            default='POSCAR_with_impurities.vasp')
    parser.add_argument('-t', '--tolerance', type=float, default=2.0)
    parser.add_argument('-n', '--nimpurities', type=int, required=True)
    parser.add_argument('-m', '--impurity_molecule_file',
            default='molecule.xyz')
    args = parser.parse_args()

    tolerance = args.tolerance
    nimpurities = args.nimpurities
    impurity_atom = ase.io.read(args.impurity_molecule_file, format='xyz')

    if len(impurity_atom) > 1:
        print "Length of impurity atoms cannnot be > 1"
        sys.exit(-1)

    atoms = ase.io.read(args.input_file, format='vasp')
    for i in range(nimpurities):
        new_atoms = get_atoms_with_impurity(atoms, impurity_atom)
        atoms = new_atoms.copy()

    ase.io.write(args.output_file, new_atoms, format='vasp', vasp5=True, direct=True)
