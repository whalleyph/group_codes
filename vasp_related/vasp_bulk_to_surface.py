#! /usr/bin/env python

import ase
import ase.io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-file', default='POSCAR.bulk.vasp',
            help='bulk crystal structure')
    parser.add_argument('-v', '--vacuum', default=7.5, type=float)
    parser.add_argument('-o', '--out-file', default='POSCAR.surface.vasp')
    parser.add_argument('-c', '--constraints', action='store_true', default=False)
    parser.add_argument('-n', '--natoms', nargs='*', type=int,
            help='list of atoms to be fixed (count starts at 1)')

    args = parser.parse_args()
    atoms = ase.io.read(args.in_file)
    atoms.center(vacuum=args.vacuum, axis=2)
    if args.constraints:
        atoms_to_be_fixed = [i-1 for i in args.natoms]
        c = ase.constraints.FixAtoms(indices=atoms_to_be_fixed)
        atoms.set_constraint(c)
    ase.io.write(args.out_file, atoms, direct=False, vasp5=True)
