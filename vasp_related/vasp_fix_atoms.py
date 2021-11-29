#! /usr/bin/env python

import sys
import ase.io
import ase
import argparse

parser = argparse.ArgumentParser(
        description='')

parser.add_argument('-i', '--input_file', default='POSCAR')
parser.add_argument('-o', '--output_file', default='POSCAR.fixed.vasp')
parser.add_argument('-n', '--natoms', nargs='*', type=int,
        help='list of atoms to be fixed (count starts at 1)')
parser.add_argument('-t', '--atom_types',
        help='Atom types that have to be fixed')

args = parser.parse_args()

atoms = ase.io.read(args.input_file, format='vasp')

if args.natoms:
    atoms_to_be_fixed = [i-1 for i in args.natoms]
elif args.atom_types:
    atoms_to_be_fixed = [atom.index for atom in atoms if atom.symbol in args.atom_types]
else:
    print 'Either atom_types or natoms have to be specified'
    sys.exit(-1)

c = ase.constraints.FixAtoms(indices=atoms_to_be_fixed)
atoms.set_constraint(c)

ase.io.write(args.output_file, atoms, format='vasp', vasp5=True, direct=True)
