#! /usr/bin/env python

import sys
import argparse
from pyspglib import spglib
import ase.io


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tolerance', type=float, default=1e-5)

    parser.add_argument('-p', '--poscar', default='POSCAR')

    args = parser.parse_args()

    atoms = ase.io.read(args.poscar, format='vasp')

    print 'Space group: %s' % (spglib.get_spacegroup(atoms, symprec=args.tolerance))

    #print spglib.get_symmetry_dataset(atoms, tolerance=args.tolerance)

    print 'Primitive unitcell'
    print spglib.find_primitive(atoms)
