#!/usr/bin/env python

import ase
import ase.io
import argparse
import numpy as np

def align_parallel_to_vector(atoms, n1, n2, vector, center):
    """aligns the bond containing the atoms n1, n2
    parallel the given vector
    """

    r = np.array([atoms.get_positions()[i-1] for i in [n1, n2]])

    r21 = r[1] - r[0]

    atoms.rotate(r21, vector, center=center)

    return atoms


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='atom count starts at 0')

    parser.add_argument('-m', '--file', help='molecule structure file',
            default='molecule.xyz')

    parser.add_argument('-n', '--atoms_in_bond',
            nargs=2, type=int, required=True,
            help='indices of atoms in the bond which will be put parallel \
to the given vector')

    parser.add_argument('-a', '--axis',
            help='axis parallel to the bond')

    parser.add_argument('-v', '--vector', nargs=3, type=float,
            help='vector parallel to the bond')

    parser.add_argument('-c', '--center', nargs=3, type=float,
            default=(0.0, 0.0, 0.0),
            help='center to be kept fixed')

    parser.add_argument('-ca', '--center-atom', type=int,
            default=0,
            help='center to be kept fixed')

    parser.add_argument('-o', '--outfile', help='output file', default='molecule.plane-aligned.xyz')

    args = parser.parse_args()

    n1, n2 = args.atoms_in_bond

    extension = args.file.split('.')[-1]
    if extension == 'xyz':
        fmt = 'xyz'
    elif extension == 'vasp':
        fmt = 'vasp'
    elif args.file in ['POSCAR', 'CONTCAR']:
        fmt = 'vasp'

    if args.axis:
        vector = args.axis
    elif args.vector:
        vector = args.vector
    else:
        print 'One of axis or vector should be specified'
        sys.exit(-1)

    atoms = ase.io.read(args.file, format=fmt)

    if args.center_atom:
        center_atom = args.center_atom - 1
        center = atoms.get_positions()[center_atom]
    else:
        center = args.center

    atoms_rotated = align_parallel_to_vector(atoms, n1, n2, vector, center=center)

    ase.io.write(args.outfile, atoms_rotated, format='xyz')
