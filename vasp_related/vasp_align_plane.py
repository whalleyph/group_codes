#!/usr/bin/env python

import ase
import ase.io
import argparse
import numpy as np

def align_parallel_to_plane(atoms, n1, n2, n3, plane, center):
    """aligns the plane containing the atoms n1, n2, n3
    parallel the given plane
            1

                2

            3
    """

    r = np.array([atoms.get_positions()[i-1] for i in [n1, n2, n3]])

    r21 = r[1] - r[0]
    r23 = r[2] - r[1]

    atoms_plane = np.cross(r21, r23)

    atoms.rotate(atoms_plane, plane, center=center)

    return atoms


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', '--file', help='molecule structure file',
            default='molecule.xyz')

    parser.add_argument('-n', '--atoms_in_plane', nargs=3, type=int,
            required=True, help='indices of atoms in plane which will be \
placed parallel the given plane')

    parser.add_argument('-c', '--center', nargs=3, type=float,
            default=(0.0, 0.0, 0.0),
            help='center to be kept fixed')

    parser.add_argument('-ca', '--center-atom', type=int,
            default=0,
            help='center to be kept fixed')

    parser.add_argument('-p', '--plane', help='vector perpendicular \
to the required plane')

    parser.add_argument('-o', '--outfile', help='output file', default='molecule.plane-aligned.xyz')

    args = parser.parse_args()

    n1, n2, n3 = args.atoms_in_plane

    extension = args.file.split('.')[-1]
    if extension == 'xyz':
        fmt = 'xyz'
    elif extension == 'vasp':
        fmt = 'vasp'
    elif args.file in ['POSCAR', 'CONTCAR']:
        fmt = 'vasp'

    atoms = ase.io.read(args.file, format=fmt)

    if args.center_atom:
        center_atom = args.center_atom - 1
        center = atoms.get_positions()[center_atom]
    else:
        center = args.center

    atoms_rotated = align_parallel_to_plane(atoms, n1, n2, n3, 'z', center=args.center)

    ase.io.write(args.outfile, atoms_rotated, format='xyz')
