#! /usr/bin/env python

import ase
import ase.io
import numpy as np
import argparse
from vasp_center_molecule import put_atom_at_center

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description="""script to run a calculation to get convergence with resepect to box size""")

    parser.add_argument('-m', '--molecule', default='molecule.xyz')
    parser.add_argument('-l', '--length-range', nargs=2, type=float, default=[10.0, 20.0])
    parser.add_argument('-dl', '--length-increment', type=float, default=1.0)
    parser.add_argument('-k', '--keep-files', nargs='*', default=['vasp.out'])
    parser.add_argument('-ca', '--center-atom', type=int, default=1)

    args = parser.parse_args()

    l = args.length_range
    dl = args.length_increment
    ca = args.center_atom
    
    atoms = ase.io.read(args.molecule)
    box_lengths = np.arange(l[0], l[1], dl)

    for i in box_lengths:
        file_tag = '%4.2f' % (i)
        atoms_in_box = put_atom_at_center(atoms, ca, i)
        ase.io.write('.'.join(['POSCAR', file_tag, 'vasp']),
                atoms_in_box, format='vasp', vasp5=True, direct=True)
