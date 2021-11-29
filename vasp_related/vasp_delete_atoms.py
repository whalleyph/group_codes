#! /usr/bin/env python

import argparse
import ase.io

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--atoms', type=int, nargs='*',
            help='list of atoms to delete (count starts at 1)')

    parser.add_argument('-i', '--input_file', default='POSCAR')
    parser.add_argument('-o', '--output_file', default='POSCAR.deleted.vasp')

    args = parser.parse_args()

    atoms = ase.io.read(args.input_file, format='vasp')

    atoms_to_delete = [i-1 for i in args.atoms]
    print atoms_to_delete

    atoms.set_constraint()
    print 'Removed all constraints ... You will have to set the constraints again'
    del atoms[atoms_to_delete]

    ase.io.write(args.output_file, atoms, format='vasp', vasp5=True, direct=True)
