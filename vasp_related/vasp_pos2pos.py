#! /usr/bin/env python

"""
Script to convert POSCAR file from
    direct coordinates -> xyz
    xyz -> direct coordinates
"""

import ase.io
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input-files', nargs='*', default=['POSCAR'])
    parser.add_argument('-o', '--output-file')
    parser.add_argument('-x2d', '--xyz2direct', action='store_true')
    parser.add_argument('-d2x', '--direct2xyz', action='store_true')
    parser.add_argument('-v4', '--vasp4format', action='store_false')

    args = parser.parse_args()

    for input_file in args.input_files:
        atoms = ase.io.read(input_file, format='vasp')

        if args.xyz2direct:
            if not args.output_file:
                output_file = input_file + 'direct.vasp'
            else:
                output_file = args.output_file
            ase.io.write(output_file, atoms, format='vasp', vasp5=args.vasp4format, direct=True)

        if args.direct2xyz:
            if not args.output_file:
                output_file = input_file + 'xyz.vasp'
            else:
                output_file = args.output_file
            ase.io.write(output_file, atoms, format='vasp', vasp5=args.vasp4format, direct=False)
