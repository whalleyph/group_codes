#! /usr/bin/env python

import ase.io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', default='POSCAR')
    parser.add_argument('-o', '--output_file')
    parser.add_argument('-r', '--translation-vector', type=float, nargs=3, required=True)

    args = parser.parse_args()

    if not args.output_file:
        output_file = args.input_file + '.translated.vasp'

    atoms = ase.io.read(args.input_file, format='vasp')

    atoms.translate(args.translation_vector)

    atoms.write(output_file, format='vasp', vasp5=True, direct=False)
