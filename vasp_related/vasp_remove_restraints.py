#! /usr/bin/env python

import ase.io
import ase
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-f', '--files', nargs='*')
    
    args = parser.parse_args()

    for f in args.files:
        atoms = ase.io.read(f, format='vasp')
        atoms.set_constraint()
        ase.io.write(f, atoms, format='vasp', vasp5=True, direct=True)
