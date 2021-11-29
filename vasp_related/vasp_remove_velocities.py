#! /usr/bin/env python

"""
Script to remove velocities from a VASP POSCAR or CONTCAR file.
"""

import argparse
import ase.io

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='POSCAR')
    parser.add_argument('-o', '--outfile', default='newPOSCAR.vasp')
    args = parser.parse_args()

    atoms = ase.io.read(args.infile, format='vasp')
    ase.io.write(args.outfile, atoms, format='vasp', direct=True, vasp5=True)
