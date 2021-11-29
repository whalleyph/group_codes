#! /usr/bin/env python

import ase.io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='POSCAR')
    parser.add_argument('-o', '--outfile')
    
    args = parser.parse_args()

    infile = args.infile
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = infile + '.xyz'

    atoms = ase.io.read(infile, format='vasp')
    ase.io.write(outfile, atoms, format='xyz')
