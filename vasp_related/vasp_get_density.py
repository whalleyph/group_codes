#! /usr/bin/env python

import ase.io
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--input-file', default='POSCAR')
    
    args = parser.parse_args()

    input_file = args.input_file

    atoms = ase.io.read(input_file, format='vasp')

    volume = atoms.get_volume()
    mass = np.sum(atoms.get_masses())   # mass in gm/mole  
    mass = mass/6.022e23                # mass in gm 

    density = mass*1e24/volume          # density in gm/cm^3

    print '%.3f gm/cm^3' % (density)
