#! /usr/bin/env python

import numpy as np
import sys
import ase.io
import argparse

def shift_origin(atoms, x, y, z):
    pos = atoms.get_scaled_positions()
    newpos = pos - np.array([x, y, z])
    atoms.set_scaled_positions(newpos)
    newpos = atoms.get_scaled_positions()
    atoms.set_scaled_positions(newpos)
    return atoms

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='Shifts origin by a given offset')

    parser.add_argument('-s', '--shift', nargs=3, type=float,
            default=[0., 0., 0.],
            help = 'offsets in scaled coordinates')
    parser.add_argument('-i', '--input-file', default='POSCAR',
        help = 'input file')
    parser.add_argument('-c', '--output-cartesian',
        action = 'store_true', default = False,
        help = 'if output should be in cartesian format')
    parser.add_argument('-o', '--output-file', default='POSCAR_shifted.vasp',
        help = 'output file')
    
    args = parser.parse_args()

    x, y, z= args.shift
    iposcar = args.input_file
    oposcar = args.output_file
    output_cartesian = args.output_cartesian
    
    atoms = ase.io.read(iposcar, format='vasp')
    atoms = shift_origin(atoms, x, y, z)

    if output_cartesian:
        ase.io.write(oposcar, atoms, format='vasp', vasp5=True)#, direct=True)
    else:
        ase.io.write(oposcar, atoms, format='vasp', vasp5=True, direct=True)
