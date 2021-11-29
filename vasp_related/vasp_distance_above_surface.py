#! /usr/bin/env python

import sys
import ase.io
import numpy as np

def get_distance_above_surface(surfaceatoms, molecule):
    highest_surface_atom_pos = max([i.z for i in surfaceatoms])
    lowest_molecule_atom_pos = min([i.z for i in molecule])

    das = lowest_molecule_atom_pos - highest_surface_atom_pos

    return das

if __name__ == '__main__':
    try:
        fposcar = sys.argv[1]
    except IndexError:
        fposcar = 'POSCAR'

    print 'Using %s' % (fposcar)

    surface_plus_molecule = ase.io.read(fposcar, format='vasp')
    surfaceatoms = [i for i in surface_plus_molecule if i.symbol in ['Ru', 'O']]
    moleculeatoms = [i for i in surface_plus_molecule if i.symbol in ['C', 'H']]
    print 'Considering atoms in molecule as:', moleculeatoms

    print get_distance_above_surface(surfaceatoms, moleculeatoms)
