#! /usr/bin/env python

"""

Script to find undercoordinated atoms

"""

import ase.io
import argparse
from vasp_get_nearest_neighbors import get_nearest_neighbors

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', default='POSCAR')
    parser.add_argument('-a', '--atom_type', required=True,
            help='atom type which should be checked for undercoordination')
    parser.add_argument('-c', '--coordination', required=True, type=int,
            help='coordination number that should be looked for')
    parser.add_argument('-z', '--z_min_max', nargs=2, type=float, default=[0.0, 1.0],
            help='zmin and zmax between which undercoordination should be searched')
    parser.add_argument('-rc', '--rcut', type=float, required=True,
            help='cut off radius within which the coordinated atoms should be searched')

    args = parser.parse_args()

    atoms = ase.io.read(args.input_file, format='vasp')
    zmin, zmax = args.z_min_max
    req_atom_type = args.atom_type
    rc = args.rcut
    ncoord = args.coordination

    print 'Looking for coordination of %s' % (req_atom_type)
    print 'Looking for coordination number of %i within radius %.2f' % (ncoord, rc)
    print 'Looking for a layer of atoms between %.2f (zmin) and %.2f (zmax)' % (zmin, zmax)

    direct_positions = atoms.get_scaled_positions()

    layer_atoms = []
    for atom, r in zip(atoms, atoms.get_scaled_positions()):
        if zmin < r[2] < zmax and atom.symbol == req_atom_type:
            layer_atoms.append(atom.index)

    print 'atom indices of atoms in the layer (count starts at 0)'
    print layer_atoms

    for req_atom in layer_atoms:
        if len(get_nearest_neighbors(req_atom, atoms, rc)[1:]) == ncoord:
            print req_atom+1,
