#! /usr/bin/env python

"""Script to replicate the same molecule multiple times in a given supercell"""

import argparse
import ase
import ase.io
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--molecule_file', default='molecule.xyz')
    parser.add_argument('-p', '--poscar_file', default='POSCAR')
    parser.add_argument('-l', '--list', nargs='*', type=float, default=[0.0, 0.0, 0.0],
            help = 'direct coordinates where the molecule should be placed')
    parser.add_argument('-n', '--ref_atom_in_molecule', type=int, default=1,
            help = 'count starts at 1')
    parser.add_argument('-o', '--output_file', default='layer.xyz',
            help = 'count starts at 1')

    args = parser.parse_args()
    ref_atom = args.ref_atom_in_molecule - 1
    coordinates = args.list

    if len(coordinates)%3 != 0:
        print 'Error'
    else:
        coord_list = []
        while coordinates:
            coord_list.append(coordinates[:3])
            coordinates = coordinates[3:]

    coord_list = np.array(coord_list)

    print coord_list
    # translate molecule to origin
    molecule = ase.io.read(args.molecule_file, format='xyz')
    positions = molecule.get_positions()
    pos_ref_atom = positions[ref_atom]
    molecule.translate(-pos_ref_atom)

    # set cell
    cell = ase.io.read(args.poscar_file, format='vasp').get_cell()
    molecule.set_cell(cell)

    print molecule.get_positions()
    print molecule

    layer = ase.Atoms()
    layer.set_cell(cell)
    for direct_coordinates in coord_list:
        new_molecule = molecule.copy()
        print 'cell', cell
        print 'direct', direct_coordinates
        cartesian_coordinates = np.dot(cell.transpose(), direct_coordinates)
        print 'cart', cartesian_coordinates
        new_molecule.translate(cartesian_coordinates)
        print new_molecule
        layer.extend(new_molecule)

ase.io.write(args.output_file, layer, format='xyz')
