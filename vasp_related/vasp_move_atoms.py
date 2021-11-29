#! /usr/bin/env python

import ase.io
import ase
import argparse
import numpy as np

def get_cartesian_coordinates(direct_coordinates, cell):
    atoms = ase.Atoms('X', positions=[(0.0, 0.0, 0.0)], cell=cell)
    atoms.set_scaled_positions(direct_coordinates)
    return atoms.get_positions()[0]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-sm', '--surface_and_molecule', default='POSCAR',
            help='surface file')
    parser.add_argument('-a', '--iatoms_to_move',
            help='group of atoms to move (count starts at 1)', type=int, nargs='*')
    parser.add_argument('-r', '--direct_coordinates_to_move_to', type=float, nargs=3, required=True)
    parser.add_argument('-i', '--ireference_atom', type=int, default=0)
    parser.add_argument('-o', '--output_file', default='movedPOSCAR.vasp')
    
    args = parser.parse_args()

    # make two copies of all the atoms
    # 1. atoms_to_move
    # 2. remaining_atoms = all_atoms - atoms_to_move
    all_atoms_1= ase.io.read(args.surface_and_molecule, format='vasp')
    all_atoms_2 = all_atoms_1.copy()
    all_atoms_3 = all_atoms_1.copy()

    iatoms_to_move = args.iatoms_to_move
    iatoms_to_move = np.array(iatoms_to_move) - 1

    # get atoms_to_move
    # delete all other atoms except atoms_to_move from surface + molecule
    all_other_atoms = list(set(range(len(all_atoms_1))) - set(iatoms_to_move))
    all_atoms_1.set_constraint()
    del all_atoms_1[all_other_atoms]
    atoms_to_move = all_atoms_1
    #print atoms_to_move

    # get remaining_atoms
    # delete atoms_to_move from all_atoms
    all_atoms_2.set_constraint()
    del all_atoms_2[iatoms_to_move]
    remaining_atoms = all_atoms_2
    #print remaining_atoms

    # get coordinates of the reference atom
    ax, ay, az = all_atoms_3.get_positions()[args.ireference_atom-1]

    # get coordinates of the point to translate to in cartesian coordinates
    bx, by, bz = get_cartesian_coordinates(args.direct_coordinates_to_move_to,
            all_atoms_3.get_cell())

    # get translation vector
    tx, ty, tz = bx-ax, by-ay, bz-az
    
    #print ax, ay, az
    #print bx, by, bz
    #print tx, ty, tz

    # move the atoms_to_be_moved
    #ase.io.write('unmoved.vasp', atoms_to_move, format='vasp', vasp5=True)
    atoms_to_move.translate([tx, ty, tz])
    #ase.io.write('moved.vasp', atoms_to_move, format='vasp', vasp5=True)

    remaining_atoms.extend(atoms_to_move)

    #print remaining_atoms
    ase.io.write(args.output_file, remaining_atoms, format='vasp', vasp5=True, direct=True)
