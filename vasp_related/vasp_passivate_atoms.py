#! /usr/bin/env python

import ase.io
import ase
import argparse
import numpy as np

deg_to_rad = (np.pi/180.)

# atom hydrogen bond length
#ah_bl = 1.45

# bond_angle_between_hydrogens
bond_angle = deg_to_rad*108.
h_bond_angle = 0.5*bond_angle

r_h2 = np.array([
    [0.0, -np.cos(h_bond_angle), -np.sin(h_bond_angle)],
    [0.0,  np.cos(h_bond_angle), -np.sin(h_bond_angle)]
    ])

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""script to passivate a given list of atoms with
        one or two hydrogen atoms""")

    parser.add_argument('-i', '--surface', default='POSCAR',
            help='file with surface in POSCAR format')
    parser.add_argument('-o', '--output_file',
            default='POSCAR.passivated.vasp',
            help='file with surface in POSCAR format')
    parser.add_argument('-l1', '--list1', required=True, nargs='*', type=int,
            help='list of atoms to passivate with one H atoms (count starts at 1)')
    parser.add_argument('-l2', '--list2', nargs='*', type=int,
            help='list of atoms to passivate with two H atoms (count starts at 1)')
    parser.add_argument('-top', '--top', default=False, action='store_true',
            help='wheter to add the atoms on top or not')
    parser.add_argument('-bl', '--bond-length', default=1.4, type=float,
            help='bond length')
    args = parser.parse_args()

    atoms = ase.io.read(args.surface, format='vasp')
    positions = atoms.get_positions()

    ah_bl = args.bond_length

    if args.top:
        r_h1 = np.array([0.0, 0.0, ah_bl])
    else:
        r_h1 = np.array([0.0, 0.0, -ah_bl])

    # still have to deal with the part when two hydrogen atoms are being added

    print args.list1
    hydrogens_to_add = []
    for a in args.list1:
        r = positions[a-1]
        # position of h1 in_place
        rh1_ip = r + r_h1
        h1_ip = ase.atom.Atom(symbol='H', position=rh1_ip)
        hydrogens_to_add.append(h1_ip)

    hydrogens_to_add = ase.Atoms(hydrogens_to_add, cell=atoms.get_cell())

    atoms.extend(hydrogens_to_add)

    ase.io.write(args.output_file, atoms, vasp5=True, direct=True)
