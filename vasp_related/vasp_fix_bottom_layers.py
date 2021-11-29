#! /usr/bin/env python

import ase
import ase.io
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--input_poscar', default='POSCAR')
    parser.add_argument('-o', '--output_poscar', default='fixedPOSCAR.vasp')
    parser.add_argument('-zmax', '--zmax', type=float, default=0.5,
            help='z in direct coordinates below which the layers should be fixed')
    parser.add_argument('-between', '--between', type=float, nargs=2,
            help='z in direct coordinates below which the layers should be fixed')
    
    args = parser.parse_args()

    if args.between:
        zmin, zmax = args.between
    else:
        zmin = 0.0
        zmax = args.zmax

    atoms = ase.io.read(args.input_poscar)

    direct_coordinates = atoms.get_scaled_positions()

    atoms_to_be_fixed = []
    for i in range(len(direct_coordinates)):
        print direct_coordinates[i][2]
        if zmin < direct_coordinates[i][2] < zmax:
            atoms_to_be_fixed.append(i)

    c = ase.constraints.FixAtoms(atoms_to_be_fixed)
    c = atoms.set_constraint(c)

    ase.io.write(args.output_poscar, atoms, vasp5=True, direct=True)
