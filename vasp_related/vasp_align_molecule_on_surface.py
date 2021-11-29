#! /usr/bin/env python

import ase.io
import argparse

def align_molecule_on_surface(molecule='molecule.xyz',
        surface = 'POSCAR.surface',
        bond_in_molecule = None,
        bond_in_surface = None,
        center = None,
        center_atom = None):

    if isinstance(molecule, type('string')):
        molecule = ase.io.read(molecule, format='xyz')
    else:
        molecule = molecule

    if isinstance(surface, type('string')):
        surface = ase.io.read(surface, format='vasp')
    else:
        surface = surface

    # get positions
    molecule_positions = molecule.get_positions()
    surface_positions = surface.get_positions()

    # get the center point to be kept fixed on rotation
    if center:
        center = center
    elif center_atom:
        center = molecule_positions[center_atom-1]
    else:
        center = [0.0, 0.0, 0.0]

    # python starts count at 0, correct the indices
    bond_in_molecule = [i-1 for i in bond_in_molecule]
    bond_in_surface = [i-1 for i in bond_in_surface]

    # get the bonds on the molecule and the surface
    molecule_bond = molecule_positions[bond_in_molecule[1]] \
            - molecule_positions[bond_in_molecule[0]]
    surface_bond = surface_positions[bond_in_surface[1]] \
            - surface_positions[bond_in_surface[0]]

    # rotate the bond of the molecule to the one in the surface
    molecule.rotate(molecule_bond, surface_bond, center=center)

    return molecule


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="""script to align a bond of a molecule with another
            bond on the surface""")
    
    parser.add_argument('-m', '--molecule', default='molecule.xyz')
    parser.add_argument('-s', '--surface', default='POSCAR.surface')
    parser.add_argument('-mb', '--bond-in-molecule', nargs=2, type=int, default=[0,1])
    parser.add_argument('-sb', '--bond-in-surface', nargs=2, type=int, default=[0,1])
    parser.add_argument('-c', '--center', nargs=3, type=float, default=[0.0, 0.0, 0.0])
    parser.add_argument('-ca', '--center-atom', type=int, default=0)
    parser.add_argument('-o', '--out-file', default='molecule.rotated.xyz')

    args = parser.parse_args()

    molecule = align_molecule_on_surface(molecule = args.molecule,
            surface = args.surface,
            bond_in_molecule = args.bond_in_molecule,
            bond_in_surface = args.bond_in_surface,
            center = args.center,
            center_atom = args.center_atom)

    ase.io.write(args.out_file, molecule, format='xyz')
