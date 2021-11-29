#! /usr/bin/env python

import ase
import ase.io
import ase.lattice.surface as als
import numpy as np
import argparse

def get_cartesian_xy(surface, ab):
    ghost = ase.Atoms('X', [[0,0,0]], cell=surface.get_cell(), pbc=surface.get_pbc())
    x, y = ab
    z = 0.0
    ghost.set_scaled_positions([x, y, z])
    xy = ghost.get_positions()[0][0:2]
    #print surface.get_cell()

    return xy

def add_molecule_to_surface(surface, molecule, xy=(0.0, 0.0), das=2.5, mol_index=1, 
        direct=False, delete_ghost_atom=True, top_atom=None):

    if isinstance(surface, type('string')):
        surface = ase.io.read(surface, format='vasp')

    if isinstance(molecule, type('string')):
        molecule = ase.io.read(molecule, format='xyz')

    #print surface
    #print molecule

    surf = surface.copy()
    mol = molecule.copy()

    if top_atom:
        print top_atom
        top_atom = top_atom - 1
        surf.adsorbate_info = {}
        surf.adsorbate_info['top layer atom index'] = top_atom
        atom = surface[top_atom]
        xy = (atom.x, atom.y)
        print 'Adding fragment on top of %i atom' % (top_atom+1)
    elif direct:
        # Get the xy-coordinates for the direct coordinates in the unit cell 
        # of the original surface POSCAR
        xy = get_cartesian_xy(surface, xy)

    #print das
    #print surf.adsorbate_info

    else:
        print 'Using default xy', xy

    #print xy
    als.add_adsorbate(surf, mol, das, position=xy, mol_index=mol_index-1)

    # delete ghost atoms 'X'
    if delete_ghost_atom:
        if 'X' in surf.get_chemical_symbols():
            atoms_to_delete = [atom.index for atom in surf if atom.symbol == 'X']
            if len(atoms_to_delete) > 1:
                del surf[atoms_to_delete]
            else:
                del surf[atoms_to_delete[0]]
    return surf

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--surface', default='POSCAR.surface')
    parser.add_argument('-m', '--molecule', default='molecule.xyz')

    parser.add_argument('-n', '--natoms', type=int, nargs='*',
            help='atom index on which the molecule has to be placed \
(index starts at 1)')

    parser.add_argument('-nm', '--natom-molecule', type=int, default=1,
            help='atom of molecule which should be placed on the adsorption site\
(index starts at 1)')

    parser.add_argument('-xy', '--xy', nargs=2, type=float,
            help='x and y cartesian coordinates on the surface')

    parser.add_argument('-ab', '--ab', nargs=2, type=float,
            help='x and y direct coordinates on the surface')

    parser.add_argument('-das', '--das', type=float,
            help='distance above the surface', required=True)

    parser.add_argument('-o', '--outfile',
            help='filename of output file', default='POSCAR')

    args = parser.parse_args()

    surface = ase.io.read(args.surface, format='vasp')
    molecule = ase.io.read(args.molecule, format='xyz')
    das = args.das
    if args.natoms:
        natoms = [i-1 for i in args.natoms]
    else:
        natoms = None

    surface_with_molecule = surface.copy()
    if natoms:
        for iatom in natoms:
            surface_with_molecule.adsorbate_info = {}
            surface_with_molecule.adsorbate_info['top layer atom index'] = iatom
            atom = surface_with_molecule[iatom]
            xy = (atom.x, atom.y)
            print args.das
            surface_with_molecule = add_molecule_to_surface(surface_with_molecule, molecule, xy, args.das, args.natom_molecule)
    elif args.ab:
        # Get the coordinates of a ghost atom in the unit cell of the original surface POSCAR
        xy = get_cartesian_xy(surface, args.ab)
        print xy
    else:
        xy = args.xy

    surface_with_molecule = add_molecule_to_surface(surface_with_molecule, molecule, xy, args.das, args.natom_molecule)

        # delete ghost atoms 'X'
        #if 'X' in surface_with_molecule.get_chemical_symbols():
        #    del surface_with_molecule[[atom.index for atom in surface_with_molecule if atom.symbol == 'X']]

    ase.io.write(args.outfile, surface_with_molecule, format='vasp', vasp5=True, direct=True)
