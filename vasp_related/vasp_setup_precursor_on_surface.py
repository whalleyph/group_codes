#! /usr/bin/env python

import ase
import ase.io
import argparse
from vasp_add_molecule_to_surface import add_molecule_to_surface
from vasp_add_molecule_to_surface import get_cartesian_xy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-m', '--molecule', default='precursor.xyz',
            help='precursor molecule to setup on the surface\n
atom 1 will be placed on top of the given atom')
    parser.add_argument('-s', '--surface', default='POSCAR.surface',
            help='POSCAR of the surface on which to set up the precursor')
    parser.add_argument('-n', '--top-surface-atom', type=int,
            help='atom on top of which the precursor should be positioned\n
atom count starts at 1')
    parser.add_argument('-ab', '--direct', nargs=2, type=float,
            help='direct coordinates at which the molecule should be placed')
    
    args = parser.parse_args()

    molecule = ase.io.read(args.molecule, format='xyz')
    surface = ase.io.read(args.surface, format='vasp')

    if args.top_surface_atom:
        top_surface_atom = args.top_surface_atom - 1
        surface.adsorbate_info = {}
        surface.adsorbate_info['top layer atom index'] = top_surface_atom
        xy = (surface[n].x, surface[n].y)
    elif args.direct:
        xy = get_cartesian_xy(surface, args.direct)
    else:
        print 'Either the top atom or direct coordinate of the adsorption \
should provided'

    orientations = 'x -x y -y z -z'.split()

    for o in orientations:
        surf = surface.copy()
        mol = molecule.copy()

