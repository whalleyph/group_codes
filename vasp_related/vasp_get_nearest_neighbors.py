#! /usr/bin/env python

import ase, ase.io
import numpy as np
import argparse

def get_images(position, n):
    images = []
    for i in range(-n, n+1):
        for j in range(-n, n+1):
            for k in range(-n, n+1):
                images.append([i+position[0], j+position[1], k+position[2]])
    return np.array(images)

def get_nearest_neighbors(iatom, atoms, rcut=5.0):
    nearest_neighbors = []
    cell = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    natoms = len(atoms)
    for i in range(natoms):
        images = get_images(positions[i], 1)
        for j in images:
            rij = j - positions[iatom]
            cartesian_coordinates = np.dot(cell.transpose(), rij)
            rij = np.linalg.norm(cartesian_coordinates)
            if rij < rcut:
                nearest_neighbors.append([atoms[iatom].symbol, atoms[i].symbol, rij, i])

    return nearest_neighbors
        


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='script to put a molecule in the center of a cubic box')
    
    parser.add_argument('-i', '--input_file', default='POSCAR',
            help='input file in POSCAR format')
    parser.add_argument('-n', '--atoms_list', nargs='*', type=int,
            help='list of atoms for whom the nearest neighbors should be computed (count starts at 1)')
    parser.add_argument('-rcut', '--rcut', default=4.0, type=float,
            help='cutoff radius within which the nearest neighbors should be calculated')

    args = parser.parse_args()

    atoms = ase.io.read(args.input_file)
    if args.atoms_list:
        atoms_list = np.array(args.atoms_list) - 1
    else:
        atoms_list = range(len(atoms))

    for i in atoms_list:
        print '====== %i ======' % (i)
        nn = get_nearest_neighbors(i, atoms, rcut=args.rcut)
        nn = sorted(nn, key=lambda atom: atom[2])
        for i in nn:
            print i
