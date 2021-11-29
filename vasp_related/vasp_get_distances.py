#! /usr/bin/env python

import ase.io
import numpy as np
import argparse

def get_distance_matrix(atoms):
    positions = atoms.get_positions()

    n = len(positions)
    distances = np.zeros((n,n))

    for i in range(len(positions)-1):
        for j in range(i+1, len(positions)):
            a = positions[i]
            b = positions[j]

            dx = b[0] - a[0]
            dy = b[1] - a[1]
            dz = b[2] - a[2]

            dr = (dx*dx + dy*dy + dz*dz)**0.5

            distances[i,j] = distances[j,i] = dr

    return distances


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Print the distance matrix and the largest and the smallest distance')
    
    parser.add_argument('-m', '--molecule', default='molecule.xyz')
    
    args = parser.parse_args()

    atoms = ase.io.read(args.molecule)

    distances = get_distance_matrix(atoms)
    sorted_distances = sorted(distances.ravel())

    for i in sorted_distances:
        if i > 1E-8:
            print 'Minimum distance between two atoms: %.2f' % (i)
            break
    print 'Maximum distance between two atoms: %.2f' % (sorted_distances[-1])

