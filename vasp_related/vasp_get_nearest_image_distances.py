#! /usr/bin/env python

import sys
from itertools import product
from functools import partial
import ase, ase.io
import numpy as np
import argparse
from multiprocessing import Pool

def get_images(position, n):
    images = [(i+position[0], j+position[1], k+position[2]) for i,j,k in product(range(-n,n+1), repeat=3)]
    return images

def get_nearest_image_distance(cell, pos, images, pair):
    ri, rj_images = pos[pair[0]], images[pair[1]]
    rij_scaled = [rj-ri for rj in rj_images]
    rij_cartesian = [np.dot(cell.transpose(), rij) for rij in rij_scaled]
    distances_to_images = [np.linalg.norm(rij) for rij in rij_cartesian]
    nearest_image_distance = np.min(distances_to_images)
    return nearest_image_distance

def get_nearest_image_distances(atoms):
    natoms = len(atoms)
    cell = atoms.get_cell()
    pos = atoms.get_scaled_positions()
    atom_images = np.array([get_images(pos[i], 1) for i in range(natoms)])
    atom_pairs = [pair for pair in product(range(len(atoms)), repeat=2) if pair[0] < pair[1]]
    dist_matrix = np.zeros([natoms,natoms])
    p = Pool(8)
    nearest_image_distances = p.map(partial(get_nearest_image_distance, cell, pos, atom_images), atom_pairs)
    for pair, rij in zip(atom_pairs, nearest_image_distances):
        i, j = pair
        dist_matrix[i,j] = dist_matrix[j,i] = rij

    return dist_matrix


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='script to put a molecule in the center of a cubic box')
    
    parser.add_argument('-i', '--input_file', default='POSCAR',
            help='input file in POSCAR format')

    args = parser.parse_args()

    atoms = ase.io.read(args.input_file)

    matrix = get_nearest_image_distances(atoms)
    print matrix

    #for row in matrix:
    #    srow = ''
    #    for el in row:
    #        srow += ' %2.2f' % (el)
    #    print srow

