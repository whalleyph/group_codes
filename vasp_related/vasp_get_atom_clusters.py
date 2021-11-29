#! /usr/bin/env python

import argparse
import cluster
import ase.io
from vasp_get_nearest_image_distances import get_nearest_image_distances

def get_clusters(atoms, rmax=1.4):
    nearest_image_disgances = get_nearest_image_distances(atoms)
    rij = lambda i,j: nearest_image_disgances[i,j]
    indices = range(len(atoms))
    cl = cluster.HierarchicalClustering(indices, rij)
    clusters = cl.getlevel(rmax)
    return clusters

def is_cluster_in_vacuum(atoms, cluster, surface):
    smin, smax = surface
    pos = atoms.get_scaled_positions()
    is_atom_in_vacuum = lambda i: smax < pos[i][2] < 1.0 or 0.0 < pos[i][2] < smin
    # returns true only when all atoms are in vacuum
    result = all([is_atom_in_vacuum(iatom) for iatom in cluster])

    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='POSCAR')
    parser.add_argument('-r', '--rmax', type=float, default=1.6)
    parser.add_argument('-s', '--surface', type=float, nargs=2, default=[0.222, 0.777])
    
    args = parser.parse_args()

    atoms = ase.io.read(args.infile, format='vasp')

    clusters = [i for i in get_clusters(atoms, args.rmax) if  len(i) > 1]

    for cluster in clusters:
        print atoms[cluster].get_chemical_symbols()
        print is_cluster_in_vacuum(atoms, cluster, surface=[0.222, 0.777])
        print
