#! /usr/bin/env python

import argparse
import ase.io
import numpy as np

parser = argparse.ArgumentParser(
        description="""Takes input POSCAR1 and POSCAR2 as input.
                Returns norm(POSCAR2 - POSCAR1) as output.""")

parser.add_argument('-p1', '--poscar1')
parser.add_argument('-p2', '--poscar2')
parser.add_argument('-o', '--modecar', default='MODECAR')

args = parser.parse_args()

atoms1 = ase.io.read(args.poscar1, format='vasp')
atoms2 = ase.io.read(args.poscar2, format='vasp')

p1 = atoms1.get_scaled_positions()
p2 = atoms2.get_scaled_positions()

atoms3 = atoms1.copy()
p = p2 - p1

# pbc
m, n = p.shape
for i in range(m):
    for j in range(n):
        if p[i,j] > 0.5:
            p[i,j] -= 1.0
        elif p[i,j] < -0.5:
            p[i,j] += 1.0

atoms3.set_scaled_positions(p)
nrm = np.linalg.norm(atoms3.get_positions())
modecar = atoms3.get_positions()/nrm

np.savetxt(args.modecar, modecar, fmt='%20.10E')
