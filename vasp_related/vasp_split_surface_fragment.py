#! /usr/bin/env python

import sys
import ase.io

if len(sys.argv) != 3:
    print 'Usage: %s <POSCAR/CONTCAR> <number of atoms in fragment>' %  (sys.argv[0])
    sys.exit()

f = sys.argv[1] # poscar or contcar
n = int(sys.argv[2]) # number of atoms in fragment

atoms = ase.io.read(f, format='vasp')

surface = atoms[:len(atoms)-n]
fragment = atoms[len(atoms)-n:]

print 'Assuming fragment at end of the configuration'
print surface
print fragment

ase.io.write('surface.vasp', surface, format='vasp')
ase.io.write('fragment.xyz', fragment, format='xyz')
