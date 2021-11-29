#! /usr/bin/env python

import sys
sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
import potcar as pot
import poscar as pos
import numpy as np

potcar = pot.potcar().read('POTCAR')
nelectrons = np.array(potcar.get_nelectrons())
iontypes = potcar.get_iontypes()
poscar = pos.poscar().read('POSCAR')
nions = poscar.get_nions()
print "Ion types"
print iontypes
print "Valence electrons of each ion type"
print nelectrons
print "No. of each ion type"
print nions
print "Total no. of electrons"
totelect = int(np.sum(nelectrons*nions))
totions = np.sum(nions)
print totelect
print
print "Recommendations for NBANDS"
print
print "Non-spin polarized"
print "0.5*NELECT + 0.5*NIONS =", int(0.5*totelect+0.5*totions)+1
print
print "Spin polarized"
print "0.6*NELECT + NIONS =", int(0.6*totelect + totions) + 1
print "0.6*NELECT + 2*NIONS =", int(0.6*totelect + 2*totions) + 1
print 
