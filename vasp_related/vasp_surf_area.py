#!/usr/bin/env python2

#written by BK.
import ase.io
import numpy as np

def surf_area(slab1):
        return np.linalg.norm(np.cross(slab1.cell[0],slab1.cell[1]))


atoms=ase.io.read("POSCAR",format="vasp")


Area = surf_area(atoms)#np.linalg.norm(np.cross(v1,v2))

print "Surface area: %.2f A^2 = %.2f nm^2"%(Area,Area/100)
