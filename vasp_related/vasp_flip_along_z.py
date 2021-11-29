#!/usr/bin/env python2

#written by BK.
import ase.io
import numpy as np
from math import pi,cos,sin
from ase.visualize import view

def surf_area(slab1):
        return np.linalg.norm(np.cross(slab1.cell[0],slab1.cell[1]))


atoms=ase.io.read("POSCAR",format="vasp")
atoms.rotate('x',180*pi,rotate_cell=True)
view(atoms)
ase.io.write("flipped.vasp",atoms,format="vasp",vasp5=True,direct=True)
#Area = surf_area(atoms)#np.linalg.norm(np.cross(v1,v2))

#print "Surface area: %.2f A^2 = %.2f nm^2"%(Area,Area/100)
