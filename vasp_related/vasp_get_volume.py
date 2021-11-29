#! /usr/bin/env python

import numpy as np
import ase, ase.io

cell = ase.io.read('POSCAR').get_cell()

a, b, c = cell[:]
b_cross_c = np.cross(b, c)
volume = np.dot(a, b_cross_c)
surface_area = np.linalg.norm(np.cross(a, b))

print 'Volume: {volume:.3f} Ang^3'.format(volume=volume)
print 'Surface Area: {surface_area: .3f} Ang^2'.format(surface_area=surface_area)
