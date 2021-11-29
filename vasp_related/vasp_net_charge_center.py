#!/usr/bin/env python

# compute local potential with dipole calculation on
#from ase.lattice.surface import fcc111, add_adsorbate
#from vasp import Vasp
import numpy as np
from ase.calculators.vasp import Vasp
from ase.calculators.vasp import VaspChargeDensity
import argparse
import os

parser = argparse.ArgumentParser(
            description='')

parser.add_argument('-f', '--file', default='./CHG')
args = parser.parse_args()

path="/".join(args.file.split("/")[0:-1])
print "Working dir: ",path

vcd = VaspChargeDensity(args.file)

cd = np.array(vcd.chg[0])

n0, n1, n2 = cd.shape

s0 = 1.0 / n0
s1 = 1.0 / n1
s2 = 1.0 / n2

X, Y, Z = np.mgrid[0.0:1.0:s0,
                   0.0:1.0:s1,
                   0.0:1.0:s2]

C = np.column_stack([X.ravel(),
                     Y.ravel(),
                     Z.ravel()])

os.chdir(path)
#print os.getcwd()
calc = Vasp(restart=True)
#atoms=ase.io.read("CONTCAR",format="vasp")
atoms = calc.get_atoms()
uc = atoms.get_cell()
real = np.dot(C, uc)

# now convert arrays back to unitcell shape
x = np.reshape(real[:, 0], (n0, n1, n2))
y = np.reshape(real[:, 1], (n0, n1, n2))
z = np.reshape(real[:, 2], (n0, n1, n2))


nelements = n0 * n1 * n2
voxel_volume = atoms.get_volume() / nelements
total_electron_charge = cd.sum() * voxel_volume

electron_density_center = np.array([(cd * x).sum(),
                                    (cd * y).sum(),
                                    (cd * z).sum()])
electron_density_center *= voxel_volume
electron_density_center /= total_electron_charge

electron_dipole_moment = -electron_density_center * total_electron_charge
"""
vcd = VaspChargeDensity('./CHG')
#x, y, z, cd = calc.get_charge_density()
n0, n1, n2 = cd.shape
nelements = n0 * n1 * n2
voxel_volume = lab.get_volume() / nelements
total_electron_charge = cd.sum() * voxel_volume

electron_density_center = np.array([(cd * x).sum(),
                                    (cd * y).sum(),
                                    (cd * z).sum()])
electron_density_center *= voxel_volume
electron_density_center /= total_electron_charge
"""

print 'electron-density center = {0}'.format(electron_density_center)
uc = atoms.get_cell()

# get scaled electron charge density center
sedc = np.dot(np.linalg.inv(uc.T), electron_density_center.T).T

# we only write 4 decimal places out to the INCAR file, so we round here.
sedc = np.round(sedc, 5)

#print sedc
print "Scaled electron charge density center: "
print "DIPOL= %.5f %.5f %.5f"%(sedc[0],sedc[1],sedc[2])
