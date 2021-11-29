import ase
import time
import numpy as np

# Seed the random number generator
np.random.seed(seed=int(time.time()))

# Bond length of O2
d = 1.21

# van der Walls radius of O
vdw = 1.52

# No. of O2 molecules in 3 dimensions
nx, ny, nz = [3,3,3]

# Translation distance between neighbouring molecules
tx, ty, tz = [1.5*vdw, 1.5*vdw, 1.5*vdw]

atoms = ase.Atoms()

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            O2 = ase.Atoms('O2', positions=[(0,0,0), (0, 0, d)])
            O2.translate([i*tx, j*ty, k*tz])
            # Rotate the molecule by a random angle
            phi, theta, psi = np.pi*np.random.random_sample(3)
            O2.rotate_euler(center='COM',phi=phi, theta=theta, psi=psi)
            atoms.extend(O2)

#print atoms.get_positions()
#ase.visualize.view(atoms)
atoms.pbc = [1, 1, ]
atoms.cell = [8, 8, 8]
ase.io.write('O2.xyz', atoms)
ase.io.write('POSCAR', atoms)

# Attach calculator
#atoms.set_calculator(ase.EMT())

## Molecular Dynamics
#dyn = ase.md.verlet.VelocityVerlet(atoms, dt=1.0*ase.units.fs)
#nblocks = 100
#for i in range(nblocks):
#    pot = atoms.get_potential_energy()
#    kin = atoms.get_kinetic_energy()
#    print '%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin)
#    dyn.run(steps=10)
#
#ase.io.write('O2_end.xyz', atoms)
