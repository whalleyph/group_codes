#! /usr/bin/env python

import os
import ase
from ase.io import vasp
from ase.calculators.vasp import Vasp
import numpy as np


# Eq. bond length
deq = 1.2856

# Range of bond lengths from -40 % deq to +40 % deq
dmin = 0.6*deq
dmax = 1.4*deq

# delta d of 5% 
deld = 0.05*deq

drange = np.arange(dmin, dmax, deld)
print drange

calc = Vasp(
        xc='PBE',
        encut=400,
        ibrion=2,
        ismear=-1,
        sigma=0.1,
        prec='Accurate',
        isif=2,
        nsw=40,
        ediff=1e-5,
        lreal='Auto',
        ispin=2,
        kpts=(1,1,1), gamma=True
        )

for d in drange:
    dir = "%6.5f" % (d)
    print dir
    ozone = ase.Atoms('O3', positions=[(0.0, 0.0, 0.0),
                                       (0.0, 0.0, d),
                                       (1.133188, 0.0, -0.607167)],
            cell=[(8.0, 0.0, 0.0), (0.0, 8.0, 0.0), (0.0, 0.0, 8.0)],
            pbc=(1, 1, 1)
            )
    constraint = ase.constraints.FixAtoms([0,1])
    ozone.set_constraint(constraint)
    calc.initialize(ozone)
    if not os.path.exists(dir):
        os.mkdir(dir)
        os.chdir(dir)
        vasp.write_vasp('POSCAR', calc.atoms_sorted, symbol_count=calc.symbol_count)
        calc.write_incar(ozone)
        calc.write_kpoints()
        calc.write_potcar()
        os.chdir('..')
    else:
        print '%s exists... Skipping' % (dir)
