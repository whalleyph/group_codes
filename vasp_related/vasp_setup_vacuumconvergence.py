from ase.lattice.surface import *
from ase.io import vasp
from ase.calculators.vasp import Vasp
import os
import numpy as np

pwd = os.getcwd()

a = 5.4524953358405277/2.
c = 8.6008007405135629/2.

vacuum = np.arange(5, 27.5, 2.5)
hvacuum = 0.5*vacuum
for (i,j) in zip(vacuum,hvacuum):
    ru_0001 = hcp0001('Ru', size=(1,1,6), a=a, c=c, vacuum=j)
    dir = "%3.1f" % i
    os.mkdir(dir)
    calc = Vasp(xc = 'PBE', 
                kpts = (12, 12, 6),
                prec = 'Accurate',
                ibrion = 2,
                ismear = 1,
                sigma = 0.1,
                lreal = 'Auto',
                encut = 400,
                nsw = 1)
    os.chdir(dir)
    calc.initialize(ru_0001)
    vasp.write_vasp('POSCAR', calc.atoms_sorted, symbol_count = calc.symbol_count, direct=True)
    calc.write_incar(ru_0001)
    calc.write_kpoints()
    calc.write_potcar()
    open('ready', 'w')
    os.chdir('..')
