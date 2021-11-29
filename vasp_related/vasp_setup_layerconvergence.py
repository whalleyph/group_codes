from ase.lattice.surface import *
from ase.io import vasp
from ase.calculators.vasp import Vasp
import os

a = 5.4524953358405277/2.
c = 8.6008007405135629/2.

nlayers = range(1, 13)
nbands = [8, 12, 16, 20, 24, 28, 34, 38, 44, 48, 52, 58]
for i, j in zip(nlayers,nbands):
    ru_0001 = hcp0001('Ru', size=(1,1,i), a=a, c=c, vacuum=7.5)
    # print ru_0001
    dir = "%2.2i" % i
    os.mkdir(dir)
    calc = Vasp(xc = 'PBE', 
                kpts = (24, 24, 1),
                gamma = True,
                prec = 'Accurate',
                ibrion = 2,
                ismear = 1,
                sigma = 0.1,
                lreal = 'Auto',
                encut = 400,
                isif = 2,
                nbands = j+4,
                nsw = 40)
    os.chdir(dir)
    calc.initialize(ru_0001)
    vasp.write_vasp('POSCAR', calc.atoms_sorted, symbol_count = calc.symbol_count, direct=True)
    calc.write_incar(ru_0001)
    calc.write_kpoints()
    calc.write_potcar()
    open('ready', 'w')
    os.chdir('..')
