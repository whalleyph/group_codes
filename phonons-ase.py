#!/usr/bin/env python

from ase.spacegroup import crystal
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo

from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.thermochemistry import HarmonicThermo


from ase import Atoms, Atom
from ase.calculators.vasp import Vasp,Vasp2
from sys import exit
from os import popen,listdir
from numpy import array

#from ase.build import molecule
#from ase.calculators.emt import EMT
#from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

import inspect

import os.path
import ase.io

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--T", type=float, nargs="*", default=[298.15],
                    help="Temperature, def: 298.15 K")
parser.add_argument("-p","--P", type=float, nargs="*",default=[101325.],
                    help="Pressure, def: 101325.0 Pa (1 atm)")

parser.add_argument("-v","--verb", action='store_true', default=False,
                    help="verbose output")

args = parser.parse_args()


x = [i for i in listdir('.') if i [0:3] == "TMP"]
if len(x)>0: cdir="TMP%d"%(len(x)+1) #x.split("TMP")
else: cdir="TMP1"

#try: 
#	atoms = ase.io.read('vasprun.xml',index=0)
#	print "vasprun.xml was read successfully."
#except:
try:
	calc = Vasp(restart=True)
	atoms = calc.get_atoms()
	print "VASP run was read succesfully from OUTCAR."
#else:
except:
	calc=Vasp2()
	calc.read_incar()
	atoms=ase.io.read("POSCAR",format="vasp")
        #try:atoms=ase.io.read("POSCAR",type="vasp",format="vasp4)
	calc.set(xc="pbe",nsw=0,ibrion=-1)
	calc.directory="."#cdir
	atoms.set_calculator(calc)
	print "VASP run could not be read, starting a new run..."

try: Ep = atoms.get_potential_energy()
#except:Ep=float(popen("grep 'energy  without entropy=' OUTCAR | tail -n 1").readlines()[0][0:-1].split()[-1]); #print Ep
except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);print Ep

print atoms
print Ep

try:
        freqs= calc.read_vib_freq()
        #print freqs
        freqs=freqs[0] #Take only the real frequencies.
except: freqs=[]


potentialenergy = Ep
calc.directory=cdir
atoms.set_calculator(calc)

# Phonon analysis
N = 1 #supercell size
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
#ph.clean(empty_files=True) #does not support empty_files option
ph.run()
ph.read(acoustic=True)
phonon_energies, phonon_DOS = ph.dos(kpts=(1, 1, 1), npts=3000,
                                     delta=5e-4)

# Calculate the Helmholtz free energy
thermo = CrystalThermo(phonon_energies=phonon_energies,
                       phonon_DOS=phonon_DOS,
                       potentialenergy=potentialenergy,
                       formula_units=1)

print "Phonon vibrations"
for T in args.T:
        F = thermo.get_helmholtz_energy(temperature=T,verbose=args.verb)
        print "T= %6.1f K F= %.4f eV"%(T,F)

#F = thermo.get_helmholtz_energy(temperature=298.15)


#using vibrations
vib = Vibrations(atoms)
#vib.clean(empty_files=True)
vib.run()
vib_energies = vib.get_energies()
#vib_freqs = vib.get_frequencies(method='standard', direction='central')

print "%d frequencies were read..."%len(vib_energies)
print "No of Imaginary freqs: ", len([i*8065.54429 for i in vib_energies if i.imag >0])
#print "Summary: ", vib.summary()

#needed for consistency for comparing DynMat and VASP results.
#vib_energies.sort(reverse=True)
vib_energies.sort() #3N-6 modes are used, so making sure the lowest ones are deleted.
if 1:
        print "Deleting modes with low frequency (v<209 cm^-1 ~ 300K)..."
        vib_new=[]
        for i in vib_energies:
                #if i >= (209/8065.54429) :vib_new.append(i) #if lower than 209 cm-1, not include mode in partition function.
                if i > (10./8065.54429) and i.imag==0 :vib_new.append(i) #if lower than 209 cm-1, not include mode in partition function.
	print "%d modes deleted"%(len(vib_energies)-len(vib_new))
        vib_energies=vib_new

print "\nIdeal Gas Thermodynamics"
thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='nonlinear',
                        symmetrynumber=1, spin=0)

for T in args.T:
        for P in args.P:
                G = thermo.get_gibbs_energy(temperature=T, pressure=P,verbose=args.verb)
                print "T= %6.1f K P= %.1f Pa G= %.4f eV"%(T,P,G)

#G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)

print "\nHarmonic Thermodynamics"
thermo = HarmonicThermo(vib_energies=vib_energies, potentialenergy=potentialenergy)#,
for T in args.T:
        F = thermo.get_helmholtz_energy(temperature=T,verbose=args.verb)
        print "T= %6.1f K F= %.4f eV"%(T,F)

#G = thermo.get_helmholtz_energy(temperature=298.15)

