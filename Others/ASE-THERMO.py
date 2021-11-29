#!/usr/bin/env python

from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from sys import exit
from os import popen
from numpy import array
import ase.io

#from ase.build import molecule
#from ase.calculators.emt import EMT
#from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo,HarmonicThermo

#import inspect
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--T", type=float, nargs="*", default=[298.15],
                    help="Temperature, def: 298.15 K")
parser.add_argument("-p","--P", type=float, nargs="*",default=[101325.],
                    help="Pressure, def: 101325.0 Pa (1 atm)")

parser.add_argument("-v","--verb", action='store_true', default=False,
		    help="verbose output")

args = parser.parse_args()


if 0:
	calc = Vasp(restart=True)
	atoms = calc.get_atoms()
else:
	atoms=ase.io.read('vasprun.xml',index=0) #works better


try: Ep = atoms.get_potential_energy()
#except:Ep=float(popen("grep 'energy  without entropy=' OUTCAR | tail -n 1").readlines()[0][0:-1].split()[-1]); #print Ep
except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);print Ep
#except:exit() #Ep = atoms.get_potential_energy() 

try: 
	freqs= calc.read_vib_freq()
	#print freqs
	freqs=freqs[0] #Take only the real frequencies.
	print "%d imaginary frequency modes skipped"%len(freqs[1])
except: freqs=[]


if len(freqs)!=0: 
	vib_energies=[]
	for fr in freqs:
               #vib_energies.append( fr / 8065.54429)  # convert to eV from cm^-1
		vib_energies.append( fr / 1000.0 )  # convert to eV from meV
else:
#Check if DYMMAT calc
	print "No fequencies found in OUTCAR. Trying to read freq.dat (DynMat output)."
	try: 
		x=open("freq.dat",'r')#.readlines()
		#freqs=[]
		vib_energies=[]
		cnt=0
		for ln in x: #skip imaginary freqs.
			cnt +=1
			if ln.split()[-1]!="1":
				#freqs.append()
				vib_energies.append( float(ln.split()[0])/8065.54429) # convert to eV from cm^-1
		print "%d imaginary frequency modes skipped..."%(cnt-len(vib_energies))
	except:
		print "No useful information for vibrational frequencies."
		exit()

print "%d modes were read..."%len(vib_energies)
#needed for consistency for comparing DynMat and VASP results.
#vib_energies.sort(reverse=True)
vib_energies.sort() #3N-6 modes are used, so making sure the lowest ones are deleted. 
if 1: 
	thr=209  #if lower than XX cm-1, not include mode in partition function. 209 cm-1 corresponds to ~ 300K.
	thr=10
	print "Deleting modes with low frequency (v<%.1f cm^-1 )..."%thr
	vib_new=[]
	for i in vib_energies:
		if i >= (thr/8065.54429) :vib_new.append(i) 
		else: None; #print "Skipping mode %f cm^-1"%(i*8065.54429)
	print "%d modes deleted"%(len(vib_energies)-len(vib_new))
        vib_energies=vib_new
	
print "Potential/electronic energy (No ZPE, 0K): %.5f"%Ep
print "\nIdeal Gas Thermodynamics"
thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=Ep,
                        atoms=atoms,
                        geometry='nonlinear',
                        symmetrynumber=1, spin=0)

for T in args.T:
	for P in args.P:
		G = thermo.get_gibbs_energy(temperature=T, pressure=P,verbose=args.verb)
		print "T= %6.1f K P= %.1f Pa G= %.4f eV"%(T,P,G)
print "\nHarmonic Thermodynamics"
thermo = HarmonicThermo(vib_energies=vib_energies,potentialenergy=Ep)#,

for T in args.T:
	F = thermo.get_helmholtz_energy(temperature=T,verbose=args.verb)
	print "T= %6.1f K F= %.4f eV"%(T,F)


#1 atm = 101325. Pa
#1 mtorr= 0.13332237 Pas

exit()
