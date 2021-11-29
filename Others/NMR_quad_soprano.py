#!/bin/env python
# Basic imports
import os, sys
import matplotlib.pyplot as plt
from sys import exit
import numpy as np

# File location
filepath = os.path.split(os.path.abspath(__file__))[0]
sys.path.insert(0, os.path.join(filepath, '..')) # This to add the Soprano path to the PYTHONPATH so we can load it without installing it
datapath = os.path.join(filepath, '')#'tutorial_data')

# Other useful imports

import glob

import numpy as np

import ase
from ase import io as ase_io

from soprano.collection import AtomsCollection
from soprano.calculate import nmr#,simpson
from soprano.calculate.nmr import NMRFlags
from soprano.properties.nmr import *

# List all files in the tutorial directory

cifs = glob.glob(os.path.join(datapath, '*.magres'))
cifs=("config_2.magres","config_3.magres","config_8.magres","config_11.magres")
#print datapath,cifs
#aColl = AtomsCollection(cifs, progress=True) # "progress" means we will visualize a loading bar
print

lfreq=213.5 #MHz #neg pos farketmiyor
#lfreq=122.39
#lfreq=900
aType="Ga" #or #71Ga
ref=1759 #ours
#ref=1502 #middlemiss
#ref=1650 #better match with exp.
#ref=0
mn=1200;mx=1800
#mn=-200;mx=200 #this covers the unshifted (CS-iso corr) range.
#mn=200;mx=500
#mn=100;mx=400 #best for ref case
#mn+=ref;mx+=ref #to correct for the ref
print mn,mx
nbins=1000
intensitiesA=[];intensitiesB=[];intensitiesC=[];peaks=[]
for cif in cifs:
    atoms=ase.io.read(cif)#cifs[0])
    myNMR=nmr.NMRCalculator(atoms)
    myNMR.set_larmor_frequency(larmor_frequency=lfreq,larmor_units='MHz',element=aType)

    myNMR.set_reference(ref, aType)#Needs to be before powder.

    myNMR.set_powder(N=32,mode='hemisphere')#to give multiple orientation to the crystal. Higher N (orientations) the more accurate  (and costly) it is. Not much diff btw sphere and hemisphere. hemisphere is the default.
    #myNMR.set_single_crystal(theta=180, phi=180) #does not work well.


    #print myNMR.get_larmor_frequency(aType)
    print myNMR._references
    print "Larmor frequency for {0:.3s}: {1:.2f} MHz".format(aType,myNMR.get_larmor_frequency(aType))

    print "Field: {0:.2f} T".format(myNMR.B)
    print "Powder orientations: {0}".format(len(myNMR._orients[0]))

    labels = atoms.get_array('labels')
    indices = atoms.get_array('indices')
    # One can also recreate the Jmol-style labels if needed
    jmol_labels = ["{0}_{1}".format(l, i) for l, i in zip(labels, indices)]
    #print "The labels are {0}".format(', '.join(jmol_labels))

    # Isotropy, Anisotropy and Asymmetry (Haeberlen convention)
    iso = MSIsotropy.get(atoms)
    aniso = MSAnisotropy.get(atoms)
    asymm = MSAsymmetry.get(atoms)

    vzz = EFGVzz.get(atoms)
    #qP = EFGQuadrupolarConstant(isotopes={'H': 2}) # Deuterated; for the others use the def.
    qP = EFGQuadrupolarConstant()
    qC = qP(atoms)/1e6 # To MHz
    #EFGaniso=EFGAnisotropy.get(atoms)
    EFGasymm=EFGAsymmetry.get(atoms)

    print 'Label\tIsotropy\tAnisotropy\tAsymmetry\tVzz\t\tCQ/Chi\t\teta_Q'
    for i, jl in enumerate(jmol_labels):
        if aType in jl:  print '{0}\t{1:.2f} ppm\t{2:.2f} ppm\t{3:.2f}\t\t{4:5.2f} au\t{5:5.2f} MHz\t{6:5.2f}'.format(jl, ref-iso[i], aniso[i], asymm[i],vzz[i], qC[i],EFGasymm[i])

 
    intensities1,freqs1=myNMR.spectrum_1d(aType, min_freq=mn, max_freq=mx, bins=nbins,freq_broad=None, freq_units='ppm', effects=NMRFlags.STATIC )
 
    intensities2,freqs2=myNMR.spectrum_1d(aType, min_freq=mn, max_freq=mx, bins=nbins,freq_broad=1, freq_units='ppm', effects=45)#NMRFlags.Q_MAS )#contains CS_ISO+Q_2_SHIFT+Q_2_ORIENT_MAS) MAS=41 by def, total=47 could be better incl. CS_ORIENT and Q_1_ORIENT as well. #45 matches well with other programs. #THIS STILL NEEDS TESTING

    intensities3,freqs3=myNMR.spectrum_1d(aType, min_freq=mn, max_freq=mx, bins=nbins,freq_broad=None, freq_units='ppm', effects=NMRFlags.CS )

    #"""
    if len(intensitiesA)==0:intensitiesA=intensities1 #needed for when using ref shift.
    else:intensitiesA+=intensities1
    if len(intensitiesB)==0:intensitiesB=intensities2
    else:intensitiesB+=intensities2
    if len(intensitiesC)==0:intensitiesC=intensities3
    else:intensitiesC+=intensities3
    """
    if len(intensitiesA)==0:intensitiesA=np.flip(intensities1,axis=0) #needed for when using ref shift.
    else:intensitiesA+=np.flip(intensities1,axis=0)
    if len(intensitiesB)==0:intensitiesB=np.flip(intensities2,axis=0)
    else:intensitiesB+=np.flip(intensities2,axis=0)
    if len(intensitiesC)==0:intensitiesC=np.flip(intensities3,axis=0)
    else:intensitiesC+=np.flip(intensities3,axis=0)
    """
    #peaks.append(peaks3)

peaks=np.array(peaks)
if 0:
    intensities1/=max(intensities1)
    intensities2/=max(intensities2)
    intensities3/=max(intensities3)
#plt.plot(ref-freqs1,intensitiesA,color="green",label="Static")
plt.plot(ref-freqs2,intensitiesB,color="blue",label="MAS")
#plt.plot(ref-freqs3,intensitiesC,color="red",label="CS")
#plt.vlines(ref-peaks,0,max(max(intensitiesA),max(intensitiesB),max(intensitiesC)),color="black",label="delta_iso")
plt.gca().invert_xaxis() #plt.gca() is Axes object

plt.legend()

plt.show()


exit()

"""
    A list of the currently available NMRFlags to be used in conjunction with
    methods that require a list of effects of interest: ::

        NMRFlags.CS_ISO     => chemical shielding, isotropic effect
                .CS_ORIENT  => chemical shielding, orientation dependent
                               effects
                .CS         => chemical shielding, everything
                .Q_1_ORIENT => quadrupolar, 1st order, orientation dependent
                               effects
                .Q_2_SHIFT  => quadrupolar, 2nd order, isotropic shift
                .Q_2_ORIENT_STATIC => quadrupolar, 2nd order, orientation
                                      dependent effects; static limit
                .Q_2_ORIENT_MAS => quadrupolar, 2nd order, orientation
                                      dependent effects; ultrafast MAS limit
                .Q_2_STATIC => quadrupolar, 2nd order, all static effects
                .Q_2_MAS    => quadrupolar, 2nd order, all ultrafast MAS
                               effects
                .Q_STATIC   => quadrupolar, all static effects
                .Q_MAS      => quadrupolar, all ultrafast MAS effects
                .STATIC     => all static effects
                .MAS        => all ultrafast MAS effects

NMRFlags = NMRFlags(
    CS_ISO=1,
    CS_ORIENT=2,
    CS=1+2,
    Q_1_ORIENT=4,  # First order quadrupolar anisotropic effects
    Q_2_SHIFT=8,
    Q_2_ORIENT_STATIC=16,
    Q_2_ORIENT_MAS=32,
    Q_2_STATIC=8+16,
    Q_2_MAS=8+32,
    Q_STATIC=4+8+16,
    Q_MAS=8+32,
    STATIC=1+2+4+8+16,
    MAS=1+8+32)
"""

    #myNMR.spectrum_1d(element, min_freq=-50, max_freq=50, bins=100,freq_broad=None, freq_units='ppm', effects=NMRFlags.CS_ISO) NMRFlags.Q_2_ORIENT_MAS 




    #intensities3,freqs3,peaks3=myNMR.spectrum_1d(aType, min_freq=mn, max_freq=mx, bins=nbins,freq_broad=None, freq_units='ppm', effects=NMRFlags.CS )

plt.plot(ref-freqs1,intensitiesA,color="green",label="Static")
plt.plot(ref-freqs2,intensitiesB,color="blue",label="MAS")
plt.plot(ref-freqs3,intensitiesC,color="red",label="CS")

#intensities,freqs,peaks=myNMR.spectrum_1d('Ga', min_freq=mn, max_freq=mx, bins=1000,freq_broad=1, freq_units='ppm', effects=NMRFlags.CS_ISO )
#freqs=ref-freqs
#mx=max(intensities)
#intensities /=mx
#plt.plot(ref-freqs,intensities,color="black",label="CS_iso")
