#!/usr/bin/env python
# -*- coding=utf-8 -*-

from sys import argv,exit

import numpy as np
from numpy import array as npa

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar, BSVasprun #, Potcar, Poscar
#from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital

from copy import deepcopy as dc
#import inspect

import argparse


if __name__ == "__main__":

    #The code is adapted from http://gvallver.perso.univ-pau.fr/?p=587.


    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", type=str, default="./vasprun.xml",
                    help="XML file of the (partial/projected) DOS calculation.")

    parser.add_argument("-at","--atomTypes", type=str, default="all",
                    help="Atom types to include  in the plotted band diagram.Default: all, a proper use: -at C,H,O,Pt")
    parser.add_argument("-r","--range", type=str, default="-3:3",
                    help="""Energy range to plot the DOS and band structures. Def: -3 to 3 eV. \nproper use: -r [-4:4] (from -4 to 4)\n-r all (the whole energy range)""")
    parser.add_argument("-ymax","--ymax", type=float,
                    help="Set the max density (yaxis) to show in the DOS plot. Def: maximum density within the given energy range.")
    parser.add_argument("-d","--decomp", action='store_true', default=False,                    help="Turn on the plotting of decomposed (partial) DOS instead of total DOS., Def: False")
    parser.add_argument("-noEf","--noEfcorr", action='store_true', default=False,                    help="Efermi is not set to zero (default: Ef lies at 0 eV)")
    parser.add_argument("-s","--shift", type=float, default=0.0,
                    help="Shifts the energies by this value (in eV). Def: No shift.")
    parser.add_argument("-hide_d","--hide_d", action='store_true', default=False,                    help="Hide d orbital contributions, as they generally overshadow the other orbitals. (Def: show d orbitals)")
    args = parser.parse_args()

    inclAtoms=[]
    if args.atomTypes=="all": 
#inclAtoms=[x for x in aTypes if x not in inclAtoms]
        for x in aTypes: 
            if x not in inclAtoms: inclAtoms.append(x)
    else: 
        inclAtoms=args.atomTypes.split(",")

    # read data
    # ---------


    # density of state
    #dosrun = Vasprun("./dos/vasprun.xml")
    #f="./vasprun.xml"
    f=args.input
    print "Reading DOS info from %s"%f
    dosrun = Vasprun(f)
    #dosrun = BSVasprun(f) #not comptaible with DOS data !!
    
    # Atomic Information (added by BK).
    aTypes=dosrun.atomic_symbols
    
    #data = Procar("./PROCAR").data

    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    #gs = GridSpec(1, 2, width_ratios=[2, 1])
    gs = GridSpec(1, 1)
    fig = plt.figure(figsize=(11.69, 8.27))
    
    #tt=""
    #for at in inclAtoms: 
    #    if at in aTypes: tt += "%s, " % at
    #tt="Band diagram for "+tt[0:-2]+" atoms"
    #fig.suptitle("Bands diagram of graphene")
    #fig.suptitle(tt)

    ax1 = plt.subplot(gs[0])
    #ax1=plt
    #ax2 = plt.subplot(gs[1])  # , sharey=ax1)

 
    # set ylim for the plot
    # ---------------------

    if args.noEfcorr:
        dosrun.efermi=0.0
        #bands.efermi=0.0

    if args.shift:
        dosrun.efermi-=args.shift
        #bands.efermi-=args.shift

    efermi=dosrun.efermi
    energies=dosrun.tdos.energies - efermi
    #print len(energies)

    if args.range=="all":
        #emin -= bands.efermi + 1
        #emax -= bands.efermi - 1
        emin=min(energies)#-1.0
        emax=max(energies)#+1.0
    else:
        w=args.range.replace('[','').replace("]","").split(":")
        try:
            emin=int(w[0])
            emax=int(w[1])
        except: 
            print "Error in energy range input. Please check your input. It should read -r -3:3 etc. Proceeding with the whole energy range." 
            #emin -= bands.efermi + 1
            #emax -= bands.efermi - 1
            emin=min(energies)#-1.0
            emax=max(energies)#+1.0


    ax1.set_xlim(emin, emax)

   
    # Process  Density of state Info
    # ----------------
    #nAtoms=float(len(dosrun.pdos))
    
    nedos=len(dosrun.pdos[0][Orbital.s][Spin.up]) #Grid size used for DOS calcultion.

    #orbs=Orbital.__members__
    #orbs=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']
    orbs=['s', 'py', 'pz', 'px', 'dxy','dyz', 'dz2', 'dxz','dx2']
  

    #Initialize data dictionary.
    data={}#Holds the total PDOS/atom data for each atom type (summed for all atoms of that type) and each orb type. Format data[atomType][orbType][PDOS at given Energy]
    noAtTypes={}
    actOrbs={} #actual orbs
    for at in aTypes:
        noAtTypes[at]=0 #holds the number of atoms of this specific type.
        if at not in inclAtoms: continue
        actOrbs[at]=[]
        actOrbs[at]=dc(orbs)

        #orbs=dosrun.pdos[i].keys()
        #print orbs
        tmp={}
        for orb in orbs:
            tmp[orb]=np.zeros(nedos)
            #tmp.append([orb,])
        data[at]=tmp


    #print data

    #Fill in the data.
    for i in range(len(dosrun.pdos)): #iterate over all atoms and get the DOS/atom.
        at=aTypes[i]
        noAtTypes[at]+=1
        if at not in inclAtoms: continue
        for orb in orbs: #loop over orbitals.
            x=Orbital.__members__["%s"%orb]
            if dosrun.pdos[i].has_key(x):
                #print "bk:",x
                data[at][orb] += npa(dosrun.pdos[i][x][Spin.up]) / noAtTypes[at]
                #if orb not in actOrbs[at]: actOrbs[at].append(orb)
            else:actOrbs[at].remove(orb)


    tDos={}#Total DOS for each atom type, summed over all available orbitals.
    for at in actOrbs.keys():
        if at not in inclAtoms: continue
        tDos[at]= np.zeros(nedos)
        #for key in actOrbs[at]:
        #print orbs
        for key in orbs:
            #print actOrbs
            #print "my key:", key
            if key not in data[at].keys(): actOrbs[at].remove(key);print "not in key",key ;continue
            #print key
            #key=Orbital.__members__["%s"%key]
            #print data[at][key]
            if np.count_nonzero(data[at][key])==0: #if all zeros then delete the entry
                print "Skipping %s - %s: zero contribution."%(at,key)
                actOrbs[at].remove(key)
                #data[at].pop(key)
                del data[at][key]
                #print data[at].keys()
            else:
                tDos[at] += data[at][key] #/  noAtTypes[at]  #NEED TO NORM??


    #DOS plot settings.
    #max_dens=max(dosrun.tdos.densities[Spin.up]/len(aTypes))
    
   # style
    ax1.set_ylabel("DOS / atom")
    ax1.set_xlabel(r"$E - E_f$   /   eV")
    ax1.grid()

    
    if not args.decomp: #Whether to plot the total DOS of each atom  instead of Partial DOS. 
        #uses tDos instead of data
        total=np.zeros(nedos)
        #print actOrbs.keys()
        for at in actOrbs.keys():
            total += tDos[at]
            ax1.plot( energies, tDos[at],
              #label="%s - %s"%(at,"total"), linewidth=2, color=(0.6, 0.6, 0.6))
                      label="%s"%(at), linewidth=2)

        
        ax1.plot( energies, total,
                  label="%s"%("total"), linewidth=2, color=(0.6, 0.6, 0.6))
 
        ax1.fill_between(energies, total,
                         0,
                         color=(0.7, 0.7, 0.7),
                         facecolor=(0.7, 0.7, 0.7))
            #max_dens=max(tDos[at])

       #For correct normalization of DOS in a given subregion (e.g. -3 ev to 3 eV).
        if 1:
            max_dens=0.0
            for i in range(len(energies)): #nedos
                if emin<= energies[i] <= emax:
                    if total[i] > max_dens: max_dens=total[i]
        else:
            ##max_dens=max(tDos[at])
            max_dens=max(total)#/len(actOrbs.keys())
            #max_dens=max(dosrun.tdos.densities[Spin.up]/len(aTypes))


        ax1.set_ylim(0, max_dens)

# fermi level at 0
        ax1.vlines(x=0, ymin=0, ymax=max_dens, color="k", lw=2,linestyle="dashed")#!!!

        ax1.legend(fancybox=True, shadow=True, prop={'size': 14})
        plt.subplots_adjust(wspace=0)

        plt.show()

    #plt.savefig(argv[0].strip(".py") + ".pdf", format="pdf")
        plt.savefig("dos.png", format="png")


        exit()


    total=np.zeros(nedos)
    ats=actOrbs.keys()
    ats.sort()

    #ats=inclAtoms

    for at in ats:
        if at not in inclAtoms: continue     
        if at=='H' or at=='He': actOrbs[at]=['s']

        for orb in actOrbs[at]:           
            if at=="C" and orb=="s":
                ax1.plot( energies, data[at][orb],
                          label="%s - %s"%(at,orb), linewidth=2,c="r")
                total += data[at][orb]
            elif at=="C" and orb=="py":
                total_pxy=(data[at]["px"]+data[at]["py"])#/2
                ax1.plot( energies, total_pxy ,
                          label="%s - %s"%(at,"px+py"), linewidth=2,c='g')
                total += total_pxy

            elif at=="C" and orb=="pz":
                ax1.plot( energies, data[at][orb],
                          label="%s - %s"%(at,orb), linewidth=2,c="b")
                total += data[at][orb]

            elif (at=="O" or at=="Pt") and orb=="pz":
                ax1.plot( energies, data[at][orb],
                          label="%s - %s"%(at,orb), linewidth=2,c="k")
                total += data[at][orb]

            elif at=="H" and orb=="s":
                ax1.plot( energies, data[at][orb],
                          label="%s - %s"%(at,orb), linewidth=2,c="m")
                total += data[at][orb]


            elif orb=="s":
                ax1.plot( energies, data[at]["s"],
                #ax1.plot( energies, (data[at]["px"]+data[at]["py"]),
                          label="%s - %s"%(at,"s"), linewidth=2)
                total += data[at]["s"]
            elif orb=="py":
                total_pxy=(data[at]["px"]+data[at]["py"])#/2
                ax1.plot( energies, total_pxy ,
                          label="%s - %s"%(at,"px+py"), linewidth=2)
                total += total_pxy

            elif   orb=='dxy':
                if args.hide_d : continue
                total_d=(data[at]["dxy"]+data[at]["dxz"]+data[at]["dyz"]+data[at]["dx2"]+data[at]["dz2"])#/5
                ax1.plot( energies,  total_d,
                          label="%s - %s"%(at,"all d"), linewidth=2)
                total += total_d

            elif orb=="px" or  orb=='dyz' or orb=='dz2' or orb== 'dxz' or orb== 'dx2':continue
            
            else: 
                ax1.plot( energies, data[at][orb],
                          label="%s - %s"%(at,orb), linewidth=2)
                total += data[at][orb]

    ax1.plot( energies, total,
              #label="%s - %s"%(at,"total"), linewidth=2, color=(0.6, 0.6, 0.6))
              label="%s"%("total"), linewidth=2, color=(0.6, 0.6, 0.6))

        
    ax1.fill_between(energies, 
                     0,total,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.7, 0.7, 0.7))
    
    #For correct normalization of DOS in a given subregion (e.g. -3 ev to 3 eV).
    if 1:
        max_dens=0.0
        for i in range(len(energies)): #nedos
            if emin<= energies[i] <= emax:
                if total[i] > max_dens: max_dens=total[i]
    else:
        ##max_dens=max(tDos[at])
        max_dens=max(total)#/len(actOrbs.keys())
        

    if args.ymax:
        try: max_dens=args.ymax
        except: 
            print "ymax value is not valid. Max density value will be used instead."
            

    ax1.set_ylim(0, max_dens)
# fermi level at 0
    ax1.vlines(x=0, ymin=0, ymax=max_dens, color="k", lw=2,linestyle="dashed")


    ax1.legend(fancybox=True, shadow=True, prop={'size': 14})
    plt.subplots_adjust(wspace=0)

    plt.show()

    #plt.savefig(argv[0].strip(".py") + ".pdf", format="pdf")
    plt.savefig("dos.png", format="png")

    #plt.clf()
    #plt.close()
    #plt.show()

    exit()
