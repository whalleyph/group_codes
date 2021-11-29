#!/usr/bin/env python
# -*- coding=utf-8 -*-

from sys import argv,exit

import argparse
import numpy as np
from numpy import array as npa

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar, BSVasprun #, Potcar, Poscar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital

from copy import deepcopy as dc

def rgbline(ax, k, e, red, green, blue, alpha=1. ):#,lw=(1.5,1.5,1.5)):
    #alpha determines the opacity of the lines.!!
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=1.5)#, linewidth=lw)  # !!!
    
    ax.add_collection(lc)


if __name__ == "__main__":

    #The code is adapted from http://gvallver.perso.univ-pau.fr/?p=587.

    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--dosfile", type=str, default="../../vasprun.xml",
                    help="XML file of the partial DOS calculation.")
    parser.add_argument("-b","--bandfile", type=str, default="./vasprun.xml",
                    help="XML file of the partial band structure calculation.")
    parser.add_argument("-p","--procar", type=str, #default="./PROCAR",
                    help="PROCAR file from the band structure calculation. Not necessary to provide a PROCAR file, as the same info is in XML file.") 
    parser.add_argument("-onlypz","--onlypz", action='store_true', default=False,                    help="This is to plot only the pz-orbital contributions from all included atom types. Default: all orbital contributions are plotted.") 
    parser.add_argument("-rgb","--rgb", action='store_true', default=False,                    help="This is to plot the rgb color scale based on the contirbutions from s,px-py and pz ortbitals from the selected list of atom types.")
    parser.add_argument("-at","--atomTypes", type=str, default="all",
                    help="Atom types to include  in the plotted band diagram.Default: all, a proper use: -at C,H,O,Pt")

    parser.add_argument("-r","--range", type=str, default="-3:3",
                    help="""Energy range to plot the DOS and band structures. Def: -3 to 3 eV. \nproper use: -r [-4:4] (from -4 to 4)\n-r all (the whole energy range)""")
    parser.add_argument("-noEf","--noEfcorr", action='store_true', default=False,                    help="Efermi is not set to zero (default: Ef lies at 0 eV)")
    parser.add_argument("-s","--shift", type=float, default=0.0,
                    help="Shifts the energies by this value (in eV). Def: No shift.")

    args = parser.parse_args()
    onlypz=args.onlypz
    rgb=args.rgb

    # read data
    # ---------

    # kpoints labels
    try: 
        #ln=open("KPOINTS","r").readlines()[0][0:-1] #Assuming HighSymm Kpoints  path is given in the title line.
        ln=open("KPOINTS","r").readlines()
        

        path=ln[0][0:-1].replace("-","").split()
        sts=ln[1][0:-1].split("!")[0].split()
        steps=[int(st) for st in sts]

        #print path
        labels=[]
        for i in path:
            if i=="G": labels.append("$\\Gamma$")
            else:labels.append("$%s$"%i)
            #else:labels.append("u'$%s$'"%i)
        #print labels
        print "HighSymm Kpoints path retrieved from KPOINTS file: ", path
    except: 
        print "HighSymm Kpoints  path could not be read from KPOINTS file. Path must be given  in the title line of KPOINTS file (e.g. M G K M, etc.)"
        exit()

        # Original code.
    #path = HighSymmKpath(mg.Structure.from_file("./opt/CONTCAR")).kpath["path"]
        path = HighSymmKpath(mg.Structure.from_file("./vasprun.xml")).kpath["path"]
        print path
        labels = [r"$%s$" % lab for lab in path[0][0:4]]
        print labels

    ####################
    # density of state #
    ####################
    #dosrun = Vasprun("./dos/vasprun.xml")
    #f="../../vasprun.xml"
    f=args.dosfile
    print "Reading DOS info from %s."%f
    dosrun = Vasprun(f)
    #dosrun = BSVasprun(f) #not comptaible with DOS data !!
    
    # Atomic Information (added by BK).
    aTypes=dosrun.atomic_symbols
    
    inclAtoms=[]
    if args.atomTypes=="all": 
#inclAtoms=[x for x in aTypes if x not in inclAtoms]
        for x in aTypes: 
            if x not in inclAtoms: inclAtoms.append(x)
    else: 
        inclAtoms=args.atomTypes.split(",")
      
    #print inclAtoms

    #########
    # bands #
    #########
    #bands = BSVasprun("./vasprun.xml").get_band_structure("./KPOINTS", line_mode=True)
    run=BSVasprun(args.bandfile,parse_projected_eigen=True) #a special Class (BSVasprun), which is optimized for reading band structures and eigenvalues. 
    bands = run.get_band_structure("./KPOINTS", line_mode=True)

    # projected bands
    if not args.procar: #If no specific procar file was requested.
        print "Reading projected band structure from %s"%args.bandfile
        data= run.projected_eigenvalues
    else:
        print "Reading projected band structure from %s."%args.procar
        data = Procar(args.procar).data

    del run
   



    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(11.69, 8.27))
    tt=""
    for at in inclAtoms: 
        if at in aTypes: tt += "%s, " % at
    tt="Band diagram for "+tt[0:-2]+" atoms"
    #fig.suptitle("Bands diagram of graphene")
    fig.suptitle(tt)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # set ylim for the plot
    # ---------------------
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))

    if args.range=="all":
        emin -= bands.efermi + 1
        emax -= bands.efermi - 1
    else:
        w=args.range.replace('[','').replace("]","").split(":")
        try:
            emin=int(w[0])
            emax=int(w[1])
        except: 
            print "Error in energy range input. Please check your input. It should read -r -3:3 etc. Proceeding with the whole energy range." 
            emin -= bands.efermi + 1
            emax -= bands.efermi - 1
        #emin = -3
        #emax =  3

    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    if args.noEfcorr:
        dosrun.efermi=0.0
        bands.efermi=0.0

    if args.shift:
        dosrun.efermi-=args.shift
        bands.efermi-=args.shift


    # Band Diagram
    # ------------

    # sum up contribution over carbon atoms

    #data is of the format: spin: nd.array accessed with (k-point index, band index, ion index, orbital index)
    #data = data[Spin.up].sum(axis=2)
    
    data=dc(data[Spin.up])
    #print len(data),len(data[0]),len(data[0][0]),len(data[0][0][0])
    a=len(data)
    b=len(data[0])
    c=len(data[0][0])
    d=len(data[0][0][0])
    #print a,b,c,d
    data2=np.zeros((a,b,d))

    total=np.zeros((a,b)) #total contirbutions fore ach band at each Kpoint

    for i in range(a): #k-points
        for j in range(b): #bands
            #print data[i][j]
            #maxs=[]
            for k in range(c): #ions
                if aTypes[k] not in inclAtoms: continue #print aTypes[k];continue
                #maxs.append(max(data[i][j][k][:])**2)
                for l in range(d):#orbitals
                    #if aTypes[k]  in inclAtoms:
                    data2[i][j][l]+=data[i][j][k][l]
                    #total[i][j]+=data[i][j][k][l]**2

            #total[i][j]+=max(maxs)

    data=dc(data2)

    del data2


    # sum up px and py contributions and normalize contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):    
            sc = data[k][b][Orbital.s.value]**2
            pxpyc = data[k][b][Orbital.px.value]**2 + \
                data[k][b][Orbital.py.value]**2
            pzc = data[k][b][Orbital.pz.value]**2


            if rgb:
                tot = sc + pxpyc + pzc
                if tot==0.0:continue
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pxpyc / tot
                contrib[b, k, 2] = pzc / tot

            else:
                contrib[b, k, 0] = sc  #/ total[k][b]
                contrib[b, k, 1] = pxpyc  #/ total[k][b]
                contrib[b, k, 2] = pzc  #/ total[k][b]
    
    # plot bands using rgb mapping
    for b in range(bands.nb_bands): 
        if rgb:
            rgbline(ax1, #if 0 no color map. Only pz beads are plotted.
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]], #original
                #[e - dosrun.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0], #!!
                contrib[b, :, 1], 
                contrib[b, :, 2],1.0) #!! Normalized contibutions from the pz orbitals.
            ax1.scatter(range(len(bands.kpoints)),
                               [e - bands.efermi for e in bands.bands[Spin.up][b]],
                               color='b',marker="o", s=contrib[b, :, 2]*10 , alpha=1.0)
        else: 
            sz=150
            ax1.plot(range(len(bands.kpoints)),
                               [e - bands.efermi for e in bands.bands[Spin.up][b]],alpha=0.3,color='k')
            ax1.scatter(range(len(bands.kpoints)),
                               [e - bands.efermi for e in bands.bands[Spin.up][b]],
                               color='b',marker="o", s=contrib[b, :, 2]*sz , alpha=1.0)
            if not onlypz:
                ax1.scatter(range(len(bands.kpoints)),
                               [e - bands.efermi for e in bands.bands[Spin.up][b]],
                               color='g',marker="o", s=contrib[b, :, 1]*sz , alpha=1.0)
                ax1.scatter(range(len(bands.kpoints)),
                               [e - bands.efermi for e in bands.bands[Spin.up][b]],
                               color='r',marker="o", s=contrib[b, :, 0]*sz , alpha=1.0)

    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=3,linestyle="dashed")#!!!

    # labels
    nlabs = len(labels)
    #keep this 0 !!!
    if 0: #len(steps)>1: #if unequi spacing was used btw HighK points. (taken from KPOINTS file). VASP actually does not work with unequi spacing, when using auto KPOINT generation, one must give explicitly the KPOINT path.
        #step=0
        #ticks=[0]
        #for st in steps: step+=int(st);ticks.append(step)
        ticks=[0]
        tsteps=sum(steps); #[tsteps += int(st) for st in steps]
        con=float(tsteps)/float(len(bands.kpoints))
        step=0
        for st in steps: step+=int(int(st)/con); ticks.append(step)
        ax1.set_xticks(ticks)
        #print ticks, len(bands.kpoints)

    else: #If fixed spacing is given in KPOINTS file, then equidistance is used.
        step = len(bands.kpoints) / (nlabs - 1)
        for i, lab in enumerate(labels):
            ax1.vlines(i * step, emin, emax, "k")
        ax1.set_xticks([i * step for i in range(nlabs)])

    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, len(bands.kpoints))

    # Density of state
    # ----------------
    #nAtoms=float(len(dosrun.pdos))
    
    nedos=len(dosrun.pdos[0][Orbital.s][Spin.up]) #Grid size used for DOS calcultion.
    #print nedos
    total_s=np.zeros(nedos)
    total_px=np.zeros(nedos)
    total_py=np.zeros(nedos)
    total_pz=np.zeros(nedos)
    others=np.zeros(nedos)

    nAtoms=0
    for i in range(len(dosrun.pdos)): #iterate over all atoms and get the DOS/atom.
        if aTypes[i] in inclAtoms: #Only include desired atom types.
            nAtoms+=1
            total_s+=npa(dosrun.pdos[i][Orbital.s][Spin.up]) #/ nAtoms
            total_px+=npa(dosrun.pdos[i][Orbital.px][Spin.up]) #/ nAtoms
            total_py+=npa(dosrun.pdos[i][Orbital.py][Spin.up]) #/ nAtoms
            total_pz+=npa(dosrun.pdos[i][Orbital.pz][Spin.up]) #/ nAtoms
            #others+=
            #Add others here !! No Need !!
        else: continue

    total_s/=nAtoms
    total_px/=nAtoms
    total_py/=nAtoms
    total_pz/=nAtoms

    total_dos=total_s+total_px+total_py+total_pz#+others

    #For correct normalization of DOS in a given subregion (e.g. -3 ev to 3 eV).
    data=[]
    for i in range(len(dosrun.tdos.energies)):
        if emin<=dosrun.tdos.energies[i]- dosrun.efermi <=emax:
            #data.append([dosrun.tdos.energies[i]- dosrun.efermi,total_dos[i]])
            data.append(total_dos[i])

    #DOS plot settings.
#    max_dens=max(dosrun.tdos.densities[Spin.up])/nAtoms
    #max_dens=max(total_dos)
    max_dens=max(data)                    
    ax2.set_yticklabels([])
    ax2.grid()
    intv=float(max_dens)/4.0
    ax2.set_xticks(np.arange(0, max_dens, intv))
    #ax2.set_xticklabels(np.arange(0, max_dens, intv))
    ax2.set_xticklabels([])

    ax2.set_xlim(1e-6, max_dens)
    ax2.hlines(y=0, xmin=0, xmax=max_dens, color="k", lw=3,linestyle="dashed")
    ax2.set_xlabel("\nDOS/atom")

    # s contribution  
    ax2.plot( total_s,
             dosrun.tdos.energies - dosrun.efermi,
             "r-", label="s", linewidth=2)

   # px+py contribution
    ax2.plot( total_px+total_py ,
             dosrun.tdos.energies - dosrun.efermi,
             "g-",
             label="(px, py)",
             linewidth=2)

   # pz contribution
    ax2.plot( total_pz,
             dosrun.tdos.energies - dosrun.efermi,
             "b-", label="pz", linewidth=2)
   # others
    ax2.plot( np.zeros((nedos)),
             dosrun.tdos.energies - dosrun.efermi,
             "k-", label="others", linewidth=2)

    # total dos
    #ax2.fill_between(dosrun.tdos.densities[Spin.up]/len(aTypes),  #This would also include all orbitals from all atoms (selected and non-selected), like d etc.
    ax2.fill_between(total_dos,
                     0,
                     dosrun.tdos.energies - dosrun.efermi,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.7, 0.7, 0.7))

    #ax2.plot(dosrun.tdos.densities[Spin.up]/len(aTypes),
    ax2.plot(total_dos,
             dosrun.tdos.energies - dosrun.efermi,
             color=(0.6, 0.6, 0.6),
             label="total DOS")

    # plot format style
    # -----------------

    ax2.legend(fancybox=True, shadow=True, prop={'size': 14})
    plt.subplots_adjust(wspace=0)

    plt.show()
    #plt.savefig(argv[0].strip(".py") + ".pdf", format="pdf")
    plt.savefig("band-dos.png", format="png")

