#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import argparse,math
import ase
import ase.io
from ase.visualize import view
from os import system,popen,chdir,getcwd,getenv,putenv,listdir
import os.path
from re import search
from sys import exit
#from ase.build import *
import ase.build
from ase import Atoms,Atom
import fractions
from copy import deepcopy
import time
#import ase.build.tools
import matplotlib.pyplot as plt
import ase.calculators.castep,ase.calculators.vasp
from ase.calculators.vasp import Vasp,Vasp2 #requires ase 3.17.0
from spglib import find_primitive,standardize_cell,get_spacegroup #, niggli_reduce #niggli_reduce from ASE-tools conflicts with that from spglib.
from ase.build import niggli_reduce

from ase.geometry import find_mic  #min dist through PBC using the minimum-image covnention

from ase.constraints import FixAtoms
#from ase.spacegroup import crystal

#from asap3.analysis import CoordinationNumbers

from ase.build import cut


try:
        from os import popen4 #as popen
except:
        from os import popen #as popen

#Surface creation and alignment tool by James P. Darby.
from tools.builder import *

def get_thickness(atoms):
        zcoords=[at.position[2] for at in atoms]
        print (max(zcoords)-min(zcoords))
        return max(zcoords)-min(zcoords)

def vasp_continue(ifDel=0):
        fl = listdir(".")
        fl.sort()
        #print fl
        flag=0
        ind=0
        for fn in fl:
                if fn[0:3]=="RUN" :
                        ind = int(fn[-1])

        run="RUN%d"%(ind+1)
        print ("Previous data stored in %s."%run)
        system('mkdir %s'%run)
        system('cp * %s 2> /dev/null '%run)
        if ifDel: system('rm WAVECAR CHGCAR')
 
def compute_rdf(dists,cellVol,r_range=None, bin_width=0.005, n_bins=None, norm_by_V=True):
    """Compute radial distribution functions for interatomic distances along a given trajectory. Adapted from mdtraj pacakge. Returns, bin centre points, g(r) and the coordination(r)

    Parameters
    ----------
    dists: list or array, shape (1,N)
    r_range : array-like, shape=(2,), optional, default=(min(dists), max(dists))
        Minimum and maximum radii.
    bin_width : float, optional, default=0.005
        Width of the bins in nanometers.
    n_bins : int, optional, default=None
        The number of bins. If specified, this will override the `bin_width`
         parameter.

    Returns
    -------
    r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radii values corresponding to the centers of the bins.
    g_r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radial distribution function values at r.
    coord : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        coordination number at r.
    """

    dists=np.array(dists)

    if r_range is None:
        r_range = np.array([min(dists), max(dists)])
    #r_range = ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
    #                      shape=(2,), warn_on_cast=False)  #This is a mdtraj property as wel !!
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('n_bins must be a positive integer')
        bins = np.linspace(r_range[0], r_range[1], n_bins)
    else:
        bins = np.arange(r_range[0], r_range[1] + bin_width, bin_width)

    #distances = compute_distances(traj, pairs, periodic=periodic, opt=opt) #Replace this with ASE !!
    g_r, edges = np.histogram(dists, bins=bins)
    r = 0.5 * (edges[1:] + edges[:-1])


    if norm_by_V:
    # Normalize by volume of the spherical shell.
    # See discussion https://github.com/mdtraj/mdtraj/pull/724. There might be
    # a less biased way to accomplish this. The conclusion was that this could
    # be interesting to try, but is likely not hugely consequential. This method
    # of doing the calculations matches the implementation in other packages like
    # AmberTools' cpptraj and gromacs g_rdf. 
    # VMD and LAMMPS also seem to do this !!!
        V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        #norm = len(pairs) * np.sum(1.0 / traj.unitcell_volumes) * V
        norm = len(dists) * np.sum(1.0 / cellVol ) * V
        g_r = g_r.astype(np.float64) / norm  # From int64.
    else:
        norm = len(dists) 
        g_r = g_r.astype(np.float64) / norm  # From int64.

    coord=np.array([np.sum(g_r[0:i]) for i in range(len(g_r))])

    return r, g_r , coord

def minDist_ASE(arr1,arr2,latt): #needs np.array inputs
    if len(arr1)==4: arr1=arr1[1:4]
    if len(arr2)==4: arr2=arr2[1:4]
    D,x=find_mic(D=np.array([arr1-arr2]),cell=latt,pbc=True) 
    return x[0]


    
def minDist(arr1,arr2,latt):
    #Finds the minimum distance from the images applying the PBC. (Minimum image convention).
    #This function works fine for ortho and non-ortho cells, however the find_mic() from ase.geometry is upto 2x faster (using a different algorithm).
    if len(arr1)==4: arr1=arr1[1:4]
    if len(arr2)==4: arr2=arr2[1:4]
    #print arr1,arr2
    minD=dist(arr1,arr2)
    #newarr2=deepcopy(arr2)
    newarr2=[0.0, 0.0, 0.0]
    for i in range(-1,2,1):
        for j in range(-1,2,1):
            for k in range(-1,2,1):
                #newarr2=deepcopy(arr2)
                newarr2=[0.0, 0.0, 0.0]
                newarr2[0]=arr2[0]+i*latt[0][0]+j*latt[1][0]+k*latt[2][0]
                newarr2[1]=arr2[1]+i*latt[0][1]+j*latt[1][1]+k*latt[2][1]
                newarr2[2]=arr2[2]+i*latt[0][2]+j*latt[1][2]+k*latt[2][2]
                currD=dist(arr1,newarr2)
                if currD<minD: minD=currD ;#print minD
    return minD
def dist(arr1,arr2):
    #Distance btw two atoms. Works with two 1x3 arrays.
    if len(arr1)==4: arr1=arr1[1:4]
    if len(arr2)==4: arr2=arr2[1:4]
    return math.sqrt((arr1[0]-arr2[0])**2+(arr1[1]-arr2[1])**2+(arr1[2]-arr2[2])**2)

def misfit(slab1,slab2,ptol=2,ifPlot=False): #Computes the misfit parameter between two slabs. (see Surf. Interface Anal. 2003, 35, 835-841. )
        #Assuimng that two slabs are well aligned, overlap area should be equalt to smaller one's surface area.  !!FIND A MORE ACCURATE WAY OF DETERMINING OVERLAPPING SURFACE.
        #align_slab_axes(slab1,slab2,ptol)
        A1=surf_area(slab1);A2=surf_area(slab2)
        if ifPlot: plot_slab_axes(slab1,slab2)
        return 1 - 2*min(A1,A2)/(A1+A2)

def surf_area(slab1):
        return np.linalg.norm(np.cross(slab1.cell[0],slab1.cell[1]))

def volume(cell):
        return np.abs(np.dot(cell[2], np.cross(cell[0], cell[1])))

def get_fu(atoms):
        aTypes=atoms.get_chemical_symbols()
        atoms={}
        #print aTypes
        for at in aTypes:
            #print ln
            #at=ln.split()[0]
            if at not in atoms:atoms[at]=1
            else: atoms[at]+=1
        

        keys=list(atoms.keys())
        keys.sort()
        content=""
        for key in keys:
                content+=key
                if atoms[key]!=1:content+=str(atoms[key])

        #Determine the formula unit.

        try:
                fu=1
                vals=list(atoms.values())
                for i in range(2,min(vals)+1):
                        fl=1
                        for j in vals:
                                if j%i!=0:fl=0
                        if fl:fu=i
        #print fu
        except: print("Error in fu determination, fu set to 1");   fu=1
    
        return fu


                

def find_prim(atoms,tol=1e-4):#using SPGlib find primitive cell of a given atoms object.
        scaled_positions= atoms.get_scaled_positions()#(wrap=True) 
        cell=(atoms.cell, scaled_positions, atoms.numbers)
        #print cell
        print("Space group of the given cell using tolerance=%f: %s"%(tol,get_spacegroup(cell,symprec=tol)))
        lattice, scaled_positions, numbers = find_primitive(cell, symprec=tol)
        #lattice, scaled_positions, numbers = standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=tol)
        atoms2=Atoms(symbols=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
        return atoms2

def grep(key,fname,n=-1): 
    #Uses grep to get the nth line with the keyword from a given file. By default the last occurence is returned.
        try:
                return popen4('grep -m %d "%s" %s '%(n,key,fname),"r")[1].readlines()[-1][0:-1] #don't take \n at the end.  #Fastest one!!!
        except: 
                return popen('grep -m %d "%s" %s '%(n,key,fname),"r").read()[0:-1]
        else:
                return ""

def call_vasp(atoms,calc=None, typ="sp",xc='PBE',name='./VASP-tmp',PP='',param='opt.param',resDir="",dipolCorr=False,dipolDir='z',KPgrid="1 1 1",KPspacing="", ifPrint=False, ENCUT=300,ifManualRun=True,FixCell=False,FixVol=False,FixList=[],hubU={},ispin=0,restart=None,nelm=150,nsw=800,etol=1E-8,ediffg=-0.01,slowConv=0,ismear=0,sigma=0.01,passivate=None,gamma=1,algo="Fast",nosymm=0,exe=None):
        import ase.io

        if not exe: exe=getenv('VASP_COMMAND') ;# print ('bk')
        if not exe: exe='vasp_std'
        exe="mpirun -np 128 ~/APPS/vasp.5.4.4/bin/vasp_std > vasp.out&"
        os.environ['VASP_COMMAND']=exe
        print("VASP command: ", exe, os.environ['VASP_COMMAND'])

        cwd=os.getcwd()


        try: ncores=int(getenv('SLURM_NTASKS'))  #ncores=int(popen('echo "$SLURM_NTASKS"',"r").read()[0:-1]) #add ones for PBS
        except:ncores=1
        if len(atoms)<20: ncores=1 #NPAR>1 can be problematic with small systems due to partioning of the bands.
 

        if calc!=None:
                #This is good for assiging the settings from a previously created Calculator object (i.e. not to repeat all the settings the second time). None of the other options will be taken into account.
                atoms.set_calculator(calc)
                #manual VASP run icin ayri bir function yap ve burda ve asagida cagir.
                chdir(calc.directory)
                E=atoms.get_potential_energy()
                chdir('..')
                return E,atoms


        wDir=name
        os.system('mkdir -p %s'%wDir)
        chdir(wDir)

        asyms=[]
        for at in atoms:  #To guarantee the same order as in POSCAR.
                if at.symbol not in asyms:asyms.append(at.symbol)

        if 1 and os.path.exists('OUTCAR'):
                print("A previous calculation found in %s"%name)
                print(getcwd())
                restart=True

        try: 
                calc = Vasp2(atoms,restart=restart, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt=None) 
                print ("Previous run was read successfully, checking if converged...")
        except: #if restarting fails
                print("Problem with restarting from the previous run. Starting a new calculation.")
                calc = Vasp2(atoms,restart=None, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt=None);restart=False

        if restart: #if restarted no need to assign the variables and use the ones read from INCAR and OUTCAR by ASE.
                if Vasp2.read_convergence(calc):#calc.converged:#read_convergence('OUTCAR'):
                        print ("Calculation has converged, reading data from vasprun.xml/OUTCAR... ")
                        try:atoms=ase.io.read('vasprun.xml',index=-1)
                        except:atoms=ase.io.read('OUTCAR',index=-1)
                else:
                        #Standard ASE impletation for VASP-restart does not read the last geom for some reason, this is a simple trick to get the CONTCAR and same settigns from the prevoious run.
                        print ("Not converged, restarting from last point")
                        atoms=ase.io.read('CONTCAR') #no energy so must be repeated.
                        atoms.set_calculator(calc)
                        vasp_continue() #copies prev. run  files into RUNX folder.
                        #atoms.write('POSCAR',format='vasp',vasp5=1)

                E=atoms.get_potential_energy()
                chdir('..')
                return E,atoms



        if typ.lower()=='geom':IBRION=2 #CG algo
        elif typ.lower()=="sp": IBRION=-1;nsw=0

        if FixCell:isif=2 #fixed cell optimization
        elif FixVol:isif=4 #cell dimensions can vary but total volume is kept fixed (good for keeping the vacuum padding).
        else:isif=3 #full cell relaxation.
        if args.vac<=2:isif=3 #this must be set to relax the cell only in the vacuum direction (i.e. c-axis), but not possible in VASP.

        #calc.initialize()
        ncore=1
        npar=int(np.sqrt(ncores))
        if ncores%npar!=0 or npar==1:npar=8
        if len(atoms)<20:ncore=16 #but NPAR should be unset/commented, this is neded for the Sub-Space-Matrix is not hermitian in the Davidson algorithm.
        ncore=16 #only 8 works for small systems like 2-atom graphene.
        calc.set(xc=xc.lower(),encut=ENCUT,prec="Normal",ediff=etol, ediffg=ediffg,sigma=sigma,reciprocal=0,algo=algo,ispin=ispin,lreal="AUTO",nelm=nelm,ibrion=IBRION,gamma=gamma,isif=isif,nsw=nsw,ismear=ismear,npar=npar) #ncore=ncore)#,npar=npar)
        calc.set(nwrite=1, lcharg=0 , lwave=0) #Wavecar could be used for next run, specially for ISPIN=2
        if ispin==2:calc.set(lwave=1) 

        calc.set(lmaxmix=6,lasph=1) #need when there is a trans metal.

        if nosymm: calc.set(isym=0)
        #TODO add support for reading INCAR file if exists.
        #calc.set(setups='recommended')

        if KPspacing: 
                calc.set(kspacing=KPspacing,kgamma=gamma)#calc.cell.kpoints_mp_spacing = str(KPspacing) #default=0.05 2*pi/A
                print ("Multiple KPOINTS: VASP standard implementation (vasp_std) is activated.")
                cmd=exe.split('vasp_')[0]+"vasp_std"
                calc.command=cmd
        else: #KPgrid given
                #calc.cell.kpoint_mp_grid = KPgrid #def='1 1 1'
                kpgrid=[int(x) for x in KPgrid.split()]
                import ase.dft.kpoints
                kpts = ase.dft.kpoints.monkhorst_pack(kpgrid) #+  [1./2./kpgrid[0],1./2./kpgrid[1],1./2./kpgrid[2]] #for placing Gamma point in the center.  #CHECK if works CORRECTLY!!
                calc.set(kpts=kpts)
                if KPgrid=='1 1 1': 
                        print ("VASP gamma-only implementation (vasp_gam) is activated.")
                        cmd=exe.split('vasp_')[0]+"vasp_gam"
                else:
                        print ("Multiple KPOINTS: VASP standard implementation (vasp_std) is activated.")
                        cmd=exe.split('vasp_')[0]+"vasp_std"

                calc.command=cmd
                #print(cmd)

        a,b,c,_,_,_=atoms.get_cell_lengths_and_angles()
        if a >35 or b>35 or c>35: calc.set(amin=0.01) #gor super large cells to avoid the charge sloshing along the long lattice vector. 
        if slowConv:
                calc.set(amix=0.1,bmix=3.0,lmaxmix=6)

        if args.vac >2 and (dipolCorr and dipolCorr.lower()!='none'):
                calc.set(ldipol = ".TRUE.") #this one switches on  the corrections to the potential and thus forces
                if dipolDir=='x':idipol=1 #IDIPOL switches on the single point energy corrections.
                elif dipolDir=='y':idipol=2
                elif dipolDir=='z':idipol=3
                calc.set(idipol = idipol) #no need to set DIPOL keyword (giving the charge center, if no charged dipole in the system (i.e. total charge=0) 

                #One needs to set the DIPOL keyword (i.e. the center of the net charge in scaled/fractional coordinates), as the automated algorithm leads to terrible convergence. It should be the centre of charge (that requires the CHG file to be created and analysed to be determined), instead  the centre of mass could be used as a good approximation. Setting DIPOL helps even the veryslow convergence when the LDIPOL=true. It even affects the total spin of the system (when ISPIN=2).

                calc.set(dipol=atoms.get_center_of_mass(scaled=1)) #Does this need to be updated for continuation run or better to stick with initial value for better covnergence ??? (CHECK)


        if 1 and len(hubU)!=0: #is not None:
                calc.set(ldau=".True.",ldauprint=0,ldautype=2,lmaxmix=6)
                ats=sorted(hubU.keys())
                #asyms=np.unique(atoms.get_chemical_symbols())
                #asyms=popen('head -n1 POSCAR',"r").read()[0:-1].split() #this doen't work always well.

                #if len(asyms)==0: print ("DFT-U Warning: Atomic info cannot be read from the head line of POSCAR")
                ldauu=[];ldaul=[]
                for at in asyms:
                        if at in ats: #continue
                                if len(hubU[at].keys())>1: print("WARNING: cannot set multiple U values for the same atom type (%s), setting to 0... "%at);ldauu.append(0);ldaul.append(0);continue
                                #print (hubU[at].keys()[0])
                                orb="".join(hubU[at].keys())#[0] indexing in dict_key format is a problem in Python3.
                                ldauu.append(float(hubU[at][orb]))
                                if orb=="p":ldaul.append(1)
                                elif orb=="d":ldaul.append(2)
                                elif orb=="f":ldaul.append(3)
                        else:ldaul.append(0);ldauu.append(0)
                calc.set(ldauj=[ 0 for i in range(len(asyms)) ])
                calc.set(ldauu=ldauu,ldaul=ldaul)
                calc.set(lasph=1)#This is essential for accurate total energies and band structure calculations for f-elements (e.g. ceria), all 3d-elements (transition metal oxides), and magnetic atoms in the 2nd row (B-F atom), in particular if LDA+U or hybrid functionals or meta-GGAs are used, since these functionals often result in aspherical charge densities.

        if passivate: #compute the coordination numbers for each atom, do statistics and then add passivating (pseudo)-Hydrogen atoms. 
                #from ase.build import molecule
                
                #Preparations
                setups={'base': 'recommended'}
                cov_radii={}
                for ln in open(os.environ['HOME']+"/covalent_radii.dat",'r'):
                        x=ln.split()
                        cov_radii[x[0]]=float(x[1])

                #Get the valance electron for each atom type from corr. POTCAR file (ZVAL)


                from ase.neighborlist import NeighborList,neighbor_list
                from ase.utils import natural_cutoffs #distance cutoff based on covalent radii.
                #ASEs neighborlist class has bugs, it counts number of coordination/neighbors for atoms located at the 0.0, 0.5 and 1.0 in fractional coordinates (a,b,c does not matter) !! Use rdf calculation instead. Actually this is a bug in natural_cutoff which cannot capture the different oxidation states of e.g. Fe. It's defined too short for capturing Fe6+.
                #natural_cutoff of ASE does not work properly for finding the orrect no of neighbours.
                cov_cutoff={}
                for asym in asyms:
                        for asym2 in asyms:
                                key=(asym,asym2)
                                if key not in cov_cutoff: cov_cutoff[key]=(cov_radii[asym]+cov_radii[asym2])*1.1
                if 0: #This does not work properly as one cannot deine the pair distances explicitly for each atom type pair, instead a radius defiend for each atom type (so same as natural_cutoffs) . USE neighbour_list instead.
                        #nl = NeighborList(cutoffs=natural_cutoffs(atoms), skin=0.3, sorted=False, self_interaction=0, bothways=1)
                        #nl = NeighborList(cutoffs=cov_cutoff, skin=0.3, sorted=False, self_interaction=0, bothways=1)
                        nl = NeighborList(cutoffs=[cov_radii[at.symbol] for at in atoms], skin=0.3, sorted=False, self_interaction=0, bothways=1)
                        nl.update(atoms)
                        #print(nl.npbcneighbors)
                        coord=[]
                        for at in atoms:   coord.append(len(nl.get_neighbors(at.index)[0]))
                else:  # bothways is on by default in this. This works well. #THIS DOES NOT WORK WELL EITHER, USE OVITO ISNTEAD (GIVES THE BEST RESULTS).
                        i,j=neighbor_list('ij',atoms,cutoff=cov_cutoff)#natural_cutoffs(atoms)) #there is a bug for determining the neighbours for atoms located at the origin. #i: index of the central atom, j: index of the neighbour.
                        coord = np.bincount(i) #Nx1 ; coutns no of occurences of each value in the input array. #one could use np.unique() for counts as well.

                        #unique, counts = np.unique(i, return_counts=True)
                        #print(unique, counts)   

                        #The dictionary of neighbour types for each atom index.
                        ntypes={}
                        tmp=[]
                        pkey=""
                        for k in range(len(i)):
                                #print(i[k],j[k])
                                key=str(i[k])
                                #if key not in ntypes: 
                                if k!=0 and key!=pkey: ntypes[pkey]=tmp;tmp=[];pkey=key#[atoms[j[k]].type]
                                elif k==0:pkey=key
                                elif k==len(i)-1: ntypes[key]=tmp
                                tmp.append(atoms[j[k]].symbol)

                        for k in sorted(ntypes.keys()): 
                                #k =int(k)
                                print (k,ntypes[k],coord[int(k)])
                        #print (sorted(ntypes))
                  
                        #unique, counts = np.unique(ntypes.values(), return_counts=True)
                        #print(unique, counts)

                print(coord,len(coord))

                #exit()
   
                #Get layers from the structure based on z coordinates. #Use 0.7 A bin_width
                zcoords=[at.position[2] for at in atoms]
                bin_width=0.5 #Looks like a reasonable binning but Determine this from the avg distances btw atoms in z direction !!  original:0.7
                n_bins=None
                r_range = np.array([min(zcoords), max(zcoords)])
                if n_bins is not None:
                        n_bins = int(n_bins)
                        if n_bins <= 0:
                                raise ValueError('n_bins must be a positive integer')
                        bins = np.linspace(r_range[0], r_range[1], n_bins)
                else:
                        bins = np.arange(r_range[0], r_range[1] + bin_width, bin_width)
                #plt.hist(zcoords,bins=bins)
                #plt.show()
                counts,bin_edges=np.histogram(zcoords,bins=bins) #the layers in terms of c-axis sepration.

                #Instead of histogram of zcoords, compare the zcoords of each atom and classify the based on being in close vicinity. More precise than histogram in determining distinct layers.
                layers={}
                tol=1e-2 #tolerance for determining if atom is in a partiular layer. in A
                #for at in atoms:
                #flag=0
                for i,zc in enumerate(zcoords):
                        at=atoms[i]
                        #if i==0: layers[zc]=[i]
                        keys=layers.keys()
                        flag=0
                        for key in keys:
                                if abs(zc-key)<tol: 
                                        layers[key].append(i);flag=1;break
                        if not flag:
                                layers[zc]=[i]

                for key in sorted(layers.keys()):                     print (key,layers[key])
                layers=[layers[key] for key in sorted(layers.keys())]
                print (layers)

                #Determine the matching top and layers for ensuring symmetric slab termination
                #Move layer determination to  a function.
                #Do also the reverse way (so matching topmost layer to bottom layer).
                tol=1e-5
                for i in range(4):
                        #flag=0
                        #at1=atoms[i]
                        lay1=layers[i]
                        for j in range(-1,-5,-1):
                                #flag=0
                                lay2=layers[j]
                                #at2=atoms
                                if len(lay1)!=len(lay2): continue #does not match
                                for k in range(len(lay1)):
                                        pos1=atoms[lay1[k]].position  ; pos2=atoms[lay2[k]].position
                                        if abs(pos1[0]-pos2[0]) > tol or abs(pos1[1]-pos2[1]) > tol: break#layers do not match.
                                #flag=1
                                print ("Layer #%d matches layer #%d"%(i,j))


                #Max occurence of coordination no for a given atom type can also be used, instead of avg coord no.
                if 1:
                    crds={};common_crd={};nty={};common_ntypes={};valance={};dangling={}
                    for aty in asyms:
                            #Determine the no of coordination for each atom.
                            crds[aty]=[coord[at.index] for at in atoms if at.symbol==aty]
                            common_crd[aty]=np.argmax(np.bincount(crds[aty]))
                            #Determine the neighbour types.
                            nty[aty]=[ntypes[str(at.index)] for at in atoms if at.symbol==aty] #Do we ned to store this??
                            unique, counts = np.unique(nty[aty], return_counts=True)
                            counts, unique =zip(*sorted(zip(counts, unique),reverse=1))
                            #print(unique, counts)
                            common_ntypes[aty]=unique[0]
                            #Get the valence electron no from the POTCAR files.
                            potfile="%s/potpaw/%s/POTCAR"%(getenv('VASP_PP_PATH'),aty)
                            with open(potfile,'r') as f:
                                    for ln in f:
                                            if  'ZVAL' in ln:  val=float(ln.split('ZVAL')[-1].split()[1]);  break
                            if val>8: val-=8
                            valance[aty]=val
                            dangling[aty]=val/common_crd[aty]
                            print (aty,common_crd[aty],common_ntypes[aty],valance[aty],dangling[aty])
                            
                    data=common_crd
                else:
                    #Get the average coord no for each atom type.
                    data={}
                    for i,at in enumerate(atoms): #Get stats on coordination no
                            #i=at.index
                            typ=at.symbol
                            if typ in data: data[typ]+=coord[i]
                            else: data[typ]=coord[i]


                    for key in asyms: #atom types; 
                            data[key]=data[key]/float(len([at for at in atoms if at.symbol==key]))
                            print (key,data[key])

                scpos=atoms.get_scaled_positions()
                undercoord=[]; overcoord=[]; normal=[]
                for i,at in enumerate(atoms):
                        if coord[i]<data[at.symbol]:undercoord.append(i)
                        elif coord[i]>data[at.symbol]:overcoord.append(i)
                        else:normal.append(i)
                        #print (i,at.symbol,at.position,scpos[i],coord[i])

                print ("Undercoordinated:",undercoord)
                print("Overcoordinated: ",overcoord)
                print("Standard coordination:",normal)
                #print(len(bin_edges))

                #Decide on which layers to passivate (top or bottom)
                if passivate.lower()=="bot":slayers=[0,1,2] #layer numbering starts from 0 
                elif passivate.lower()=="top":slayers=[-1,-2,-3]

                #Switch to layers determined earlier instead of using zcoords histogram.
                Hcnt=0
                for sl in slayers:
                        layer=layers[sl]
                        for i in layer:
                                if i not in undercoord: continue
                                at=atoms[i]
                                crd=coord[i]

                                print (at,crd, sl)

                                Hcnt+=1
                                offset=(cov_radii[at.symbol]+cov_radii['H'])*1.0
                                if passivate=='top': offset=-offset
                                atoms.extend(Atom("H",(at.position[0],at.position[1],at.position[2]-offset)))
                                #Determine the amont of missing e-.
                                missing=[];zval=0.0
                                c1=deepcopy(ntypes[str(i)])
                                c2=deepcopy(common_ntypes[at.symbol])

                                for c in c2:
                                        flag=0
                                        for cc in c1:
                                                if c ==cc: flag=1;c1.remove(cc);break
                                        if not flag: missing.append(c);zval+=dangling[at.symbol]
                                print ("Missing atoms: %s No of e- on pseudohydrogen: %.1f "%(missing,zval))
                                setups[len(atoms)-1]='H%.1f'%zval #when applied the added pseuodo-H atom is moved to the top of the atom list. 
                                #this can be done alltogether to prevent fractioning of the pseudohydrogens and limit number of atom enetries in POSCAR.
                                #print (setups)


                """
                #Using zcoords histogram.
                #Decide on which layers to passivate (top or bottom)
                if passivate.lower()=="bot":slayers=[1,2,3] #bin numbering starts from 1 
                elif passivate.lower()=="top":slayers=[len(bin_edges)-1, len(bin_edges)-2,len(bin_edges)-3]

                Hcnt=0
                for i in undercoord: #can be combined with the above atoms loop.
                        at=atoms[i]
                        crd=coord[i]
                        for j,be in enumerate(bin_edges): #Zcoord histogram.
                                if at.position[2] <be:
                                        if j  in slayers: #This check is needed for avaoiding passivation of  a defect in the middle of the slab rather than the surface, and passivating only the target surface.
                                                #TODO:check fro clashes with already-existing  atoms.
                                                print (at,crd, j)
                                                Hcnt+=1
                                                offset=(cov_radii[at.symbol]+cov_radii['H'])*1.0
                                                if passivate=='top': offset=-offset
                                                atoms.extend(Atom("H",(at.position[0],at.position[1],at.position[2]-offset)))
                                                #Determine the amont of missing e-.
                                                missing=[];zval=0.0
                                                c1=deepcopy(ntypes[str(i)])
                                                c2=deepcopy(common_ntypes[at.symbol])
                                                
                                                for c in c2:
                                                        flag=0
                                                        for cc in c1:
                                                                if c ==cc: flag=1;c1.remove(cc);break
                                                        if not flag: missing.append(c);zval+=dangling[at.symbol]
                                                print ("Missing atoms: %s No of e- on pseudohydrogen: %.1f "%(missing,zval))
                                                setups[len(atoms)-1]='H%.1f'%zval #when applied the added pseuodo-H atom is moved to the top of the atom list.
                                                #print (setups)
                                        else:
                                                print ('other surface:',at,crd, j)
                                        break
                """
                print (setups)
                calc.set(setups=setups)
                print ("\nAdded %d (pseudo)hydrogens to saturate the dangling bonds."%Hcnt)
                print(atoms)
                #atoms=ase.build.sort(atoms)
                #atoms.sort()
                if args.view:view(atoms)

                if dipolCorr and dipolCorr.lower()!='none':  calc.set(dipol=atoms.get_center_of_mass(scaled=1))


        atoms.set_calculator(calc)      
        

        if args.dry: #do not run VASP calculation.
                chdir('..')
                return 0.0,atoms

        #TODO: Call VASP manually for doing the geom opt and read the output for the final enerfy and geom.
        #atoms.calc.initialize(atoms) #this does not create the input files for a manual run.
        atoms.get_potential_energy()
        #atoms.calc.set(command="");atoms.get_potential_energy();atoms.calc.set(command=exe)


        #no need to run it manually as first get_energy does run the geom opt as well, depending on the IBRION.
        #system(exe) 

        #import ase.io
        #atoms=ase.io.read('OUTCAR',format="vasp",index=-1)
        atoms=ase.io.read('vasprun.xml',index=-1)
        #atoms=ase.io.read(wDir+'/OUTCAR',index=-1)
        
        chdir('..')
        return atoms.get_potential_energy(),atoms
        

def call_castep(atoms,calc=None, typ="sp",PP='',wDir='./CASTEP-tmp',name='try',param='opt.param',resDir="",dipolCorr=False,dipolDir='z',KPgrid="1 1 1",KPspacing="", ifPrint=False,ifDryRun=False,ENCUT=0,ifManualRun=True,FixCell=False,FixList=[],hubU={},slowConv=0):
    #For applying constraints (FixList) atom numbering starts from 1.

        
    #exe="mpirun -n 4 castep";PP_path='/rscratch/bk393/pspots/CASTEP'
    #exe="mpirun -n 20 castep";PP_path='/u/fs1/bk393/pspots/CASTEP'

    #system("export CASTEP_COMMAND='%s'"%exe)
    #system("export CASTEP_COMMAND='mpirun -n 4 castep'")
    #system("export CASTEP_PP_PATH='%s'"%PP_path)

    exe=popen('echo "$CASTEP_COMMAND"',"r").read()[0:-1]
    PP_path=popen('echo "$CASTEP_PP_PATH"',"r").read()[0:-1]

    if calc!=None:
            #This is good for assiging the settings from a previously created Calculator object (i.e. not to repeat all the settings the second time). None of the other options will be taken into account.
            atoms.set_calculator(calc)
            #manual CASTEP run icin ayri bir function yap ve burda ve asagida cagir.
            return atoms 

    calc = ase.calculators.castep.Castep()

    #Assign the environment variables.
    calc._castep_command=exe
    calc._castep_pp_path=PP_path
    
    # include interface settings in .param file
    calc._export_settings = True

    # reuse the same directory
    calc._directory = wDir
    calc._rename_existing_dir = False
    calc._label = name

    
    if param:
        #Read paramters from param file input.
        calc.merge_param(param)
    else:        
        # Use default param settings (depreceated)
        calc.param.xc_functional = 'PBE'
        calc.param.cut_off_energy = 100 #500
        calc.param.num_dump_cycles = 0
        calc.param.geom_method = "lbfgs"
        calc.param.geom_max_iter= 10
        calc.param.write_cell_structure=True
        calc.param.spin_polarised=False
        calc.param.opt_strategy="speed"
        calc.param.mix_history_length=20
        calc.param.max_scf_cycles=100
        calc.param.calculate_stress=True
        #calc.param.finite_basis_corr=0
    
    # Cell settings
    #
    calc.cell.symmetry_generate=True
    calc.cell.snap_to_symmetry=True
    if KPspacing:calc.cell.kpoints_mp_spacing = str(KPspacing) #default=0.05 2*pi/A
    else: 
            calc.cell.kpoint_mp_grid = KPgrid #def='1 1 1'
            kpgrid=[float(x) for x in KPgrid.split()]
            calc.cell.kpoints_mp_offset='%.4f %.4f %.4f'%(1./2./kpgrid[0],1./2./kpgrid[1],1./2./kpgrid[2]) #for placing Gamma point in the center.

    
    #calc.cell.fix_com = False
    if FixCell: calc.cell.fix_all_cell = True
    if len(FixList)!=0:
            #c = FixAtoms(indices=[atom.index for atom in atoms if atom.symbol == 'Cu'])
            #c = FixAtoms(indices=FixList)
            #atoms.set_constraint(c) #This only work if the CASTEP is called by ASE (not for manual runs).

            str1=[]
            i=1
            for at in FixList:
                    atom=atoms[at]
                    #for j in range(1,4):
                    str1.append("%d %s %d %.8f %.8f %.8f"%(i,atom.symbol,atom.index,1,0,0))
                    str1.append("%d %s %d %.8f %.8f %.8f"%(i+1,atom.symbol,atom.index,0,1,0))
                    str1.append("%d %s %d %.8f %.8f %.8f"%(i+2,atom.symbol,atom.index,0,0,1))
                    i+=3

            calc.cell.ionic_constraints=str1 #a list object needed as input.
            calc.cell.snap_to_symmetry=False
            calc.cell.symmetry_generate=False

    if len(hubU)!=0: #is not None:
            ats=sorted(hubU.keys())
            #str1="%BLOCK HUBBARD_U \n  eV"
            str2=["eV"]
            asyms=atoms.get_chemical_symbols()
            for at in ats:
                    if at not in asyms: continue
                    str1="%3s "%at
                    for orb in sorted(hubU[at].keys()):
                            str1+="%2s: %s "%(orb,hubU[at][orb])
                    str2.append(str1)
            #str1+='\n %ENDBLOCK HUBBARD_U'
            #print(str1)
            calc.cell.hubbard_u=str2
     
    #This overwrites the task paramter from the param input.
    if typ.lower()=="sp":    calc.param.task = 'SinglePoint'
    elif typ.lower()=="geom":calc.Task = 'GeometryOptimization'
    
    if dipolCorr: #def: No dipole corrections. 
        if dipolCorr.upper()=="SC": calc.param.dipole_correction= "SELFCONSISTENT"
        elif dipolCorr=="static": calc.param.dipole_correction= "static"
        else: calc.param.dipole_correction= "None"
        calc.param.dipole_dir=dipolDir #options: x,y,z and a (only energy-corr)

    if slowConv:
            calc.param.mix_charge_amp=0.1
        
    
    #calc.initialize()#Creates all needed input in the _directory. (Not needed for actual run.)
    atoms.set_calculator(calc)  #Set for the previously created interface

    if ENCUT!=0: calc.param.cut_off_energy=ENCUT
    
    if PP!="":#otherwise on-the-fly(OTF) is used as default
        fnc=str(calc.param.xc_functional).split("\n")[1].upper()
        #print fnc
        PP=PP.upper()
        #print PP.lower().find(str(calc.param.xc_functional).lower())
        if PP != "OTF" and  PP.find(fnc)== -1:
                print("There is a problem with the pseudopotential choice. \nSelected PP does not match with XC functional being used: ",PP,fnc)
                exit()
        elif PP=="OTF": None #default is OTF anyway.
        else: atoms.calc.set_pspot(PP)  #This automatically sets the pseudo-potential for all present species to <Species>_<library>.usp. Make sure that _castep_pp_path is set correctly in the shell.


                            

    # Or Read all settings from previous calculation
    if resDir != "": #Does not read from checkpoint file for some reason, needs to be checked !!
        # Reset to CASTEP default
        #atoms.calc.param.task.clear()
        atoms = ase.io.castep.read_seed('%s/%s' % (wDir,name))
        calc.param.reuse=True
        calc.param._try_reuse=True
        calc._export_settings = True
        #print atoms.calc,"bora\n\n",calc
        
    if ifPrint:print (calc) #prints calculation summary.

    # necessary for tasks with changing positions
    # such as GeometryOptimization or MolecularDynamics (This option does not work as deisgnated !!! The atomic coords are not updated at the end of the geom. opt. unlike the energy.
    calc._set_atoms = True
    atoms.calc._set_atoms = True
    
    
    # Check for correct input
    if not ifDryRun:
            if ifManualRun: #If a manual run of the CASTEP is needed.
                    #str1="%s %s/%s"%(exe,wDir,name)
                    str1="%s %s"%(exe,name)  
                    #Add here the run3 option !
                    print("Running ",str1)

                    calc._copy_pspots=True
                    calc.initialize()

                    chdir(wDir)
                    
                    system(str1) #PSPOT missing in the folder
                    #x=parseCASTEP("%s/%s.geom"%(wDir,name),atoms=atoms)
                    task=str(atoms.calc.param.task).split()[-1]
                    print(task)
                    if task=='SinglePoint' : #use <seed>.castep file
                            x=parseCASTEP("%s.castep"%(name),atoms=atoms)
                    elif task=='GeometryOptimization': #use <seed>.geom file.
                            try:x=parseCASTEP("%s.geom"%(name),atoms=atoms) #in case no ptimisation is done (that could happen for 2D systems, where the slab is equivalent to the bulk str, so no geom opt is needed/possible.)
                            except:x=None

                            if x == None: x=parseCASTEP("%s.castep"%(name),atoms=atoms);#print ("bk",atoms)

                            if x[-2]==False: print("parseCASTEP: WARNING: Geometry optimization in %s.geom is not converged!"%name)
                    else:
                            print("parseCASTEP: ERROR: Calculation type is not supported.")
                            x=None
                    chdir("..")
                    return x
                    
            else: #CASTEP calculation is not done here. It will be called in the main script, when trying to reach the attributes, e.g. atoms.get_potential_energy().
                    return atoms
    else:
            if calc.dryrun_ok():

                    return atoms
            else:
                    print("CASTEP run: Found error in input")
                    print((calc._error))
                    return None
            
def parseCASTEP(fname,atoms=None):  #TODO: fix wrong atomic order (if the atoms object is ordered).
        #Read the CASTEP output to retrieve final 0K energy, atomic coords (Cart. and fractional), and forces.
        bohr2ang=0.52917721
        Ha2eV=27.211386
        atP2GPa=0.29421912e5 #P in atomic units to GPa.
        E=0.0;H=0.0 #Total Energy and Enthalpy
        h=[] #unit cell vectors
        s=[] # stress tensor
        xyz=[] #Cartesian atomic coords.
        forces=[] #forces in au (i.e. Ha/bohrradius).
        ifConv=0 #if the geometry is converged.
        fract=[]
        ids=[] #atomic id's

        print("Parsing ",fname)
        if fname.split(".")[-1]=="geom":
                try: tSteps=len(popen("grep '<-- c' %s"%fname).readlines())-1
                except:tSteps=0

                flag=0
                for ln in open(fname,'r'):#.readlines(): #readlines is outdated and slow.
                        ln=ln[0:-1]
                        if not flag and search("<-- c",ln):
                                if int(ln.split()[0])==tSteps:
                                        #print 'here'
                                        flag=1
                                        if search("T   T   T   T",ln): ifConv=1

                        elif flag:
                                if search("<-- E",ln):E=float(ln.split()[1])*Ha2eV; H=float(ln.split()[1])*Ha2eV
                                elif search("<-- h",ln): h.append([float(i)*bohr2ang for i in ln.split()[0:3]])
                                elif search("<-- S",ln): s.append([float(i)*atP2GPa for i in ln.split()[0:3]])
                                elif search("<-- R",ln): 
                                        x= ln.split(); 
                                        xyz.append([float(i)*bohr2ang for i in x[2:5]]); 
                                        ids.append(x[0]); 
                                elif search("<-- F",ln): forces.append([float(i)*(Ha2eV/bohr2ang) for i in ln.split()[2:5]])
                                
                                
       
        elif fname.split(".")[-1]=="castep":
                tSteps=0;latFlag=0;forces=[];s=[];h=[];latFlag=0;fract=[];ids=[]
                for ln in open(fname,'r'):#.readlines(): #readlines is outdated and slow.
                        ln=ln[0:-1]
                        #if ln=="":continue
                        if len(ln)<=2: continue
                        elif "Unit Cell" in ln: #for getting only the last one
                                tSteps+=1
                                forces=[];s=[];h=[];latFlag=0;fract=[];ids=[]
                        elif search("Final free energy \(E-TS\)    =",ln):H=float(ln.split("=")[1].split()[0])
                        elif search("Final energy, E",ln): E=float(ln.split("=")[1].split()[0])
                        elif search("\*  x ",ln) or search("\*  y ",ln) or search("\*  z ",ln): s.append([float(i) for i in ln.split()[2:5]]) #Stress tensor already in GPa.
                        elif len(ln)>2 and ln[1]=="\*" and len(ln.split())==7:
                                forces.append([float(i) for i in ln.split()[2:5]]) #already in eV/A
                        elif "            x" in ln and len(ln.split())==7:
                                fract.append([float(i) for i in ln.split()[3:6]])
                                ids.append(ln.split()[1])
                        elif "Real Lattice" in ln: latFlag=1;continue
                        elif "Lattice parameters" in ln: latFlag=0
                        elif  "Geometry optimization completed successfully." in ln: ifConv=1

                        if latFlag:  
                                if ln=="":latFlag=0;continue
                                h.append([float(i) for i in ln.split()[0:3]])
                                

                #Assuming it is a SP calculation (see call_CASTEP() function), so initial positions are not changed.
                
                #if len(h)==0: print("Hey !!"); h=atoms.get_cell()#.tolist()  #if no Lattice infoin castep file (not likely).
                
                #xyz=atoms.get_positions()#.tolist()
                
        if len(ids)==0: 
                print(ids);print ("parseCASTEP: WARNING: No chemical symbol info in %s."%fname); return None
        elif len(h)==0: 
                print ("parseCASTEP: WARNING: No cell info in %s."%fname); return None
        elif len(xyz)==0 and len(fract)==0: 
                print ("parseCASTEP: WARNING: No Cartesian or fractional coordinate info in %s."%fname); return None
        
                                                
        if atoms != None:
                #print(ids,h,atoms.calc,xyz,fract)
                #print (atoms.get_chemical_symbols())
                atoms.set_chemical_symbols(ids)
                atoms.set_cell(np.array(h))
                if len(xyz)!=0:atoms.set_positions(np.array(xyz))
                elif len(fract)!=0:atoms.set_scaled_positions(np.array(fract))

                #atoms.set_positions(np.array(xyz))
                
                
                #atoms2=Atoms(symbols=ids,cell=h,pbc=True,calculator=atoms.calc)#positions=xyz)

        else:
                atoms=Atoms(symbols=ids,cell=h,pbc=True,calculator=None)#positions=xyz)
                #return atoms

        
        #if len(xyz)!=0:atoms2.set_positions(np.array(xyz))
        #elif len(fract)!=0:atoms2.set_scaled_positions(np.array(fract))

        #Energy, Enthalpy, stress tensor, Cart coords, fractional coords, atoms object.
        return E, H, s, fract, forces, ifConv, atoms
     
def parseCASTEP_ASE(fname,atoms=None):
        atoms=ase.io.read(fname,index=-1)
        return atoms.get_potential_energy(),atoms #this gets the final enthalpy rather than E-0.5TS
   
def conv_layers(atoms,ifPlot=False,ifPrim=False):#layer convergence test (Input atoms with a calc object). 
        #print "Convergence of E/atom vs. #layers"
        Etol=1e-2 #eV/atom
        Ftol=5e-2 #eV/Angstroem
        Estr=0.1 #GPa

        #Initial values
        E=[0]; F=[0]; S=[0]
        atoms_orig=atoms.copy()
        calc=atoms.get_calculator()
        
        #find the primitive cells to reduce comp. efforts.
        if ifPrim: atoms=find_prim(atoms);atoms.set_calculator(calc)
        
        atoms.center(vacuum=args.vac, axis=2)
        nAt=atoms.get_number_of_atoms()
        E.append(atoms.get_potential_energy()/nAt)
        i=1;layers=[1]
        while abs(E[i]-E[i-1]) > Etol:
                layers.append(1+1*i)
                atoms=atoms_orig.copy()
                atoms=atoms.repeat((1,1,layers[-1]))
                atoms.center(vacuum=args.vac, axis=2)
                atoms.set_calculator(calc)
                #if args.view:xview(atoms)
                nAt=atoms.get_number_of_atoms()
                E.append(atoms.get_potential_energy()/nAt)
                print("Iter. #%d, #layers: %d, #atoms: %d "%(i,layers[-1],nAt))
                print("deltaE: %.3e eV/atom; target: %.3e eV."%(abs(E[i]-E[i-1]),Etol))
                i += 1
        print("conv_layers: E/atom converged to %.2f eV with %d layers."%(Etol,layers[-1]))
        if ifPlot: #Do plotting of E/atom vs. #layers
                plt.plot(layers,E[1:], 'ko-')
                plt.xlabel('Number of layers')
                plt.ylabel('Energy per atom (eV/atom)')
                #plt.savefig('layer-conv.png')
                plt.show()
                
        return layers[-1],E[-1]*nAt,atoms               

def get_interface_energy(slab1,slab2,Ebulk1,Ebulk2,dist=1.0,convLayers=False): #Input is slab1,slab2 as Atoms object with a calculator property assigned.
        if not slab2.has(calc):
                calc=slab1.get_calculator()
                slab2.set_calculator(calc)
        
        if convLayers: #A layer convergence test is done here!.
                print("Layer convergence test is switched on.")
                #slab1.set_calculator(calc)
                slab1=call_castep(slab1,typ="SP",dipolCorr='SC',name='slab1',ENCUT=ecut,KPgrid='1 1 1',PP=pp) #Use SP for efficiecy.
                calc=slab1.get_calculator()       
                print("Working on slab 1.")
                nl1,Eslab1,slab1=conv_layers(slab1)#,ifPrim=1)  #Vacuum layer is added automatically within the function.
                fu1=get_fu(slab1)
                
                #slab2=call_castep(slab2,typ="SP",dipolCorr='SC',name='slab2',ENCUT=ecut,KPgrid='1 1 1',PP=pp)
                
                slab2.set_calculator(calc)
                slab2.calc._label="slab2"    #change name for slab2
                print("Working on slab 2.")                
                nl2,Eslab2,slab2=conv_layers(slab2)#,ifPrim=1)

                fu2=get_fu(slab2)
                print(nl1,nl2,Eslab1,Eslab2,fu1,fu2)

        else:
                slab1_vac=slab1.copy();slab2_vac=slab2.copy();
                slab1_vac.center(vacuum=args.vac, axis=2)
                slab2_vac.center(vacuum=args.vac, axis=2)
                atoms=call_castep(slab1_vac,typ="SP",dipolCorr='SC',name='slab1',ENCUT=ecut,KPgrid='1 1 1',PP=pp)
                atoms=slab1_vac
                atoms.set_calculator(calc)
                Eslab1=atoms.get_potential_energy()
                fu1=get_fu(atoms)
                
                atoms=call_castep(slab2_vac,typ="SP",dipolCorr='SC',name='slab2',ENCUT=ecut,KPgrid='1 1 1',PP=pp)
                atoms=slab2_vac
                atoms.set_calculator(calc)
                atoms.calc._label="slab2"
                Eslab2=atoms.get_potential_energy()
                fu2=get_fu(atoms)

        ase.io.write("slab1.cell",slab1.repeat((1,1,1)),format='castep-cell')
        ase.io.write("slab2.cell",slab2.repeat((1,1,1)),format='castep-cell')

        Ws1=(Eslab1-fu1*Ebulk1)/2/surf_area(slab1)/0.01 #A2 to nm2
        Ws2=(Eslab2-fu2*Ebulk2)/2/surf_area(slab2)/0.01 #A2 to nm2
        
        print(('%s: %s eV' % ('Ebulk 1', Ebulk1)))
        print(('%s: %s eV' % ('Ebulk 2', Ebulk2)))
        print(('%s: %s eV' % ('Eslab 1', Eslab1)))
        print(('%s: %s eV' % ('Eslab 2', Eslab2)))
        print(('%s: %.2f eV/nm^2' % ('Wsurf 1', Ws1)))
        print(('%s: %.2f eV/nm^2' % ('Wsurf 2', Ws2)))

        #exit()
        
        #Interface before alignment.
        print("Creating the interface with given slabs.")
        int1=ase.build.stack(slab1, slab2, axis=2, maxstrain=None, distance=dist,cell=None,reorder=True)
        int1.center(vacuum=args.vac,axis=2)
        if args.view:view(int1)
        ase.io.write("interface1.cell",int1.repeat((1,1,1)),format='castep-cell')
        atoms=call_castep(slab1_vac,typ="SP",dipolCorr='SC',name='int1',ENCUT=500,KPgrid='1 1 1',PP=pp)
        #atoms=int1
        #atoms.set_calculator(calc)
        Eint=atoms.get_potential_energy()

        #Wad=(Ws1+Ws2-Eint)/surf_area(int1)/0.01 #A2 to nm2 #check the formula Ea isntead of Wsurf?? This is wrong
        #print('Wad before alingment: %.2f eV/nm^2'%Wad)

        Wad=(Eslab1+Eslab2-Eint)/surf_area(int1)/0.01 #A2 to nm2 #check the formula Ea isntead of Wsurf??
        
        print(('W_ad before alingment: %.2f eV/nm^2\n'%Wad))
        return Eslab1,Eslab2,Ws1,Ws2,int1#surf_area(int1)



def slab_aligner(slab1,slab2,L,Lmax,Lstep,ptol,thickness):
        #find and repeat slabs as specified,
        #One should use only 1 layer for a layer convergence test.
        choice = find_commensurate_supercell(slab1,slab2,Lmax,Lstep,ptol)
        crep = np.ceil(abs(thickness/np.dot(slab1.cell[2],(0,0,1)))) #This does not work properly, doesn't give the thickness desired by the user!!
        #crep=1
        #crep=  np.ceil((args.thickness/get_thickness(slab1)))
        slab1 = cut_cell(slab1,choice[2],choice[3],(0,0,crep))
        slab1 = square_slab(slab1)
        crep = np.ceil(abs(thickness/np.dot(slab2.cell[2],(0,0,1))))
        #crep=1
        #crep=  np.ceil((args.thickness/get_thickness(slab2)))
        slab2 = cut_cell(slab2,choice[4],choice[5],(0,0,crep))
        slab2 = square_slab(slab2)

        #rotate slab2 so that it is alligned with slab1
        ase.build.rotate(slab2,slab1.cell[2],slab2.cell[2],slab2.cell[0],slab1.cell[0])

        #confirm that atom densities are the same as at the start
        atom_density_check(atoms1,slab1)
        atom_density_check(atoms2,slab2)
        
        return slab1,slab2


##################################################
#                  MAIN METHOD                   #
##################################################
               
if __name__== '__main__':
        #read in arguments
        parser = argparse.ArgumentParser()
        parser.add_argument("-i1", "--infile1", default="Li.cell")
        parser.add_argument("-o", "--outfile")
        parser.add_argument("-i2", "--infile2", default="Al.cell")
        parser.add_argument("-m1","--miller1", default=(1,0,0), nargs="+", type=int)
        parser.add_argument("-m2","--miller2", default=(1,0,0), nargs="+", type=int)
        parser.add_argument("-msd","--max_slab_dimension",default=50)
        #parser.add_argument("-th","--thickness",default=7,type=float)
        parser.add_argument("-pt","--percentage_tolerance",default=4)
        parser.add_argument("-conv", "--convLayers",action="store_true",default=False,help="To run a convergence test for slab(s) generated.")

        parser.add_argument("-t","--type",choices=['s','i'],help="s: check different slabs of given bulk structure(s), i: create interface from the slabs with -m1 and -m2 Miller indices from input bulk structures.")

        parser.add_argument("-hubU", "--hubU", type=str,nargs="*",help="For defining the Hubbard U parameters (in eV) for specific atom and orbital types, e.g. -hub Fe d 2.7, Sm f 6.1")

        parser.add_argument("-prog","--prog",choices=['castep','vasp'],default="vasp",help="Code to be used: castep or vasp (def)")

        parser.add_argument("-vac","--vac",default=8,type=float,help='Size of vacuum padding in z direction (only for slabs), def: 8 A on both sites')
        parser.add_argument("-sep","--sep",default=2.0,type=float,help='Seperation between the surface slabs forming the interface, def: 2.5 A')

        parser.add_argument("-th","--thickness",default=7,type=float)
        parser.add_argument("-cr1","--creps1",default=None,type=int)
        parser.add_argument("-cr2","--creps2",default=None,type=int)
        parser.add_argument("-xc", "--xc", default="PBE")

        parser.add_argument("-dry", "--dry",action="store_true",default=False,help="Make a dry run, do not run CASTEP/VASP calculation. Good for checking the slabs and the interfacing algorithm.")

        parser.add_argument("-pass", "--pas",action="store_true",default=False,help="To passivate the top and bottom surface of the interface slab based on dangling bonds.")

        parser.add_argument("-view", "--view",action="store_true",default=False,help="To view the generated structures.")

        args = parser.parse_args()

        #Parse the Hubbard_U info and make the dict
        hubU={}
        if args.hubU:
                U=" ".join(args.hubU)
                xxx=U.split(',')
                for xx in xxx:
                        x=xx.split()
                        if not x[0] in hubU:  hubU[x[0]]={x[1]:x[2]}
                        else: hubU[x[0]][x[1]]=x[2]
                print ("Hubbard U values to be used: ",hubU)

        infile1 = args.infile1
        miller1 = tuple(args.miller1)
        infile2= args.infile2
        miller2 = tuple(args.miller2)

        #2 is going on the bottom so m2 ->-m2 to get right orientation of the surface (in James' new version)
        #miller2 = [-x for x in miller2] #It's better than using the atoms.rotate function, as the latter could change the surface slabs alignment. However this have clashin atoms issue with the ASE stack function.  Bu halen daha tam calismiyor, emin omak lazim !!!


        if args.outfile:
                outfile = args.outfile
        else:
                outfile ="interfaces.out"

        system("date >> %s"%outfile)
        outf=open(outfile,'a')
        

        #tolerances
        Lmax = args.max_slab_dimension 	#max length of cell sides
        L = 5 # initial length
        Lstep = 5
        #T = 1.5	#thickness of slabs
        T=args.thickness
        ptol = args.percentage_tolerance	#percentage tolerances for angle and length matching

        miller_list=[(1,1,1),(1,1,0),(1,0,0),(0,1,0)]

        #Common DFT-related parameters.
        KP=0.10  #KP spacing
        ecut=300 #cutoff energy in eV (convergence prob. with lower cutoffs).
        #pp="00PBE" #pseudopt to use in CASTEP calcs.Def (pp=""): OTF
        pp="OTF"
        dirr="./CASTEP-tmp"

        ecut= 520 #for vasp

        #read in atoms and construct slab, need to repeat atoms to make view work
        print("Reading data from %s and %s."%(infile1,infile2))
        """#Not needed to identify the input file type, ASE is smart enough to detect."""
        if infile1.split(".")[1]=="cell":
                atoms1 = ase.io.read(infile1, format='castep-cell')
        else: #infile1.split(".")[1]=="cif":
                atoms1 = ase.io.read(infile1)
        if infile2.split(".")[1]=="cell": 
                atoms2 = ase.io.read(infile2, format='castep-cell')
        else: # infile2.split(".")[1]=="cif":
                atoms2 = ase.io.read(infile2)            

        initT=time.time()

        print("Structure 1: %s with fu=%d"%(atoms1.get_chemical_formula(),get_fu(atoms1)))
        print("Structure 2: %s with fu=%d"%(atoms2.get_chemical_formula(),get_fu(atoms2)))

        if 0: atoms1=find_prim(atoms1); atoms2=find_prim(atoms2) #Whether to use primitive cells of input structures.

        #Delete calculation files from previous CASTEP run.
        #system("rm -f %s/*"%dirr)
        system("rm -f interface.cell  interface0a.cell  slab1.cell  slab1_aligned.cell  slab1_aligned_opted.cell  slab2.cell  slab2_aligned.cell  slab2_aligned_opted.cell")
        print()


        ########################
        # Do bulk calculations #
        ########################
        #Check if data from a previous run is available (not to repeat the same calcs for bulk).
        fn1="%s/bulk1"%dirr; fn2="%s/bulk2"%dirr
        x1=None;x2=None
        if os.path.exists(fn1+".castep"):
                print ("%s was located, reading data..."%fn1)
                x1=parseCASTEP_ASE(fn1+".geom",atoms1)
                if x1==None:
                        x1=parseCASTEP_ASE(fn1+".castep",atoms1)

        if os.path.exists("bulk1/OUTCAR"):
                print ("%s was located, reading data..."%fn1)
                atoms1=ase.io.read("bulk1/OUTCAR",index=-1)
                x1=[atoms1.get_potential_energy(),atoms1]

        if x1==None:  # whether to compute bulk energies/structures
                print("Computing bulk 1 energy.")
                if args.prog=='castep':x1=call_castep(atoms1,typ="geom",dipolCorr='None',name='bulk1',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU) #normally use K-point spacing.
                elif args.prog=='vasp':x1=call_vasp(atoms1,typ="geom",dipolCorr='None',name='bulk1',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU,slowConv=0,gamma=1,xc=args.xc)

        if os.path.exists(fn2+".castep"):
                print ("bulk1 was located, reading data...")
                x2=parseCASTEP_ASE(fn2+".geom",atoms2)
                if x2==None:
                        x2=parseCASTEP_ASE(fn2+".castep",atoms2)

        if os.path.exists("bulk2/OUTCAR"):
                print ("bulk2 was located, reading data...")
                atoms2=ase.io.read("bulk2/OUTCAR",index=-1)
                x2=[atoms2.get_potential_energy(),atoms2]

        if x2==None:  
                print("Computing bulk 2 energy.")
                if args.prog=='castep':x2=call_castep(atoms2,typ="geom",dipolCorr='None',name='bulk2',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU) #normally use K-point spacing.
                elif args.prog=='vasp':x2=call_vasp(atoms2,typ="geom",dipolCorr='None',name='bulk2',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU,slowConv=0,gamma=1,xc=args.xc) #normally use K-point spacing.


        atoms1=x1[-1]
        Ebulk1=x1[0]

        fu1=get_fu(atoms1)
        sa1=surf_area(atoms1)
        #Ebulk1 /= fu1 #Get the bulk energy per formula unit.
        Ebulk1 /= len(atoms1) #Get the bulk energy per atom.

        atoms2=x2[-1]
        Ebulk2=x2[0]

        fu2=get_fu(atoms2)
        #Ebulk2 /= fu2
        Ebulk2 /= len(atoms2) #Get the bulk energy per atom.
        sa2=surf_area(atoms2)

        #print(x1,x2)
        del x1,x2

        if args.type=="s":
                #str1=""
                min_Ws=1e8; min_mil=();min_slab=None
                for mil in miller_list:
                        slab1 = make_slab(mil,atoms1,repeat=(1,1,1),square=False)
                        slab1.center(vacuum=args.vac, axis=2)
                        if args.prog=='castep':x=call_castep(slab1,typ="sp",dipolCorr='sc',name='slab1',ENCUT=ecut,KPgrid='4 4 1',PP=pp)
                        elif args.prog=='vasp':x=call_vasp(slab1,typ="sp",dipolCorr='sc',name='slab1',ENCUT=ecut,KPgrid='4 4 1',PP=pp,xc=args.xc)

                        slab1=x[-1]
                        Eslab1=x[0]
                        Ws1=(Eslab1-len(slab1)*Ebulk1)/2/surf_area(slab1)/0.01 #A2 to nm2
                        str1='%s (%s_%d%d%d): %.2f eV/nm^2\n' % ('Wsurf 1', slab1.get_chemical_formula(),mil[0],mil[1],mil[2],Ws1)
                        if Ws1<min_Ws: min_Ws=Ws1; min_mil=mil;min_slab=slab1
                        outf.writelines(str1)
                        print (str1)

                if args.convLayers:#Run layer thickness
                        #slab1 = make_slab(min_mil,atoms1,repeat=(1,1,1),square=False)
                        slab1=min_slab
                        nl1,Eslab1,slab1=conv_layers(slab1)#,ifPrim=1)  #Vacuum layer is added automatically within the function.

                        #BURDA!!!
                exit()

        #######################################################################
        #Create the intial slabs with given Miller indices (before alignment).#
        #######################################################################
        print("\nCreating the initial slabs.")
        print (miller1,miller2)
 
        #TODO: the thickness argument apparently has no effect on the created surface slabs, but only on the interface.
        if args.creps1: creps1=args.creps1
        else:creps1=1
        if args.creps2: creps2=args.creps2
        else:creps2=1

        slab1 = make_slab(miller1,atoms1,repeat=(1,1,creps1),square=0) #TODO:replace with creps
        slab2 = make_slab(miller2,atoms2,repeat=(1,1,creps2),square=0)

        #TODO:try ASE cut function instead of make_slab !! But that also uses ase.build.cut in join8.py

        #view(slab2)
        #slab1.rotate(a=180,v='x',rotate_cell=1,center = (0, 0, 0)) #Ithis doesn't have any effect here. #This actually works rather than rotating slab2, nope only rotaes the same interface upside-down.
        #view(slab2)

        """
       #move this to slab aligner.
        if args.creps1: creps1=args.creps1
        else:                 creps1= int((args.thickness/get_thickness(slab1))) #SLAB1 is not defined yet.
        if creps1==0:creps1=1
        if args.creps2: creps2=args.creps2
        else:                 creps2= int((args.thickness/get_thickness(slab2)))
        if creps2==0:creps2=1

        #slab1=slab1.repeat((1,1,creps1));         slab2=slab2.repeat((1,1,creps2))

        slab1 = make_slab(miller1,atoms1,repeat=(1,1,creps1),square=False) #TODO:replace with creps
        slab2 = make_slab(miller2,atoms2,repeat=(1,1,creps2),square=False)
        """

        if 0: #to calculate energies/structures of the initial slabs (before alignment).
                print("Pre-optimizing the initial slabs (before alignment).")
                if 1: #to add vacuum to the slabs (needed) 
                        slab1.center(vacuum=args.vac, axis=2)
                        slab2.center(vacuum=args.vac, axis=2)

                if args.prog=='castep':x=call_castep(slab1,typ="geom",dipolCorr='sc',name='slab1-pre',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=0)
                elif args.prog=='vasp':
                        if args.pas: pas='bot' 
                        else: pas=None
                        x=call_vasp(slab1,typ="geom",dipolCorr='sc',name='slab1-pre',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=0,FixVol=0,xc=args.xc,passivate=pas)

                slab1=x[-1]
                Eslab1=x[0]

                if args.prog=='castep':x=call_castep(slab2,typ="geom",dipolCorr='sc',name='slab2-pre',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=0)
                elif args.prog=='vasp':
                        if args.pas:pas='top'  
                        else: pas=None 
                        x=call_vasp(slab2,typ="geom",dipolCorr='sc',name='slab2-pre',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=0,FixVol=0,xc=args.xc,passivate=pas)

                slab2=x[-1]
                Eslab2=x[0]



        #Not originally here. (Helps increase the overlap between surfaces. i.e. lower lattice misfit).
        niggli=0
        if niggli: niggli_reduce(slab1);niggli_reduce(slab2)

        if 0: slab1=find_prim(slab1);slab2=find_prim(slab2) #does not work.

        if 1: #to add vacuum to the slabs (for demonstration)
                slab1_vac=slab1
                slab2_vac=slab2
                slab1_vac.center(vacuum=args.vac, axis=2) #TODO: add vacuum
                slab2_vac.center(vacuum=args.vac, axis=2)

                if args.prog=="castep":ase.io.write("slab1.cell",slab1_vac.repeat((1,1,1)),format='castep-cell');  ase.io.write("slab2.cell",slab2_vac.repeat((1,1,1)),format='castep-cell')
                elif args.prog=="vasp":ase.io.write("slab1.vasp",slab1_vac.repeat((1,1,1)),format='vasp',vasp5=1);  ase.io.write("slab2.vasp",slab2_vac.repeat((1,1,1)),format='vasp',vasp5=1)

        else:
                if args.prog=="castep":ase.io.write("slab1.cell",slab1.repeat((1,1,1)),format='castep-cell');  ase.io.write("slab2.cell",slab2_vac.repeat((1,1,1)),format='castep-cell')
                elif args.prog=="vasp":ase.io.write("slab1.vasp",slab1.repeat((1,1,1)),format='vasp',vasp5=1);  ase.io.write("slab2.vasp",slab2_vac.repeat((1,1,1)),format='vasp')


        print("\nMisfit (mu) of slabs 1 and 2 (before alignment): %.2f%%"%(misfit(slab1,slab2)*100))#,ifPlot=1)
        print
        #exit()


        ######################
        # Alignment of Slabs #
        ######################
        print("\nAligning the two slabs...")
        slab1,slab2=slab_aligner(slab1,slab2,L,Lmax,Lstep,ptol,T)
        print("\nMisfit (mu) of slabs 1 and 2 (after alignment): %.2f%%"%(misfit(slab1,slab2)*100))#,ifPlot=1)

        if args.view: view(slab1);view(slab2)
        if 0:#Flip the slab2 (top slab) upside-down.still the cell angles cannot be preserved !! using -h,k,-l could solve this. Or can be moved to after optimizations of the slabs (pasivation settings should be changed, from top to bottom)
                slab2.rotate(a=180,v='x',rotate_cell=0,center = (0, 0, 0)) #burda yapinca stack calismiyor.
                slab2.center()
                if args.view: view(slab2)
                a,b,c,alpha,beta,gamma=slab2.get_cell_lengths_and_angles()
                #print (a,b,c,alpha,beta,gamma)
                slab2.rotate(a=alpha/2,v='z',rotate_cell=0,center = (0, 0, 0)) #burda yapinca stack calismiyor.
                slab2.center()

                #slab2.center(vacuum=args.vac,axis=2)
                if args.view:view(slab2)

        if args.prog=="castep":ase.io.write("slab1_aligned.cell",slab1.repeat((1,1,1)),format='castep-cell');  ase.io.write("slab2_aligned.cell",slab2.repeat((1,1,1)),format='castep-cell')
        elif args.prog=="vasp":ase.io.write("slab1_aligned.vasp",slab1.repeat((1,1,1)),format='vasp',vasp5=1);  ase.io.write("slab2_aligned.vasp",slab2.repeat((1,1,1)),format='vasp',vasp5=1)

        if 1:
        #Interface before optimizing the individual slabs.
                if 1: #to delete the vacuum padding in the slabs (needed for correct stacking !!) This should not be done for the no-vacuum calculations !!!
                        slab1.center(vacuum=0, axis=2)
                        slab2.center(vacuum=0, axis=2)
                interface=ase.build.stack(slab1, slab2, axis=2, maxstrain=None, distance=args.sep,cell=None,reorder=1)  #using 0 distance btw slabs gives CASTEP error.
                interface.center(vacuum=args.vac, axis=2) #Vacuum on both sides. For dipole corrections at least 8A vacuum is needed.
                if args.view:view(interface)

                if (interface.get_cell()[2][2]<0): #VASP can't work with upside-down cells (with negative z-coordiantes)
                        print ("Cell is extends towards -z directiom, rotating about x-axis by 180 degree to fix it")
                        interface.rotate(a=180,v='x',rotate_cell=1,center = (0, 0, 0))
                        if args.view:view(interface)

                if args.prog=="castep": ase.io.write("interface0a.cell",interface,format='castep-cell')
                elif args.prog=="vasp": ase.io.write("interface0a.vasp",interface,format='vasp',vasp5=1)
                #ase.io.write("interface0a.cif",interface,format='cif')

        if 1: #to calculate energies/structures of the actual slabs (after alignment).
                print("\nOptimizing the slabs 1 and 2 (after alignment).")

                if 1: #to add vacuum to the slabs (needed)
                        slab1.center(vacuum=args.vac, axis=2)
                        slab2.center(vacuum=args.vac, axis=2)

                if os.path.exists("slab1/OUTCAR"):
                  print ("slab1-aligned was located, reading data...")
                  slab1=ase.io.read("slab1-aligned/OUTCAR",index=-1)
                  x=[slab1.get_potential_energy(),slab1]
                  
                else:
                  if args.prog=='castep':
                        x=call_castep(slab1,typ="geom",dipolCorr='sc',name='slab1-aligned',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=True,hubU=hubU)#,FixList=[1,2])
                        #x=call_castep(slab1,typ="geom",dipolCorr='sc',name='slab1-aligned',ENCUT=ecut,KPspacing=KP,PP=pp,FixCell=1,hubU=hubU)#,FixList=[1,2]) #Optimizer TPSD or FIRE can be used for fixed cell opts.
                  elif args.prog=='vasp': 
                        if args.pas:pas='bot' 
                        else: pas=None
                        x=call_vasp(slab1,typ="geom",dipolCorr='sc',name='slab1-aligned',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=True,hubU=hubU,FixVol=0,xc=args.xc,passivate=pas,nosymm=1)#,FixList=[1,2])

                slab1=x[-1]
                Eslab1=x[0]

                if os.path.exists("slab2/OUTCAR"):
                  print ("slab2-aligned was located, reading data...")
                  slab2=ase.io.read("slab2-aligned/OUTCAR",index=-1)
                  x=[slab2.get_potential_energy(),slab2]
                  
                else:
                  if args.prog=='castep':
                        x=call_castep(slab2,typ="geom",dipolCorr='sc',name='slab2-aligned',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=True,hubU=hubU)#,FixList=[1,2])
                        #x=call_castep(slab2,typ="geom",dipolCorr='sc',name='slab2-aligned',ENCUT=ecut,KPspacing=KP,PP=pp,FixCell=1,hubU=hubU)#,FixList=[1,2]) #Optimizer TPSD or FIRE can be used for fixed cell opts.
                  elif args.prog=='vasp':
                        if args.pas:pas='top'  
                        else: pas=None 
                        x=call_vasp(slab2,typ="geom",dipolCorr='sc',name='slab2-aligned',ENCUT=ecut,KPgrid='1 1 1',PP=pp,FixCell=True,hubU=hubU,FixVol=0,xc=args.xc,passivate=pas,nosymm=1)#,FixList=[1,2])

                slab2=x[-1]
                Eslab2=x[0]


                #Compute the surfafce energies.
                #Ws1=(Eslab1-fu1*Ebulk1)/2/surf_area(slab1)/0.01 #A2 to nm2
                #Ws2=(Eslab2-fu2*Ebulk2)/2/surf_area(slab2)/0.01 #A2 to nm2
                Ws1=(Eslab1-len(slab1)*Ebulk1)/2/surf_area(slab1)/0.01 #A2 to nm2
                Ws2=(Eslab2-len(slab2)*Ebulk2)/2/surf_area(slab2)/0.01 #A2 to nm2

                str1=''
                str1+='%s: %s eV\n' % ('Ebulk 1', Ebulk1)
                str1+='%s: %s eV\n' % ('Ebulk 2', Ebulk2)
                str1+='%s (%s_%d%d%d): %.2f eV/nm^2\n' % ('Wsurf 1', slab1.get_chemical_formula(),miller1[0],miller1[1],miller1[2],Ws1)
                str1+='%s (%s_%d%d%d): %.2f eV/nm^2\n' % ('Wsurf 2', slab2.get_chemical_formula(),miller2[0],miller2[1],miller2[2],Ws2)

                print(str1)
                outf.writelines(str1)

        
        if args.prog=="castep":
                ase.io.write("slab1_aligned_opted.cell",slab1,format='castep-cell')
                ase.io.write("slab2_aligned_opted.cell",slab2,format='castep-cell')
        elif args.prog=='vasp':
                ase.io.write("slab1_aligned_opted.vasp",slab1,format='vasp',vasp5=1)
                ase.io.write("slab2_aligned_opted.vasp",slab2,format='vasp',vasp5=1)
 
        if args.view:view(slab1);view(slab2)

        if args.vac>2: #to delete the vacuum padding in the slabs (needed for correct stacking !!)
                slab1.center(vacuum=0, axis=2)
                slab2.center(vacuum=0, axis=2)
        #Create the interface.
        interface=ase.build.stack(slab1, slab2, axis=2, maxstrain=None, distance=args.sep,cell=None,reorder=True)  #using 0 distance btw slabs gives CASTEP error.
        #interface = wrap_coords(interface)
        #interface.wrap()
        interface.center(vacuum=args.vac, axis=2) #Vacuum on both sides. For dipole corrections at least 8A vacuum is needed.

        #print (interface.get_cell()[2]) 
        #a,b,c,alpha,beta,gamma=interface.get_cell_lengths_and_angles()
        #print (a,b,c,alpha,beta,gamma)
        if args.view:view(interface)

        if (interface.get_cell()[2][2]<0): #VASP can't work with upside-down cells (with negative z-coordiantes)
                print ("Cell is extends towards -z directiom, rotating about x-axis by 180 degree to fix it")
                interface.rotate(a=180,v='x',rotate_cell=1,center = (0, 0, 0))
                if args.view:view(interface)

        #ase.io.write("interface.cell",interface,format='castep-cell')

        if args.prog=="castep": ase.io.write("interface.cell",interface,format='castep-cell')
        elif args.prog=="vasp": ase.io.write("interface.vasp",interface,format='vasp',vasp5=1)

        #ase.io.write("interface.cif",interface,format='cif')

        niggli=1
        if niggli: print ("Niggli reduce interface...");niggli_reduce(interface)
        #if niggli: interface=niggli_reduce(interface)
        if args.view:view(interface)

        if 1:
                print("\n Single point run for the interface geometry.")
                if args.prog=="castep":
                        #x=call_castep(interface,typ="sp",dipolCorr='SC',name='interface',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU)
                        x=call_castep(interface,typ="sp",dipolCorr='SC',name='interface-pre',ENCUT=ecut,PP=pp,KPgrid='1 1 1',hubU=hubU)                
                elif args.prog=="vasp": x=call_vasp(interface,typ="sp",dipolCorr='SC',name='interface-pre',ENCUT=ecut,PP=pp,KPgrid='1 1 1',hubU=hubU,FixVol=0,xc=args.xc,nosymm=1)

                interface=x[-1]
                Eint=x[0]


                Wad=(Eslab1+Eslab2-Eint)/surf_area(interface)/0.01 #A2 to nm2 #check the formula Ea isntead of Wsurf??
                str1='W_ad before alingment: %.2f eV/nm^2\n'%Wad
                print(str1)
                outf.writelines(str1)

        if 1:
                print("\nOptimizing the final interface geometry.")
                if args.prog=="castep": x=call_castep(interface,typ="geom",dipolCorr='SC',name='interface',ENCUT=ecut,PP=pp,KPspacing=KP,hubU=hubU)
                elif args.prog=="vasp": x=call_vasp(interface,typ="geom",dipolCorr='SC',name='interface',ENCUT=ecut,PP=pp,KPgrid='1 1 1',hubU=hubU,FixVol=1,xc=args.xc,nosymm=1) #fix vol works for some cases (YIG) but not for some (Al2O3//Si)

                interface=x[-1]
                Eint=x[0]


                Wad=(Eslab1+Eslab2-Eint)/surf_area(interface)/0.01 #A2 to nm2 #check the formula Ea isntead of Wsurf??
                str1='W_ad after alingment: %.2f eV/nm^2\n'%Wad
                print(str1)
                outf.writelines(str1)

                if args.prog=='castep': ase.io.write("interface_opted.cell",interface,format='castep-cell')
                elif args.prog=='vasp': ase.io.write("interface_opted.vasp",interface,format='vasp',vasp5=1)
                ase.io.write("interface_opted.cif",interface,format='cif')
                      

        outf.close()

        print ("Elapsed time: %.2f sec."%( time.time()-initT))

        exit()


##################
# Excluded Parts #
##################


