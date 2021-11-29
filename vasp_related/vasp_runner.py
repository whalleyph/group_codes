#!/usr/bin/env python3
import ase.io
from ase.visualize import view
from os import system,popen,chdir,getcwd,getenv,putenv,listdir,environ
import os.path,glob
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
import argparse

import os.path
from ase.phonons import Phonons



from subprocess import Popen,PIPE # as popen # check_output
from sys import stdout,version_info
def Popen4(cmd,timeout=None):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate(timeout=timeout)
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err


def vasp_continue(ifDel=0):
        fl = listdir(".")
        fl.sort()
        #print fl
        flag=0
        ind=0
        for fn in fl:
                if fn[0:3]=="RUN" :
                        ind = int(fn[3:])

        run="RUN%d"%(ind+1)
        print ("Previous data stored in %s."%run)
        system('mkdir %s'%run)
        system('cp * %s/ 2> /dev/null '%run)
        system('cp CONTCAR POSCAR') #this was missing in the previous version !!
        if ifDel: system('rm WAVECAR CHGCAR')


def grep(key,fname,n=-1): 
    #Uses grep to get the nth line with the keyword from a given file. By default the last occurence is returned.
        try:
                return popen4('grep -m %d "%s" %s '%(n,key,fname),"r")[1].readlines()[-1][0:-1] #don't take \n at the end.  #Fastest one!!!
        except: 
                return popen('grep -m %d "%s" %s '%(n,key,fname),"r").read()[0:-1]
        else:
                return ""



def call_vasp(atoms,calc=None, typ="sp",xc='PBE',name='./VASP-tmp',PP='',resDir="",dipolCorr=False,dipolDir='z',KPgrid="1 1 1",KPspacing="", ifPrint=False, ENCUT=300,ifManualRun=True,FixCell=False,FixVol=False,FixList=[],hubU={},ispin=0,restart=None,nelm=150,nsw=800,etol=1E-8,ediffg=-0.01,slowConv=0,ismear=0,sigma=0.01,passivate=None,gamma=1,algo="Fast",nosymm=0):
        import ase.io

        #exe=popen('echo "$VASP_COMMAND"',"r").read()[0:-1]
        exe=getenv('VASP_COMMAND')
        if not exe: exe='vasp_std'
        print("VASP command: ", exe)

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
                calc = Vasp(atoms,restart=restart, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt='vasp.out') 
                print ("Previous run was read successfully, checking if converged...")
        except: #if restarting fails
                print("Problem with restarting from the previous run. Starting a new calculation.")
                calc = Vasp(atoms,restart=None, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt='vasp.out');restart=False

        if restart: #if restarted no need to assign the variables and use the ones read from INCAR and OUTCAR by ASE.
                if Vasp.read_convergence(calc):#calc.converged:#read_convergence('OUTCAR'):
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


        if 1 and len(hubU)!=0: 
            #TODO: USE ASE implementation instead!!
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
                else:  # bothways is on by default in this. This works well.
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

def call_vasp_v2(fname,exe=None):

    if not exe: exe=getenv('VASP_COMMAND') ; #print ('bk')
    if not exe: exe='vasp_std'
    
    os.environ['VASP_COMMAND']=exe
    print("VASP command: ", exe)#, os.environ['VASP_COMMAND'])

    cwd=os.getcwd()

    seed=fname.split('.')[0]
    try:chdir(seed)
    except:None
    
    system('vasp_make_potcar -xc potpaw_PBE.54 &> /dev/null')

    #try:
    #       atoms = ase.io.read('vasprun.xml',index=0)
    #       print "vasprun.xml was read successfully."
    #except:
    try:
        try:calc = Vasp(restart=True)
        except:calc=None
        atoms = calc.get_atoms()
        print ("VASP run was read succesfully from OUTCAR.")
        if Vasp.read_convergence(calc): print('Geom opt already converged...')
        else:print('Geom opt not converged running a continuation job...')

        #else:
    except:
        print ("VASP run could not be read, starting a new run...")
        calc=Vasp()
        #calc = Vasp2(atoms,restart=restart, directory="./", label='vasp', ignore_bad_restart_file=False, command=exe, txt=None)
        calc.read_incar(filename='INCAR')
        if os.path.exists("./OUTCAR"):   vasp_continue()
        atoms=ase.io.read("POSCAR",format="vasp")
        #try:atoms=ase.io.read("POSCAR",type="vasp",format="vasp4)
        calc.set(xc="pbe",ibrion=2)
        calc.directory="."#cdir
        atoms.set_calculator(calc)


    #try: 
    Ep = atoms.get_potential_energy()
    
    #except:Ep=float(popen("grep 'free  energy   TOTEN' OUTCAR | head -n 1").readlines()[0][0:-1].split()[-1]);#print Ep  
    #TODO:use popen4 instead

    chdir(cwd)

    return Ep,atoms


parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                    help='an integer for the accumulator')


parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')
#parser.add_argument('-ot','--otype', type=str,required=True, help='Output file type')
#parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: input name is used as root.')
parser.add_argument('-it','--itype', type=str,required=False, help='Input file type. Def: determined automatically from the extension.')

parser.add_argument('-t', '--tol',type=float, default=1e-4,help="The symmetry tolerance.")


parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')

parser.add_argument('-dry','--dryrun', default=False,action='store_true', help='Dry run: do not call VASP, just gather information from VASP run folders. Def: No')

parser.add_argument('-np', '--nprocs',type=int, default=32,help="No of processes to start for each VASP calculation throuh mpirun. Def:32")

parser.add_argument('-exe','--exe', type=str,required=False, default='vasp_std',help='Vasp exectuable. Def: vasp_std')

parser.add_argument('-ph','--phonons', default=False,action='store_true',help='Do a phonon calculation if the geometry is converged.')
parser.add_argument('-ph_sc','--phonons_supercell', type=int,nargs=3,default=(2,2,2),help='Size of the supercell to be uused in the phonon calculation.Def: No supercell used, i.e. 1 1 1')
parser.add_argument('-ph_del', '--phonons_delta',type=float, default=0.05,help="Stepsize to use in the finite difference method to compute force constants. Def:0.05 A")

parser.add_argument('-ph_np','--phonons_noplot', default=False,action='store_true',help='Do a phonon calculation if the geometry is converged.')

args = parser.parse_args()

initT=time.time()

from ase.spacegroup import Spacegroup,get_spacegroup


exe="mpirun -np %d %s"%(args.nprocs,args.exe)

if 0:
  files=Popen4("""grep "reached required accuracy - stopping structural energy minimisation" C-*/vasp.out -L | awk -F/ '{print $1}' """) [0]

  print('List of calculations with non-converged geometry: ', files)

  cwd=os.getcwd()
  for f in files:
    #os.chdir(f)
    #if os.path.exists("%s/%s-final.res"%(f,f)): system('cp %s/%s-final.res %s.res'%(f,f,f))
    #else:
      atoms=ase.io.read('%s/XDATCAR'%f, index=-1)
      ase.io.write('%s.res'%f,atoms,format='res')

    #os.chdir(cwd)
  #exit()

system ('mkdir -p ./inputs')
cwd=os.getcwd()
cnt=0
import random
#for inpf in random.choices(args.inpf):
for inpf in args.inpf:
    print (inpf)

    if not os.path.isfile(inpf):
        print('%s does not exists, skipping...'%inpf)
        continue

    atoms=ase.io.read(inpf)
    seed=inpf.split('.')[0]
    system('mkdir -p %s'%seed) 
    atoms.write('%s/POSCAR'%seed,format='vasp',vasp5=1)
    system('cp -f INCAR  POTCAR KPOINTS %s/'%(seed))
    system('cp %s inputs/; mv %s %s/'%(inpf,inpf,seed))

    chdir(seed)

    #if Popen4(""" grep 'reached required accuracy - stopping structural energy minimisation' """)[0]: 
    if not args.dryrun:   
        try:Ep,atoms=call_vasp_v2(inpf,exe)
        except:print('Problem with %s, skipping...'%inpf);chdir(cwd);continue
    else: 
      try: Ep=atoms.info['energy']; 
      except:Ep=0.0; #atoms=ase.io.read(inpf)
      if Ep==None: Ep=0.0


    #atoms.write('%s/%s-final.res'%(seed,seed),format='res')
    atoms.info['energy']=Ep #atoms.get_total_energy()
    #atoms.calc.results['energy']=Ep #atoms.get_total_energy()
    #atoms.calc.set('energy')=Ep #atoms.get_total_energy()
    #print(dir(atoms.calc))
    #print(dir(atoms))
    #print (atoms.energy)
    atoms.info['name']=inpf

    #ASE does not get the Pressure right when restarting
    #P=atoms.info.get('pressure')
    #if P is None: 
    P=float(Popen4("""grep pressure OUTCAR | tail -1 | awk '{print ($4+$9)*0.1}' """)[0][0]) #kB to GPa
    atoms.info['pressure']=P

#grep "free  energy   TOTEN  = " $seed/OUTCAR  | tail -1 | awk '{print " VASP: Final Enthalpy     = "$5}' >> $seed.castep

    #if P is not None: atoms.info['pressure']=P
    #else: atoms.info['pressure']=0.
    SG=atoms.info.get('spacegroup')
    if SG is None: 
        try:SG=Popen4('symm -cl %s'%seed,timeout=30)[0][0]  #part of the airss package #now only for cluster, also support the crystal mode.
        except:SG=''
        print (inpf,SG)#,get_spacegroup(atoms, symprec=args.tol))

    if SG is None or SG == '':   SG=str(get_spacegroup(atoms, symprec=args.tol).symbol.replace(' ',''))
    atoms.info['spacegroup']=SG.replace('(','').replace(')','')
    atoms.info['times_found']=1
    atoms.info['name']=seed
    #print(atoms.info.keys())

    #print(atoms.info['spacegroup'])
    #print(atoms.info['energy'])
    
    print('%s: Ep=%.5f P= %.2f GPa SG=%s'%(inpf,Ep,P,atoms.info['spacegroup']))


    myres=ase.io.res.Res(atoms)
    myres.energy=Ep
    myres.write_file('%s.res'%(seed), write_info=0,  significant_figures=6) #works!

    #other ways
    #ase.io.res.write_res('%s.res'%(seed), atoms, write_info=0, write_results=1, significant_figures=6) #worksbut energy is written as None even though the value is stored.    
    #ase.io.write('%s-final.res'%(seed), atoms,format='res')


    #Do the phonons calculations (if requested by the user and if the geometry has converged)
    if args.phonons: 
        print('Doing phonons calculations as requested')
        calc=atoms.get_calculator()
        if not Vasp.read_convergence(calc): print('Geometry has not converged, skipping the phonons calculation...')
        else:
            
            #try:
            #    from ase.dft.kpoints import get_special_points
            #    path = get_special_points(atoms.cell,eps=args.tol)#.special_path.split(',')[0]

            #    print ('High symmetry path was determined for the input structure as: ',path)
            #except: raise Error; #continue

            if os.path.exists('phonon'): 
                print('Previous phonons calc was found in the %s/phonon directory, will read existing force data and do the missing points.'%seed)
                for f in glob.glob("phonon/*.json",recursive=0):
                    if os.path.getsize(f)==0: system('rm %s'%f); print (f,' is empty, deleted')#delete the empty files to recalculate

            calc.directory=calc.directory+"/Phonon-tmp/"
            print(calc.directory)
            calc.set(nsw=0,ibrion=-1,algo='fast') #do a single point at each phonon step
            calc.set(lwave=1) #save wavecar for restarting electronic wavefunction

            #N=1 #TODO: get this from user
            #ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
            ph = Phonons(atoms, calc, supercell=args.phonons_supercell, delta=args.phonons_delta)
            ph.run()

            # Read forces and assemble the dynamical matrix
            method='Frederiksen' #'standard'
            #ry:
            ph.read(acoustic=True,method=method)
            #xcept:ph.read(acoustic=True,method=method)
            #ph.clean()

            path = atoms.cell.bandpath(None, npoints=100) #if special_path is set to None,then it will be automatically determined.

            bs = ph.get_band_structure(path)
            dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)
            print(path)

            if args.phonons_noplot: continue

            # Plot the band structure and DOS:
            fig = plt.figure(1, figsize=(7, 4))
            ax = fig.add_axes([.12, .07, .67, .85])

            #emax = 0.035 
            emax=max(dos.get_energies()) #TODO: get this from user
            bs.plot(ax=ax, emin=0.0, emax=emax)

            dosax = fig.add_axes([.8, .07, .17, .85])
            dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color='grey',
                               edgecolor='k', lw=1)

            dosax.set_ylim(0, emax)
            dosax.set_yticks([])
            dosax.set_xticks([])
            dosax.set_xlabel("DOS", fontsize=18)

            fig.savefig('%s_phonon.png'%seed)
            #plt.show()



    chdir(cwd)
