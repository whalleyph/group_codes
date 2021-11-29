#!/usr/bin/env python
from __future__ import print_function
import argparse,re
from sys import argv,exit,version_info
import numpy as np
from os import system,popen#,popen4 #as popen  #popen4 is not available in Python3
from subprocess import Popen,PIPE # as popen # check_output
from time import time
import os.path
from ase import Atoms
import ase.io
import numpy.linalg as la
from copy import deepcopy as dc

from spglib import find_primitive,standardize_cell,get_spacegroup #, niggli_reduce #niggli_reduce from ASE-tools conflicts with that from spglib.
from ase.spacegroup import Spacegroup,get_spacegroup

#from scipy.stats import linregress
#from timeit import timeit
#from pandas import rolling_mean
#import pandas as pd


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def cum_mov_avg(x):
    #print(x,len(x))
    arr=[]
    CMA=x[0]#0.
    for ind,val in enumerate(x):
        #ind+=1
        CMA = (CMA*ind+val)/(ind+1)
        arr.append(CMA)
    return arr

def Popen4(cmd):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate()
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err


def grep(key,fname,n=-1): 
    #Uses grep to get the nth line with the keyword from a given file. By default the last occurence is returned.
    try:
        #return popen4('grep -m %d "%s" %s '%(n,key,fname),"r")[1].readlines()[-1][0:-1] #don't take \n at the end.  #Fastest one!!! #Not available in python3
        return Popen4('grep -m %d "%s" %s 2>/dev/null'%(n,key,fname))[0]
    except:
        return ""

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def sort_dict_by_val(data,index=0):
    #index: index to use if values of the dictionary is not scalar.
    keys=list(data.keys())
    if len(keys)==1: return keys

    skeys=sorted(keys)#sorted keys
    for i in range(len(keys)):
        key=keys[i]
        val=data[key]
        try: val=val[index]
        except:None
        flag=0
        for j in range(len(skeys)):
            key2=skeys[j]
            val2=data[key2]
            try: val2=val2[index]
            except:None
            if val <= val2 and i!=0: skeys.insert(j,key);flag=1;break
        if not flag and key not in skeys:skeys.append(key)#not to miss the largest element.

    return skeys

def Euclid(*integers): #for finding the gretest common denominator (gcd).
    if type(integers[0])==list or type(integers[0])==tuple:
        integers = integers[0];return sum(integers)
    if len(integers) == 1:return integers
    elif len(integers) == 2:
        a,b = integers[0],integers[1]
        while b != 0:
            a,b = b, a%b
        return a
    else:
        return Euclid(integers[0],Euclid(integers[1:]))
    
def getAtoms(filen):#Reads the CASTEP cell files to get the atomic content.
    try:inpf=open(filen)
    except: print ("getAtoms: %s could not be found.");return ""
    #lines=inpf.readlines()
    atoms={}
    flag=False
    key1=""; key2=""
    ext=filen.split("/")[-1].split(".")[1]
    if ext=="cell": key1="%block positions_frac";key2="%endblock positions_frac"
    elif ext=="castep": key1="x----";key2="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    
    for ln in inpf.readlines():
        if not flag and key1 in ln.lower():  flag=True;continue
        if flag and key2 in ln.lower() :flag=0;break
        if flag:
            #print ln
            ln=ln.replace("x","")
            at=ln.split()[0]
            if not at in atoms:atoms[at]=1
            else: atoms[at]+=1
        
    keys=atoms.keys()
    #keys.sort()
    keys=sorted(keys)
    content=""
    for key in keys:
        content+=key
        if atoms[key]!=1:content+=str(int(atoms[key]))
        
    inpf.close()

    fu=getfu(atoms)

    reduced,red_dict=getChemForm(atoms)
    return content,fu,reduced,red_dict

def getChemForm(atoms): #can be integrated back to above for efficiency.
    #Get the reduced chemcial formula of the system.
    keys=atoms.keys()
    #keys.sort()
    keys=sorted(keys)
    fu=getfu(atoms)
    reduced=""
    red_dict={}
    for key in keys:
        red_dict[key]=atoms[key]/fu
        #print atoms[key]
        #reduced[key]=int(atoms[key])/fu
        reduced += key
        if atoms[key]/fu!=1:reduced+=str(int(atoms[key]/fu))

    return reduced,red_dict

def calcEf_orig(atoms,cpots,H):
    #Calculate the formation energy per atom usign the chem. potentials provided.
    if len(cpots)==0: print ("Error in calcEf, empty chempots dictionary !!!"); return 0.0

    Ef=H
    total=sum(atoms.values())

    for key in atoms.keys():
        Ef -= atoms[key]*cpots[key]
    return Ef/total #Ef per atom
 
def calcEf(atoms,cpots,H):
    #Calculate the formation energy per atom usign the chem. potentials provided.
    if len(cpots)==0: print ("Error in calcEf, empty chempots dictionary !!!"); return 0.0

    #input atoms are taken in the reduced form (i.e. per fu)

    Ef=H  #is this H/atom or total H ??? correct: H / fu
    total=sum(atoms.values())
    #fu=getfu(atoms) #should be fu=1 

    akeys=list(atoms.keys())
    pkeys=list(cpots.keys())

    if set(akeys).issubset(set(pkeys)):# the case monoatomic references (i.e. constituent atoms)
        for key in akeys:
            Ef -= atoms[key]*cpots[key]  #In this case cpots should be in E/fu form.

        return Ef/total #Ef per atom  (and per fu) (both are equal for monoatomic systems.

    else: #the case of binary/unary species as the reference.
        #ADD HERE THE ADDITIONAL PRODUCTS as USER INPUT.
        reduced,red_dict=getChemForm(atoms)#this only works with a list of atoms.
        a=[];b=[]
        rkeys=red_dict.keys()
        #rkeys.sort() #reduced_dict of atoms.
        rkeys=sorted(rkeys)

        #print "Refs:",pkeys,rkeys 
        for key in rkeys:#over atom types
            tmp=[]
            b.append(red_dict[key])
            for pkey in pkeys: #over ref species.
                cnt=0
                at=Atoms(pkey)#use ASE built-in function to get the atom content (symbols) from the chem. formula. #TODO: write an own function using RE.                
                for sy in at.get_chemical_symbols():
                    if sy == key: cnt+=1
                tmp.append(cnt)
            a.append(tmp)

        a=np.array(a);b=np.array(b)
        #Check if there is an all-zero entry of coeffs. #??CORRECT??
        for ax in a:
            if sum(ax)==0:
                #print "Ef (%s) can't be computed  with given ref. chem. potentials."%(reduced)
                return 999.999999
        try:
            m=la.lstsq(a,b)[0]
            #print m
        except:
            #print "Ef (%s) can't be computed  with given ref. chem. potentials."%(reduced)
            return 999.999999

        for i in range(len(pkeys)): #over ref species.
            pkey=pkeys[i]
            Ef -= m[i]*cpots[pkey]#*len(at)   

        return Ef/total

def getReduced(atoms):
    #Determine the formula unit and reduced formula.
    if isinstance(atoms,str):
        atoms=Atoms(atoms)
    elif isinstance(atoms,ase.Atoms):
        atoms=atoms

    ats={}
    aTypes=atoms.get_chemical_symbols()
    for at in aTypes:
        if at in ats: ats[at]+=1
        else: ats[at]=1

    #ats=re.findall(r'([A-Z][a-z]*)(\d*)', stoich) #gives [('H', '2'), ('S', ''), ('O', '4')]

    try:
        fu=1
        vals=list(ats.values())
        for i in range(2,min(vals)+1):
            fl=1
            for j in vals:
                if j%i!=0:fl=0
            if fl:fu=i
    except: print ("Error in fu determination, fu set to 1");   fu=1

    reduced=""
    keys=sorted(ats.keys())
    for key in keys:
        reduced += key
        if ats[key]/fu!=1:reduced+=str(int(ats[key]/fu)) 

    return reduced,fu  
    

def getfu(atoms):
    #Determine the formula unit.
    #print Euclid(vals)
    try:
        fu=1
        vals=atoms.values()
        for i in range(2,min(vals)+1):
            fl=1
            for j in vals:
                if j%i!=0:fl=0
            if fl:fu=i
        #print fu
    except: print ("BLA  Error in fu determination, fu set to 1");   fu=1
    return fu

def readAirssOut(inpf): #To read the airss.out (i.e. stdout created by the AIRSS code).
    spg="";sym="";tspin=0.0;mult="";fname='';cnt=0;H=0.0
    conv=[] #converged attempts
    lowestE={};data={}
    fus={} #holds the formula unit info for each composition.

    for ln in open(inpf).readlines():
        if re.search("ENFORCED SYMMETRY",ln):sym=(" ").join(ln.split()[2:])
        elif re.search("Space group ",ln):spg=ln.split()[2:]
        elif re.search("Multiplicities: ",ln):mult=ln.split()[1].replace("[","")
        elif re.search ("'finalTotalSpin':",ln):
            try:tspin=float(ln.split()[1][0:-1])
            except:
                try:tspin=float(ln.split()[3][0:-1])
                except:tspin=0.0
        elif re.search("seed",ln) and re.search("_final'",ln) and re.search("converged",ln):#Final step for the current attempt.

            cnt+=1 #no Attempts (converged or non-converged.)
            fname=ln.split("'seed':")[1].split("'")[1]+".cell"
            ifConv=ln.split("'converged': ")[1].split(",")[0]

            #print ifConv
            if ifConv=="True": # converged geometry opt.
                conv.append(fname.split('.')[0])#fname w/out extension.
                #print conv[-1]
                x=getAtoms(fname)
                atoms=x[0]
                fu=x[1]
                H=float(ln.split("'enthalpy': ")[1].split(",")[-1].replace("]","").replace("[","").replace("{","").replace("}","")) #take the last item on the enthlpy list.

                #take H from .castep file, as airss.out can be mistaken.

                try:pres=float(popen('grep "Pressure:" %s| tail -n 1 '%(conv[-1]+".castep")).read()[0:-1].split(":")[1].split()[0])
                except:pres=0.0

                try:vol=float(popen('grep "Current cell volume =" %s| tail -n 1 '%(conv[-1]+".castep")).read()[0:-1].split("=")[1].replace("A**3",""))
                except:vol=0.0

                if 1: #Use SG info from AIRRS output (i.e. intended/designated space group)
                    key=(", ").join([spg[0][0:-1],spg[1],sym])
                else: #Use the actual space group read from CASTEP output.
                    spg_new=popen('grep "Space group" %s'%(conv[-1]+".castep")).read()[0:-1].replace(" ","").split(",")[0].split("=")[1].split(":")
                    key=(", ").join(["#"+spg_new[0],spg_new[1]])


                dt=[conv[-1],mult,tspin,H,pres,vol]

                if not atoms in fus:fus[atoms]=fu
                if not atoms in data:data[atoms]={key:[dt]};
                else:
                    if not key in data[atoms]:data[atoms][key]=[dt]
                    else: data[atoms][key].append(dt)

                #Update the lowest energy info for each space group.
                if not atoms in lowestE:lowestE[atoms]={key:[H,dt]}
                else:
                    if not key in lowestE[atoms]:lowestE[atoms][key]=[H,dt]
                    else:
                        if H < lowestE[atoms][key][0]:lowestE[atoms][key]=[H,dt]

    return [data,lowestE,conv,cnt,fus]

def readASEfiles(fileList,chempots={},ifSep=False):

    dt=[fname,H,mult,tSpin,pres,vol,fu,nIons,Ef]
    return [data,lowestE,conv,cnt]#,fus,noAtoms,reds]

    
def readInpfs(logList,chempots={},ifSep=False):
    conv=[] #converged attempts
    cnt = len(logList)
    lowestE={};data={}

    for log in logList:
        if log=="":continue
        #if grep("Geometry optimization completed successfully.",log,1)=="": continue #Not converged geometry, skipping.   #Now  filtering (searching for *_final.castep) is done at the very beginning of the main code (that is clearly faster), yet it might not filter out the geometries with non-converged SCF/geometry (depending on how the filter is set below). The second check does bring in only minor additional cost (<10%).
        
            
        tSpin=0.0; tCharge=0.0; vol=0.0;pres=0.0;mult="";H=0.0;E=0.0;nIons=0
        ext=log.split('.')[-1]
        if ext=='castep':
            for ln in open(log):#.readlines(): #readlines is outdated and slow.
                if re.search("net spin",ln):tSpin=float(ln.split(":")[-1])
                elif re.search("net charge",ln):tCharge=float(ln.split(":")[-1])
                elif re.search("Pressure:",ln):pres=float(ln.split(":")[1].split()[0])
                elif re.search("Current cell volume =",ln):vol=float(ln.split("=")[1].replace("A**3",""))
                elif re.search("Space group",ln):spg_new=ln.replace(" ","").split(",")[0].split("=")[1].split(":")
                elif re.search("Final free energy \(E-TS\)    =",ln):H=float(ln.split("=")[1].split()[0])
                elif re.search("Final energy, E",ln): E=float(ln.split("=")[1].split()[0])
                elif re.search("Total number of ions in cell =",ln):nIons=int(ln.split("=")[1])
                #elif re.search("Final bulk modulus =",ln):bulkMod=int(ln.split("=")[1])

                #elif re.search("",ln):
        else:
            ASEatoms=ase.io.read(log)
            try:tSpin=ASEatoms.get_magnetic_moment()
            except:tSpin=0.
            tCharge=0.
            try:pres=ASEatoms.info['pressure']
            except:pres=-999.
            vol=ASEatoms.get_volume()
            #SG=str(get_spacegroup(ASEatoms, symprec=args.tol).symbol.replace(' ',''))
            #print(SG)           

            #ASEatoms.info['spacegroup']=SG.replace('(','').replace(')','')
            SG=get_spacegroup(ASEatoms, symprec=args.tol)
            #spg_new=SG.replace('(','').replace(')','')
            #spg_new="%d, %s"%(SG.no,SG.symbol.replace(' ',''))
            spg_new=[SG.no,SG.symbol.replace(' ','')]
            #print(spg_new)
            ASEatoms.info['spacegroup']=spg_new
            H=ASEatoms.info['energy']; E=H
            nIons=len(ASEatoms)
            


        if H==0.0 and E==0.0: continue #Do not include logfiles with scf error or other errors.
        conv.append(log)

        #key=(", ").join(["#"+spg_new[0],spg_new[1]])
        key="#%s, %s"%(spg_new[0],spg_new[1])
        #print('mykey:',key)

        if ext=='castep': atoms,fu,reduced,atoms_dict=getAtoms(log)
        else: 
            atoms=ASEatoms.get_chemical_formula(mode='hill', empirical=False)
            reduced=ASEatoms.get_chemical_formula(mode='hill', empirical=True)
            reduced2,fu=getReduced(ASEatoms)
            
            
            #print(atoms,reduced,reduced2,fu)

            #Get reduced number of each atom type as dict
            ats={}
            aTypes=ASEatoms.get_chemical_symbols()
            for at in aTypes:
                if at in ats: ats[at]+=1
                else: ats[at]=1

            atoms_dict={}
            for at in ats.keys():
                atoms_dict[at]=ats[at]#/fu
        #print(atoms,reduced,ats,atoms_dict,fu,aTypes)

        H /= fu
        E /= fu

        if len(chempots)!=0:
            #Decide whether to use E or H.
            Ef=calcEf(atoms_dict,chempots,H)#reduced atomlist is input. H must be used for compatibility with the rest of the output types, and also used in MATADOR.
            H=Ef
        else:Ef=0.0

        if not ifSep:  atoms=reduced
        
        fname=log.replace(".%s"%ext,"")
        dt=[fname,H,mult,tSpin,pres,vol,fu,nIons,Ef]


        if not atoms in data:data[atoms]={key:[dt]};
        else:
            if not key in data[atoms]:data[atoms][key]=[dt]
            else: 
                flag=0
                for j in range(len(data[atoms][key])):
                    d=data[atoms][key][j]
                    if H < d[1]: #this is for addition of new entry in ascending enthlapy order. 
                        flag=1
                        data[atoms][key].insert(j,dt)
                        break
                if not flag: data[atoms][key].append(dt) #not to miss the entry with highest enthalpy.


        #lowestE stores the lowest enthalpy info for each space group. (This has become reduntant with energy-ordered data).
    for atoms in data.keys():
        lowestE[atoms]={}
        for key in data[atoms].keys():
            #print atoms,key,data[atoms][key]
            lowestE[atoms][key]=[data[atoms][key][0][1],data[atoms][key][0]] #add the lowest entalphy element. i.e. [H,dt]


    return [data,lowestE,conv,cnt]#,fus,noAtoms,reds]
 

def boltz_dist(energies,T=298.15,omega=[]):#Return the occupation probabilities of the configurations at a given temperature based on their energies.
    kb= 8.6173303*10**-5 #Boltzmann constant (eV/K).
    if len(omega)==0:#If the degeneracies are not given explciitly.
        omega=[1 for E in energies]

    if 1: #Get relative energies.  Doesn't really matter as long as you use the normalised factors.
        mn=min(energies)
        energies=[i-mn for i in energies]
    probs=[]
    for E in energies:
        probs.append(math.exp(-E/kb/T))
    #Normalise    
    Z=sum(probs) #i.e. partition fnc
    probs=[Pn/Z for Pn in probs]

    #Configurational statistics as given in R. Grau-Crespo et al. J.Phys: COndens. Matter, 19,2007,256201
    E_avg=sum([energies[i]*probs[i] for i in range(len(energies))])
    print ("\nAverage energy of the sytem in configurational equilibirum,  E=%.5f eV"%E_avg)

    F=-kb*T*np.log(Z)
    print ("Configurational free energy in the complete space, F=%.5f eV"%F)

    S= (E_avg-F)/T
    print ("Configurational entropy in the complete space, S=%.5f eV/K"%S)

    Smax=kb*np.log(len(energies))
    print ("Upper limit of config. entropy, Smax= %.5f eV/K"%Smax)

    #Now count in the degenaricies of the configs.
    Sm=[kb*T*np.log(om) for om in omega] #degeneracy entropy

    #for i,E in enumerate(energies):
    Em_bar=[energies[i]-T*Sm[i] for i in range(len(energies))]

    Pm=[np.exp(-Em_bar[i]/kb/T) for i in range(len(energies))] 
    Z_bar=sum(Pm)
    #Pm_bar=[(1/Z)*np.exp(-Em_bar[i]/kb/T) for i in range(len(energies))] #reduced  probability for an independent config.
    Pm_bar=[P/Z_bar for P in Pm]

    E_avg=sum([Em_bar[i]*Pm_bar[i] for i in range(len(energies))])

    F=-kb*T*np.log(Z_bar)
    print ("Configurational free energy in the reduced config. space, F=%.5f eV"%F)

    S= (E_avg-F)/T
    print ("Configurational entropy in the reduced config. space, S=%.5f eV/K"%S)

    return Pm_bar#,E_avg

#############
#Main script#
#############
parser = argparse.ArgumentParser()
parser.add_argument("-i","--inpf", type=str, nargs="*", default=None, #default="airss.out",
                    help="Input file to read the data from.")

parser.add_argument("-pa","--pa", action='store_true', default=False,
                    help="Report Enthalpies and volumes per atom.")
parser.add_argument("-s","--summarise", action='store_true', default=False,
                    help="Report only the ground state entry for each stoichometry randomly.")
#parser.add_argument("-ss","--sortsumm", action='store_true', default=False,
#                    help="Report only the ground state entry for each stoichometry in ascending order of energy. (Slower).")

parser.add_argument("-all","--allCastep", action='store_true', default=False,
                    help="Use all *.castep files in cwd and subdirectories. Def: Only *_final.castep used.")
                    
parser.add_argument("-sp","--singlepoint", action='store_true', default=False,
                    help="Consider the single point type calculations in CASTEP files. Def: Only converged structures from geometry optimsation are used.")
                    
parser.add_argument("-t","--top", type=int,default=5,
                    help="How many top enetries to show for each stoichometry.Def: 5")

parser.add_argument("-d","--details", action='store_true', default=False,
                    help="Reports detailed information for upto -t/--top entries for each sapce group. Default: false")

parser.add_argument("-hideSG","--hsg", type=int,default=0,
                    help="Hide entries with symmetry lower than the given space group. Default: All entries shown.")

parser.add_argument("-sep","--separate", action='store_true', default=False,
                    help="To separate the different fu entries of the same composition (e.g. X2Y4 besides XY2. By default all entries with same composition but diff fu s will be grouped.")

parser.add_argument("-ef","--compEf", action='store_true', default=False,
                    help="Report formation energy per atom  instead of enthalpy per fu/atom. This option requires chemical potential of the constitutents atoms, see -mu keyword for details.")

parser.add_argument("-mu","--chempots", type=str, nargs='*',
                    help="Chemical potentials of the constitutent atoms for computing the formation energy per atom for comparing different stoichometries.Energies should be per fu. \nUsage: -mu Li=value Si=value ...")

parser.add_argument("-aref","--atomicRef", action='store_true', default=False,
                    help="Convert the input polyatomic reference energies into chemical potential of the atomic constituents. Def: False, so binary/ternary species could be used as chem. pot. ref. See also -mu tag.")

parser.add_argument("-ma","--molar", type=str,
                    help="Run a partial molar property analysis for a chosen atom type. e.g. -ma Li Def: False ")

parser.add_argument("-save","--save", action='store_true', default=False,
                    help="Save the geometries for the top T number of geometries with lowest enthalpies (no mater which space group), where T can be set by -t/--top flag. Default: No geometry saved. Geometry files will be stored in separate folders for each composition under SAVED-GEOM folder.")

#parser.add_argument("-save2","--save_v2", action='store_true', default=False,
#                    help="Save the geometries with lowest energy from each of top T Space Group, where T can be set by -t/--top flag. Default: No geometry saved. Geometry files will be stored in separate folders for each composition under SAVED-GEOM folder.")

parser.add_argument("-ow","--overwrite", action='store_true', default=False,
                    help="Overwrite the folder where saved geometries are kept. Def: No overwriting. Caution: enabling this will delete the previous directory.")

parser.add_argument("-ec","--ecut", type=float,
                    help="The upper limit of formation energy for the structures to be saved. Def: All geometries are saved.")

parser.add_argument("-plot","--plot", action='store_true', default=False,
                    help="Plots the 2D/3D phase diagram using ASE tools. Requires -ef and -mu keywords. Def: False ")

parser.add_argument("-dim","--dim", type=int,choices=[2,3],default=None,
                    help="Force 2D or 3D plot of the phase diagram (by ASE). Default: The default in ASE is used.")

parser.add_argument("-decomp","--decomp", action='store_true', default=False,
                    help="Decomposes the binaries/ternaries/quaterneries into its constitutents using ASE tools. Requires -ef and -mu keywords. Def: False ")

parser.add_argument('-tol', '--tol',type=float, default=1e-3,help="The symmetry tolerance. Def: 1e-3")

args = parser.parse_args()

args.save_v2=args.save #Delete this when fixed the save2 flag.

#Necessary assignments.
if args.molar:args.summarise=True
if args.summarise:args.sortsumm=True


initT=time()
#data,lowestE,conv,cnt,fus=readAirssOut(args.inpf)

sfold="SAVED-GEOM"
if args.save or args.save_v2:
    if os.path.exists(sfold):
        if args.overwrite:
            system('rm -rf ./%s; mkdir -p ./%s'%(sfold,sfold))
        else: print ("%s dir already exists, terminating..."%sfold) ; exit()
    else:
        system(' mkdir -p ./%s'%(sfold))

if args.ecut:
    if not args.compEf:
        print ("Energy cutoff -ec is only compatible with formation energies, so -ef must be switched on. Terminating...");exit()
    elif not args.save:
        print ("Energy cut off (-ec) only works for geometry saving, not for filtering out the displayed geometries !!!")

#Get the chem. pot. for comp. Ef.
cpots={};cpots_orig={}
if args.compEf: 
    if not args.chempots: print ("For computing formation energy (Ef) per atom, you need to provide the chemical potential of the consitutent atoms. Use -mu option to provide the potentials."); exit()
    else:
        for pot in args.chempots: 
            x=pot.split("=")
            try:
                if x[0] in cpots: cpots[x[0]].append(float(x[1]))
                else:cpots[x[0]]=[float(x[1])]
            except: print ("Error in decoding "),pot;exit()

        for pot in cpots.keys():  cpots[pot]=np.mean(cpots[pot])

if args.atomicRef: #convert  multiatomic-species reference energies into individual atomic chemical potentials. It does not affect the already monoatomic references. #essentially gives the same result as the case without -aref (since the same procedure is followed).
    print ("Extracting atomic chemical potentials from %d input reference species: %s"%(len(cpots),", ".join(cpots.keys())))

    gatIDs=[] #general list of atom symbols covered by the reference species.
    for pot in cpots.keys():
        #pot=cpots[pot]
        at=Atoms(pot)#use ASE built-in function to get the atom content (symbols) from the chem. formula. #TODO: write an own function using RE.                
        for sy in at.get_chemical_symbols():
            if sy not in gatIDs: gatIDs.append(sy)
    #print (gatIDs)

    b=[];a=[] #Ax=b; results array
    for pot in cpots.keys(): #over ref species.
        val=cpots[pot]
        b.append(float(val))
        catIDs={}

        at=Atoms(pot)#use ASE built-in function to get the atom content (symbols) from the chem. formula. #TODO: write an own function using RE.                
        for sy in at.get_chemical_symbols():
            if sy in catIDs: catIDs[sy]+=1
            else: catIDs[sy]=1
        a.append([catIDs[ID] if ID in catIDs else 0 for ID in gatIDs ])

    a=np.array(a);b=np.array(b)
    m=la.lstsq(a,b)[0]
    #create new ref. cpot dict.
    cpots_orig=dc(cpots);del cpots;cpots={}
    for i in range(len(gatIDs)):
        cpots[gatIDs[i]]=m[i]
    print ("List of derived atomic chemical potentials: \n", cpots)
    print()

if args.plot:  
    from ase.phasediagram import PhaseDiagram
    refs=[(x,0.0) for x in cpots.keys()]
    if len(cpots_orig)>0: refs.extend([(x,0.0) for x in cpots_orig.keys()])
    terns=[]#ternaries

if args.inpf: logList=args.inpf #if input files given by the user explicitly.
else: #by default read CASTEP log files, i.e. *.castep files
    if args.singlepoint: 
        #logList=[x[0:-1] for x in popen4('grep -l "Geometry optimization completed successfully." `find ./ -type f -name "*.castep"`','r')[1].readlines()] #get rid of \n at each line. 
        logList=Popen4('grep -l "single point energy" `find ./ -type f -name "*.castep"`')[0]
        if len(logList)==0:print ("No *.castep was found. Terminating...");exit()

    elif args.allCastep: 
        #logList=[x[0:-1] for x in popen4('grep -l "Geometry optimization completed successfully." `find ./ -type f -name "*.castep"`','r')[1].readlines()] #get rid of \n at each line. 
        logList=Popen4('grep -l "Geometry optimization completed successfully." `find ./ -type f -name "*.castep"`')[0]
        if len(logList)==0:print ("No *.castep was found. Terminating...");exit()

    else: #check only *_final.castep (faster)
        logList=Popen4('grep -l "Geometry optimization completed successfully." `find ./ -type f -name "*_final.castep"`')[0]
        if len(logList)==0:
            print ("No *_final.castep was found. Switching to *.castep...")
            logList=Popen4('grep -l "Geometry optimization completed successfully." `find ./ -type f -name "*.castep"`')[0]

print ('%d files to analyse.'%len(logList))
data,lowestE,conv,cnt=readInpfs(logList,cpots,args.separate)


keys1=data.keys()
#keys1.sort()
keys1=sorted_nicely(keys1)#sort keys in alphanumerical order in a proper way.

print ("\n\nTotal of %d attempts were found, %d of which have a converged geometry/SCF (for Geometry opt/Single Points). \nThey can be grouped into %d different stoichiometries formula: "%(cnt,len(conv),len(keys1)),(", ").join(keys1))


if args.compEf: lbl="Ef / atom"
elif args.pa: lbl='H / atom'
else: lbl='H / fu'

sort_summ={} #summary sorted wrt enthalpies.

if args.summarise:
    print ("\n%-22s %-13s %-8s %-12s %-4s %-6s %-8s %-10s %-24s"%('Stoichiometry','Space Group','#Entries',lbl,'fu','TSpin','Pressure','Volume/fu','Name of Lowest-H Entry'))
    print (22*'-',13*'-',8*'-',12*'-',4*'-',6*'-',8*'-',10*'-',24*'-' )


for key in keys1: #compositions
    keys2=data[key].keys() #space groups.

    #Sort in ascending Enthalpy.
    vals={}
    for k in lowestE[key].keys():#space groups.
        noAtoms=lowestE[key][k][1][7]
        fu=lowestE[key][k][1][6]
        if args.compEf:    vals["%.12f"%(lowestE[key][k][0])]=k 
        elif args.pa:     vals["%.12f"%(lowestE[key][k][0]/noAtoms*fu)]=k 
        else:           vals["%.12f"%(lowestE[key][k][0])]=k
        
    nk=[float(x) for x in vals.keys()]
    #nk.sort()
    nk=["%.12f"%k for k in sorted(nk)]
    new_keys=[vals[k] for k in nk]

    if args.save:
        system('mkdir -p ./%s/%s'%(sfold,key))

    """
    if args.save_v2: #This copies the top entries from each of the first args.tol SG entries. Does not work.
        system('mkdir -p ./%s/%s'%(sfold,key))
        #print key
        cnt=1
        #tmp_H
        for ky in new_keys: #space groups
            fn=lowestE[key][ky][1][0].replace(".cell","") #filename
            system('cp %s.res %s.cell ./%s/%s/ 2>/dev/null '%(fn,fn,sfold,key))
            cnt+=1
            if cnt > args.top:break
            """

        
    if args.summarise:
        min_key=vals[str(nk[0])]
        k=min_key
        noAtoms=lowestE[key][k][1][7]
        fu=lowestE[key][k][1][6]
        if args.compEf: None
        elif args.pa: lowestE[key][k][0] /= noAtoms/fu
        else:None #lowestE[key][k][0] /= fu
        noEntry=0 # no of str for each stoichometry.
        for j in data[key].keys(): noEntry+= len(data[key][j])

        if args.sortsumm:
            sort_summ[key]=[lowestE[key][k][0],key,k,noEntry,lowestE[key][k][0],fu,lowestE[key][k][1][3],lowestE[key][k][1][4],lowestE[key][k][1][5]/fu,lowestE[key][k][1][0] ]
        else:
            print ("%-22s %-15s %-6d %-13.6f %-3d %-5.2f %8.4f %8.3f   %-25s"%(key,k,noEntry,lowestE[key][k][0],fu,lowestE[key][k][1][3],lowestE[key][k][1][4],lowestE[key][k][1][5]/fu,lowestE[key][k][1][0] ))

        if args.save:
            #print lowestE[key][k][0]
            if args.ecut == None or float(lowestE[key][k][0]) <= args.ecut: 
                fn=lowestE[key][k][1][0].replace(".cell","") #filename
                system('cp %s.res %s.cell ./%s/%s/ 2>/dev/null '%(fn,fn,sfold,key))

        if args.plot:None
            #refs.append((key,lowestE[key][k][0]))        

        continue  #this needs to be here!!

    print ("\n Stoichiometry: %s"%(key))

    print ("%-13s %-8s %-12s %-4s %-10s %-8s %-10s %-24s"%('Space Group','#Entries',lbl,'fu','TotalSpin','Pressure','Volume/fu','Name of Lowest-H Entry'))
    print (13*'-',8*'-',12*'-',4*'-',10*'-',8*'-',10*'-',24*'-' )

    cnt2=1
    for k in new_keys: #space groups
        #print (k)
        noAtoms=lowestE[key][k][1][7]
        fu=lowestE[key][k][1][6]
        #print k,fu,noAtoms,lowestE[key]
        if cnt2 > args.top: break
        if int(k.split(',')[0].split("#")[1]) <= args.hsg: continue
        if args.compEf: None
        elif args.pa: lowestE[key][k][0] /= noAtoms/fu
        else:None

        if args.details:
            for j in range(0,args.top):
                try: # scaling not done for j==0, Scaling Not needed for j==0 as using lowestE.
                    fu=data[key][k][j][6]
                    noAtoms=data[key][k][j][7]
                    if args.compEf: None
                    elif args.pa: data[key][k][j][1] /= noAtoms/fu
                    else: None
                    if j==0:
                        print ("%-15s %-6s %-13.6f %-5d %-7.2f %8.4f %8.3f   %-25s"%(k,len(data[key][k]),data[key][k][j][1],fu,data[key][k][j][3],data[key][k][j][4],data[key][k][j][5]/fu,data[key][k][j][0]))
                    else: 
                        print ("%-15s %-6s %-13.6f %-5d %-7.2f %8.4f %8.3f   %-25s"%("","",data[key][k][j][1],fu,data[key][k][j][3],data[key][k][j][4],data[key][k][j][5]/fu,data[key][k][j][0]))            
                except: break

                if args.plot: refs.append((key,data[key][k][j][1]))

                if args.save:
                    if args.ecut == None or (float(data[key][k][j][1]) <= args.ecut):
                        fn=data[key][k][j][0].replace(".cell","") #filename
                        system('cp %s.res %s.cell ./%s/%s/ 2>/dev/null '%(fn,fn,sfold,key))
                
        else:
            print ("%-15s %-6d %-13.6f %-5d %-7.2f %8.4f %8.3f   %-25s"%(k,len(data[key][k]),lowestE[key][k][0],fu,lowestE[key][k][1][3],lowestE[key][k][1][4],lowestE[key][k][1][5]/fu,lowestE[key][k][1][0] ))
            if args.plot: refs.append((key,lowestE[key][k][0]))

            if args.save:
                if args.ecut == None or float(data[key][k][0][1]) <= args.ecut:
                    fn=data[key][k][0][0].replace(".cell","") #filename 
                    system('cp %s.res %s.cell ./%s/%s/ 2>/dev/null '%(fn,fn,sfold,key))
        cnt2+=1

"""
#xb=boltz_dist(energies,args.boltz)
    if args.boltz != 0: #If Boltzmann temp is non-zero.
        #This is already done in the function.
        mn=min(energies)
        energies=[i-mn for i in energies]
        xb=boltz_dist(energies,args.boltz)
        if args.EperAtom:        print "\nConfig ID, Rel. Energy/atom [eV], Boltzmann Weight (at %.2f K)"%args.boltz
        else:     print "\nConfig ID, Rel. Energy [eV], Boltzmann Weight (at %.2f K)"%args.boltz
        print "%s"%"-"*60
        for i in range(len(energies)):
            print "%-15s %6.3f  %.4f"%(  args.inpf[i].replace(".magres","").split("./")[-1],energies[i],xb[i])
    else:
        xb=[1.0/len(energies) for i in range(len(energies))] #If no Boltzmann distribution then all configs are weighted equally.
"""

if args.summarise and args.sortsumm:
    keys=sort_dict_by_val(sort_summ,index=0)

    for key in keys:
        print ("%-22s %-15s %-6d %-13.6f %-3d %-5.2f %8.4f %8.3f   %-25s"%(sort_summ[key][1],sort_summ[key][2],sort_summ[key][3],sort_summ[key][4],sort_summ[key][5],sort_summ[key][6],sort_summ[key][7],sort_summ[key][8],sort_summ[key][9]))
        
        if args.plot: refs.append((key,sort_summ[key][4]))

#molar=1
if args.molar:#Report the (partial) molar properties against the mole fraction of a slected atom type.
    atype=args.molar
    print ("\n Writing partial molar properties for %s atoms in molar.dat ..."%atype)
    outf=open("molar.dat",'w')
    outf2=open("molar_str.dat",'w')
    outf.write("#x_%-3s  %3s_x   %s (eV)   mu_%-3s(eV)  V_total (A^3) V_molar (A^3)  V_molar2 (A^3)  Bulk Modulus (GPa) Stoichiometry \n"%(atype,atype,lbl,atype))
    outf2.write("#x_%-3s  %3s_x  a, b, c, alpha, beta, gamma   total V  \n"%(atype,atype))

    ext=logList[0].split('.')[-1]
    #keys=sorted_nicely(sort_summ.keys())#sort keys in alphanumerical order in a proper way.
    keys=sort_dict_by_val(sort_summ,index=0)
    dE=[]; dN=[] ; data=[]
    for ind,key in enumerate( keys):
        fu=sort_summ[key][5]
        #atoms=Atoms(key,cell=[1,1,1]).repeat((fu,1,1))#.sort()
        fname=sort_summ[key][9]
        #atoms=ase.io.read(sort_summ[key][9]+".castep",format='castep-castep')
        atoms=ase.io.read(fname+"."+ext)

        nX=float(len([atom.index for atom in atoms if atom.symbol==atype]))
        syms=list(atoms.get_chemical_symbols())
        #print(syms)
        #gcd=Euclid([syms.count(sym) for  sym in np.unique(syms) if sym!=atype]) #greatest common
        #gcd=np.gcd.reduce([syms.count(sym) for  sym in np.unique(syms) if sym!=atype]) #greatest common
        #gcd=Euclid(np.unique(syms,return_counts=1)[1])
        gcd=np.gcd.reduce(np.unique(syms,return_counts=1)[1])
        #print(fname,fu,gcd)
        if gcd==0:gcd=1

        ntot=float(len(atoms))
        vol_fu=sort_summ[key][8]
        
        if args.compEf: E=sort_summ[key][4] #Ef/atom
        elif args.pa: E=sort_summ[key][4]#*ntot #H/atom
        else: E=sort_summ[key][4]*fu  #H/fu, sosal total H
        mu_X=0.#E/nX
        try:bulkMod=float(grep("Final bulk modulus =",sort_summ[key][9]+".castep")[-1].split()[-2]) #GPa #take the last entry only
        except:bulkMod=0. #TODO: take from ASEatoms object?? Not available

        dE.append(E);dN.append(nX)

        #Get a,b,c,alpha,beta,gamma
        if 1:
            a,b,c,alpha,beta,gamma=atoms.get_cell_lengths_and_angles()
        else:
            a,b,c,alpha,beta,gamma=0,0,0,0,0,0

        data.append([nX/ntot,nX/gcd,E,mu_X,vol_fu*fu,vol_fu*fu*nX/ntot,bulkMod,atoms.get_chemical_formula(mode='hill'),E,nX,a,b,c,alpha,beta,gamma])
        #outf.write("%5.3f %7.3f %13.6f %13.6f %8.3f %8.3f %8.3f %8.3f  #%s\n"%()

    data=sorted(data,key=lambda x: x[0]) #sort by the molar fraction of X

    #mu_X=np.gradient(np.array(dE),np.array(dN))

    dE=[];dN=[];dV=[]
    for dt in data:
        dE.append(dt[8]); dN.append(dt[9]);dV.append(dt[4])

    mtype=['none','scale_by_nX','rm','cma'] #averaging/mean type: none, cumulative moving avg., running/rolling mean
    mtype=mtype[0]

    #np.gradient compute slope using the diff wrt to the first element. (not a progressive one).
    #mu_X=np.gradient(np.array(dE),np.array(dN),edge_order=1) #chem pot of X
    #V_bar=np.gradient(np.array(dV),np.array(dN)) #partial molar volume of X

    mu_X=(np.ediff1d(dE)/np.ediff1d(dN)).tolist()
    V_bar=(np.ediff1d(dV)/np.ediff1d(dN)).tolist()
    #mu_X.insert(0,0.); V_bar.insert(0,0.)

    if mtype=='scale_by_nX':
        mu_X=(np.ediff1d(dE)/np.ediff1d(dN)/dN[1:]).tolist() #mu_X per X atom.
        V_bar=(np.ediff1d(dV)/np.ediff1d(dN)/dN[1:]).tolist()
        mu_X.insert(0,0.); V_bar.insert(0,0.)

    elif mtype=="rm":
        N=2
        if 1:
            for i in range(N):mu_X.insert(0,0.); V_bar.insert(0,0.)
            mu_X=running_mean(mu_X,N)
            V_bar=running_mean(V_bar,N)
        else: #use pandas implementation (gives the same results).
            import pandas as pd
            mu_X.insert(0,0.); V_bar.insert(0,0.)
            mu_X=pd.Series(mu_X)
            mu_X=mu_X.rolling(N).mean()
            V_bar=pd.Series(V_bar)
            V_bar=V_bar.rolling(N).mean()

    elif mtype=='cma':
        #mu_X.insert(0,0.); V_bar.insert(0,0.)
        cummean = lambda x:  x.cumsum()/np.arange(1, len(x)+1) #does the same as cum_mov_avg
        #mu_X=cum_mov_avg(mu_X)
        #V_bar=cum_mov_avg(V_bar)
        mu_X=cummean(np.array(mu_X)).tolist()
        V_bar=cummean(np.array(V_bar)).tolist()
        mu_X.insert(0,0.); V_bar.insert(0,0.)

    else: #no averaging just the actual derivatives.
        mu_X.insert(0,0.); V_bar.insert(0,0.)


    #print (np.array(dE),np.array(dN))
    #print (mu_X, len(mu_X))
    #print (np.array(dV),np.array(dN))
    #print (V_bar, len(V_bar))

    for ind,mu in enumerate(mu_X):
        data[ind][3]=mu
        data[ind][5]=V_bar[ind]

    data=sorted(data,key=lambda x: x[0]) #sort by the molar fraction of X

    #Write output
    for ind,key in enumerate( keys):
        #outf.write("{d[0]:5.3f} {d[1]:7.3f} {d[2]:13.6f} {0:13.6f} {d[4]:8.3f} {d[5]:8.3f} {d[6]:8.3f} {d[7]:8.3f}  #{d[8]:s} \n".format(mu_X[ind],d=data[ind]))
        outf.write("{d[0]:5.3f} {d[1]:7.3f} {d[2]:13.6f} {d[3]:13.6f} {d[4]:8.3f} {d[5]:8.3f} {d[6]:8.3f}  #{d[7]:s} \n".format(d=data[ind]))
        outf2.write("{d[0]:5.3f} {d[1]:7.3f} {d[10]:8.3f} {d[11]:8.3f} {d[12]:8.3f} {d[13]:8.3f} {d[14]:8.3f} {d[15]:8.3f} {d[4]:8.3f} #{d[7]:s} \n".format(d=data[ind]))

    outf.close()

print ("Elapsed time: %.2f sec."%( time()-initT))

if args.ecut:ecut=args.ecut
else: ecut=999.
if args.decomp:
    pd = PhaseDiagram([x  for x in refs if x[1] < ecut])
    for i in [x for x in refs if x[1] < ecut]:
        try: print ("%s with Ef = %.2f eV decomposes into: "%(i[0], i[1])); pd.decompose(i[0])
        except:  print ("Error in decompostion of %s"%i[0])

if args.plot:
    print
    #print refs
    print ("Filtering out the entries with high Ef...")
    sel=[x  for x in refs if x[1] < ecut or x[1]==0.0]
    #if ecut<0: sel.extend([])
    pd = PhaseDiagram(sel)

    pd.plot(ax=None,dims=args.dim,show=True)#,only_label_simplices=False, only_plot_simplices=False) #ax:dimension for projection.



#now copy all the lowest energy structures into subfolders of compositions and space groups for following CASTEp calculations.


exit()


##########
# BACKUP #
##########
#logList=[x[0:-1] for x in popen4('find ./ -type f -name "*.castep"','r')[1].readlines()] #get rid of \n at each line.   


"""
    if args.save: #This copies the top args.tol number of entries, no mater what space group they belong to.
        system('mkdir -p ./%s/%s'%(sfold,key))
        #print key
        cnt=1
        all_H={}
        for ky in new_keys: #space groups
            #fn=data[key][ky][0]
            all_H
        for ky in     
            fn=lowestE[key][ky][1][0].replace(".cell","") #filename
            system('cp %s.res %s.cell ./%s/%s/ 2>/dev/null '%(fn,fn,sfold,key))
            cnt+=1
            if cnt > args.top:break
"""

""" #This approach might be faster but not compatible with new Popen due to too long cmd arguments.
    #files=[x[0:-1] for x in popen4('find ./ -type f -name "*.castep"','r')[1].readlines()]
    #files=popen4('find ./ -type f -name "*.castep"','r')[1].readlines()
    files=Popen4('find ./ -type f -name "*.castep"')[0]
    print (files)
    if len(files)==0:print ("No *.castep was found. Terminating...");exit()
    #logList=[x[0:-1] for x in popen4('grep -l "Geometry optimization completed successfully." `printf "%s"`'%''.join(files),'r')[1].readlines()] #get rid of \n at each line. 
    logList=Popen4('grep -l "Geometry optimization completed successfully." `printf "%s"`'%''.join(files))[0]
"""

"""
    #files=popen4('find ./ -type f -name "*_final.castep"','r')[1].readlines()
    files=Popen4('find ./ -type f -name "*_final.castep"')#[0]
    if len(files)==0:print ("No *_final.castep was found. Switching to *.castep...");files=Popen4('find ./ -type f -name "*.castep"')[0]
    #logList=[x[0:-1] for x in popen4('grep -l "Geometry optimization completed successfully." `printf "%s"`'%''.join(files),'r')[1].readlines()] #get rid of \n at each line.
    logList=Popen4('grep -l "Geometry optimization completed successfully." `printf "%s"`'%''.join(files))[0]
"""

"""
    #new_keys=[vals[str("%.12f"%k)] for k in nk]
    new_keys=[]
    for k in nk:
        if "%.12f"%k in vals:  		new_keys.append(vals[str("%.12f"%k)])
        else: print("key=%s not found !!"%k)
    """
