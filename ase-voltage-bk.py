#!/usr/bin/env python3
#Written by Bora Karasulu, requires ASE.
from __future__ import print_function

import ase.io,time
import argparse,os.path,re
from ase import Atoms
from ase.build import sort
#import ase.build
from copy import deepcopy as dc
import numpy as np
from os import environ,system
from copy import deepcopy
from ase import Atoms,Atom
from ase.neighborlist import neighbor_list
import matplotlib
matplotlib.use('TkAgg') #needed as OVITO uses the default QT5 backend
import matplotlib.pyplot as plt
from sys import exit,stdout,argv,version_info

from ase.neighborlist import NeighborList
#from time import time
import time
import numpy.linalg as la
from subprocess import Popen,PIPE # as popen # check_output


from scipy.constants import physical_constants

FARADAY_CONSTANT_Cpermol = physical_constants['Faraday constant'][0]
HARTREE_TO_EV = physical_constants['Hartree energy in eV'][0]
Cperg_to_mAhperg = 2.778e-1
C_TO_mAh = Cperg_to_mAhperg
BOHR_TO_ANGSTROM = physical_constants['Bohr radius'][0] * 1e10
RY_TO_EV = physical_constants['Rydberg constant times hc in eV'][0]
KBAR_TO_GPA = 0.1
eV_PER_ANGSTROM_CUBED_TO_GPa = 160.21776
AVOGADROS_NUMBER = physical_constants['Avogadro constant'][0]
ANGSTROM_CUBED_TO_CENTIMETRE_CUBED = 1e-24
ELECTRON_CHARGE = physical_constants['elementary charge'][0]
KELVIN_TO_EV = physical_constants['kelvin-electron volt relationship'][0]

#From MATADOR
def get_binary_grav_capacities(x, m_B):
    """ Returns capacity in mAh/g from x/y in A_x B_y
    and m_B in a.m.u.
    """
    x = np.array(x)
    if m_B != 0:
        return x * FARADAY_CONSTANT_Cpermol * Cperg_to_mAhperg / m_B
    return float('NaN')

def get_generic_grav_capacity(concs, elements):
    """ Returns gravimetric capacity of
    <elements[0]> in mAh/g of matador doc.
    """
    tmp_concs = np.array(concs, copy=True)
    # if no Li, capacity = 0...
    # tmp_concs /= np.min(concs)
    x = tmp_concs[0]
    if x == 0:
        return 0.0
    masses = dict()
    m_B = 0
    for elem in elements:
        masses[elem] = get_molar_mass(elem)
    for ind, elem in enumerate(elements):
        if ind == 0:
            continue
        else:
            m_B += masses[elem] * tmp_concs[ind]
    Q = get_binary_grav_capacities(x, m_B)
    return Q

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
        return Popen4('grep -m %d "%s" %s '%(n,key,fname))[0][0]
    except:
        return ""

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


def calcEf(atoms,cpots,H):
    #Calculate the formation energy per atom usign the chem. potentials provided.
    if len(cpots)==0: print ("Error in calcEf, empty chempots dictionary !!!"); return 0.0

    #input atoms are taken in the reduced form (i.e. per fu)

    Ef=H  #is this H/atom or total H ??? correct: H / fu
    total=sum(atoms.values())
    #fu=getfu(atoms) #should be fu=1 

    akeys=atoms.keys()
    pkeys=cpots.keys()

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

def Euclid(*integers): #for finding the gretest common denominator (gcd).
    if type(integers[0])==list or type(integers[0])==tuple:
        integers = integers[0]
    if len(integers) == 2:
        a,b = integers[0],integers[1]
        while b != 0:
            a,b = b, a%b
        return a
    else:
        return Euclid(integers[0],Euclid(integers[1:]))
    
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
    except: print ("Error in fu determination, fu set to 1");   fu=1
    return fu

def sorted_version( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def sort_dict_by_val(data,index=0):
    #index: index to use if values of the dictionary is not scalar.
    keys=list(data.keys()) # python3 requires a list construction.
    if len(keys)==1: return keys

    skeys=[keys[0]]#sorted keys  
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
            #print(skeys)#,data[skeys[0]])
        bot=data[skeys[0]][index]
        top=data[skeys[-1]][index]
        if not flag and key not in skeys:#not to miss the last element.
            if val>=top: skeys.append(key)
            elif val<=bot: skeys.prepend(key)
    return skeys

parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')

parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')

parser.add_argument("-mu","--chempots", type=str, nargs='*',
                    help="Chemical potentials of the constitutent atoms for computing the formation energy per atom for comparing different stoichometries.Energies should be per fu. \nUsage: -mu Li=value Si=value ...")

parser.add_argument("-aref","--atomicRef", action='store_true', default=False,
                    help="Convert the input polyatomic reference energies into chemical potential of the atomic constituents. Def: False, so binary/ternary species could be used as chem. pot. ref. See also -mu tag.")

parser.add_argument("-at","--atype", type=str,default='Li',
                    help="Run a partial molar property analysis for a chosen atom typeAtomic type to be used in the voltage calcuations. e.g. -at Li  ")

parser.add_argument("-c","--comp", type=str,required=1,
                    help="Chemical composition to be analysed. No default")

parser.add_argument("-ef","--compEf", action='store_true', default=False,
                    help="Report formation energy per atom  instead of enthalpy per fu/atom. This option requires chemical potential of the constitutents atoms, see -mu keyword for details.")

parser.add_argument("-ec","--ecut", type=float,
                    help="The upper limit of formation energy for the structures to be saved. Def: All geometries are saved.")

parser.add_argument("-plot","--plot", action='store_true', default=False,
                    help="Plots the 2D/3D phase diagram using ASE tools. Requires -ef and -mu keywords. Def: False ")

parser.add_argument("-dim","--dim", type=int,choices=[2,3],default=None,
                    help="Force 2D or 3D plot of the phase diagram (by ASE). Default: The default in ASE is used.")

parser.add_argument("-decomp","--decomp", action='store_true', default=False,
                    help="Decomposes the binaries/ternaries/quaterneries into its constitutents using ASE tools. Requires -ef and -mu keywords. Def: False ")

parser.add_argument("-vol","--vol", action='store_true', default=False,
                    help="Also plot the volume change vs. capacity graph. Def: False ")


args = parser.parse_args()
initT=time.time()


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

    print(cpots)
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

elems=[]
for inpf in args.inpf:
    #print (inpf)
    atoms=ase.io.read(inpf)
    elems.extend([at for at in atoms.get_chemical_symbols()])

elems=np.unique(elems)
#print(elems)
try:cp=[cpots[el] for el in elems]
except: cp=None

from glob import glob
from matador.hull import QueryConvexHull
from matador.scrapers.castep_scrapers import res2dict
#%matplotlib inline
cursor = [res2dict(inpf)[0] for inpf in args.inpf]
#print(cursor)
#hull = QueryConvexHull(cursor=cursor, elements=elems)

hull = QueryConvexHull(cursor=cursor, no_plot=0, voltage=1, volume=args.vol,summary=1,chempots=cp,composition=args.comp,labels=1,hull_cutoff=0.5) #species=elems) 
#, hull_cutoff=7.5e-2,kpoint_tolerance=0.03) #subcmd='volume' , query=None , chempots=['','']
#hull = QueryConvexHull(cursor=cursor, no_plot=0, subcmd='volume',summary=1)#,species=elems) 

#hull.volume_curve()#prints the cell volumes
#plot_volume_curve()
#plotting.plot_voltage_curve(hull)

#print('Filtering down to only ternary phases... {}'.format(len(hull.hull_cursor)))
#hull.hull_cursor = [doc for doc in hull.hull_cursor if len(doc['stoichiometry']) == 3]
#print('Filtering unique structures... {}'.format(len(hull.hull_cursor)))
#uniq_list, _, _, _ = list(get_uniq_cursor(hull.hull_cursor[1:-1], debug=False))
#cursor = [hull.hull_cursor[1:-1][ind] for ind in uniq_list]
#print('Final cursor length... {}'.format(len(cursor)))
#print('over {} stoichiometries...'.format(len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in cursor]))))
#print([doc['stoichiometry'] for doc in cursor])

#polished_hull = QueryConvexHull(cursor=cursor, subcmd='hull', no_plot=1, composition=['Li3P:Li2S'],summary=True, spin=0,details=False,labels=False,source=False,loose=False,intersection=True,hull_cutoff=0.05)  #composition=['LiPS'], 'Li2S:P2S5', 'Li2S:Li3P'


exit()

print ("Elapsed time: %.2f sec."%( time.time()-initT))


