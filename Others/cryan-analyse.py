#!/bin/env python3

#Required packages: pip install soprano sklearn dscribe ovito seaborn p_tqdm --user  #also pyqt5

from subprocess import Popen,PIPE,DEVNULL,run,CompletedProcess,call # as popen # check_output
import shlex
from sys import exit,stdout,argv,version_info
import matplotlib,glob
matplotlib.use('TkAgg') #needed as OVITO uses the default QT5 backend
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
import numpy as np
from copy import deepcopy as dc
from ase.neighborlist import neighbor_list
import ase.io
from os import system,cpu_count
import argparse,time
import os.path

from ovito.io import import_file,export_file
from ovito.io.ase import *
from ovito.vis import Viewport

import ase.data #for covalent radii

from ovito.modifiers import * #ConstructSurfaceModifier,VoronoiAnalysisModifier,ExpressionSelectionModifier,ColorCodingModifier

from tqdm import tqdm # to show progress bars for the for loops, but makes them slower.
from spglib import find_primitive,standardize_cell,get_spacegroup #, niggli_reduce #niggli_reduce from ASE-tools conflicts with that from spglib.

import scipy.stats as stats
import seaborn as sns
import pickle

from soprano.collection import AtomsCollection
import glob

from sklearn.preprocessing import normalize
from matplotlib.gridspec import GridSpec

#Symm groups do not work as a Gene, as the values are not numeric (string), 
#so as a workaround we define a dictionary for all unique symmetry groups corresponding to the list of sturtcures.
from soprano.analyse.phylogen import Gene
#from soprano.analyse.phylogen.Gene import parsegene_symm
#from soprano.properties.basic import CalcSymm
from soprano.properties import AtomsProperty


from p_tqdm import * #p_imap, p_map,p_umap, p_uimap #u:unordered, i:iterator #pip install p_tqdm --user     #https://pypi.org/project/p-tqdm/
from collections import defaultdict


from subprocess import Popen,PIPE # as popen # check_output
from sys import stdout,version_info
def Popen4(cmd,timeout=None):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate(timeout=timeout)
    #out=proc.stdout.readline() #does not work.


from ase.io import write,read
from ase.io import pov
from os import system
#from ase.utils import hsv #no need
def ase_save_fig(atoms,fname,coords,radius=1.1,view='top'):#options: top, ortho, side
    #TODO: add energy to filename
    #Make colors for each atom
    #colors = hsv(atoms.positions[:, 0])
    cls=[]
    for crd in coords: #colors in rgb
        if crd==2: cls.append([0,0,1]) #sp: blue
        elif crd==3:cls.append([0,1,0]) #sp2: green
        elif crd==4:cls.append([1,0,0]) #sp3: red
        else: cls.append([0,0,0]) #black
    #colors=None #use coord no for colouring
    colors=cls
    
    # Textures
    #tex = ['jmol',] * 288 + ['glass',] * 288+ ['ase3',] * 288 + ['vmd',] * 288
    tex = ['jmol',] * 288 +  ['ase3',] * 288 + ['vmd',] * 288

    #found using ASE-GUI menu 'view -> rotate'
    if view=='top':rotation = '0x, 0y, 0z' # top (z) view
    elif view=='xside':rotation = '-90x, -90y, 0z' # xside view
    elif view=='yside':rotation = '90x, 0y, 90z' # yside view
    elif view=='ortho':rotation = '-15x, 45y, -10z' # ortho view

    
    # keywords
    kwargs = { # Keywords that exist for eps, png, and pov
    'rotation': rotation,
    'colors': colors,
    #'radii': None,
    'radii'         : .65, # float, or a list with one float per atom
    'show_unit_cell': 0,   # 0, 1, or 2 to not show, show, and show all of cell
    }

    extra_kwargs = { # For povray files only
    'display'      : False, # Display while rendering
    'pause'        : False, # Pause when done rendering (only if display)
    'transparent'  : False, # Transparent background
    #'canvas_width' : 1024,  # Width of canvas in pixels #Def: None
    'canvas_width' : 640,  # Width of canvas in pixels #Def: None
    'canvas_height': None,  # Height of canvas in pixels #Def: None
    #'camera_dist'  : 50.,   # Distance from camera to front atom
    'camera_dist'  : 7.5,   # Distance from camera to front atom
    'image_plane'  : None,  # Distance from front atom to image plane
                            # (focal depth for perspective)
    #'camera_type'  : 'perspective', # perspective, orthographic, ultra_wide_angle
    'camera_type': 'orthographic',  # perspective, ultra_wide_angle
    'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
    'area_light'   : [(2., 3., 40.) ,# location
                      'White',       # color
                      .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
    'background'   : 'White',        # color
    'textures'     : tex, # Length of atoms list of texture names #Def: None
    'celllinewidth': 0.05, # Radius of the cylinders representing the cell
    'transmittances': None,  # transmittance of the atoms

    # use with care - in particular adjust the camera_distance to be closer
    'depth_cueing': 1,  # fog a.k.a. depth cueing
    'cue_density': 1e-3,  # fog a.k.a. depth cueing #def: 5e-3
    'bondlinewidth': 0.10,  # radius of the cylinders representing bonds
    #'bondatoms': pov.get_bondpairs(atoms, radius=1.1) ,  # [[atom1, atom2], ... ] pairs of bonding atoms  #Def: []
    'exportconstraints': False,  # honour FixAtoms and mark relevant atoms?
    }

    kwargs.update(extra_kwargs)

    # Make the raytraced image
    #To draw bonds between atom pairs based on (covalent radii*radius).
    if 1:  bp=[ [p[0],p[1]] for p in pov.get_bondpairs(atoms, radius=radius)] #radius =scaling factor for the covalent radii,Def: 1.1
    else: bp=[]
    write('%s.pov'%fname, atoms, run_povray=True, bondatoms=bp,**kwargs) #works
    #write('%s.png'%fname, atoms, **kwargs) #kwargs does not work works
    #atoms.write('%s.png'%fname) #works but tranparent background
    
    #system('convert -delay 40 -loop 100 `ls -v neb-movie-*.png` neb.gif')
    system('rm -f %s.ini %s.pov'%(fname,fname))
        

    #povray  +I neb-movie.1.pov +Oneb-movie.1.png +W1024 +H768 +V -D +FN +Q9 +P +UD +UL +UV +A +AM2 +UA
    


from matplotlib import colors
def mycmap(c='C0',ifRev=0): #input python color code, def: C0: blue
    x=colors.LinearSegmentedColormap.from_list(
        'incr_alpha', [(0, (*colors.to_rgb(c),0)), (1, c)])
    if ifRev: return x.reversed() #data mapping points must start with x=0 and end with x=1 @We use a reversed cmap, as wwe want to higher-energy structures to have white color.
    else: return x

#from scipy.stats import linregress
#from sklearn.linear_model import LinearRegression
def linear_fit(x,y,method='lstsq'):
    
    if method=='lstsq': # compute by linear least-squares fit. This one is slightly faster.
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y,rcond=None)[0] #m: slope is 2d*MSD

    elif method=='linreg1': #compute by Linear regression.
        reg = LinearRegression(fit_intercept=0, normalize=False).fit(x, y) #from sklearn
        m=reg.coef_
        c=reg.intercept_
    elif method=='linreg2': #another linear regression
        m,c, r_value, p_value, std_err=linregress(x,y) #from scipy.stats
        #r2_vals.append(r_value**2); #p_vals.append(p_value); #stderrs.append(std_err)
        
    elif method=='polyfit': #compute by polynomial fit
                lfit=np.polyfit(x,y,1)
                m=lfit[0]
                c=0.
    return m,c

def boltz_dist(energies,T=298.15,omega=[],ifNorm=0):#Return the occupation probabilities of the configurations at a given temperature based on their energies.
    kb= 8.6173303*10**-5 #Boltzmann constant (eV/K).
    if len(omega)==0:#If the degeneracies are not given explciitly.
        omega=[1 for E in energies]
        
    energies=np.array(energies)

    if 1: #Get relative energies.  Doesn't really matter as long as you use the normalised factors.
        mn=np.min(energies)
        #energies=[i-mn for i in energies]
        energies-=mn
    probs=[]
    for E in energies:   probs.append(np.exp(-E/kb/T))

    if ifNorm:
        #Normalise
        Z=sum(probs) #i.e. partition fnc
        probs=[Pn/Z for Pn in probs]

        #Configurational statistics as given in R. Grau-Crespo et al. J.Phys: COndens. Matter, 19,2007,256201
        E_avg=sum([energies[i]*probs[i] for i in range(len(energies))])
        #print ("\nAverage energy of the sytem in configurational equilibirum,  E=%.5f eV"%E_avg)

        F=-kb*T*np.log(Z)
        #print ("Configurational free energy in the complete space, F=%.5f eV"%F)

        S= (E_avg-F)/T
        #print ("Configurational entropy in the complete space, S=%.5f eV/K"%S)

        Smax=kb*np.log(len(energies))
        #print ("Upper limit of config. entropy, Smax= %.5f eV/K"%Smax)

        #Now count in the degenaricies of the configs.
        Sm=[kb*T*np.log(om) for om in omega] #degeneracy entropy

        Em_bar=[energies[i]-T*Sm[i] for i in range(len(energies))]

        Pm=[np.exp(-Em_bar[i]/kb/T) for i in range(len(energies))]
        Z_bar=sum(Pm)
        #Pm_bar=[(1/Z)*np.exp(-Em_bar[i]/kb/T) for i in range(len(energies))] #reduced  probability for an independent config.
        Pm_bar=[P/Z_bar for P in Pm]

        E_avg=sum([Em_bar[i]*Pm_bar[i] for i in range(len(energies))])

        F=-kb*T*np.log(Z_bar)
        #print ("Configurational free energy in the reduced config. space, F=%.5f eV"%F)

        S= (E_avg-F)/T
        #print ("Configurational entropy in the reduced config. space, S=%.5f eV/K"%S)

        return Pm_bar#,E_avg

    else: return probs

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
        return Popen4('grep -m %d "%s" %s '%(n,key,fname))[0]
    except:
        return ""


def getVol_OVITO(data,radius=1.85):
    """
    # Set up the Voronoi analysis modifier.
    voro = VoronoiAnalysisModifier(
        compute_indices = 0,
        use_radii = True,
        edge_threshold = 0.1
    )
    #pipeline.modifiers.append(voro)
    #pipeline.modifiers.append(ConstructSurfaceModifier(radius = radius))
    """

    try:
        #data.apply(voro)
        data.apply(ConstructSurfaceModifier(radius = radius,identify_regions=1))
    except: return None; #10**-10

    return data.attributes['ConstructSurfaceMesh.filled_volume']

########
# Main method
########
    
parser = argparse.ArgumentParser(description='Script for analysing set of predicted cluster structures (res from AIRSS or any type supported by ASE).')

parser.add_argument('-t', '--type',nargs='*',required=1, help="""Determines the type of analysis. Default: all (i.e -t all). Available options are: 
'coord', 'energy', 'dos', 'shape',  'order', 'rings', 'rdf', 'bondangle','all' """) #,choices=['',"all"]  #default=['all']

parser.add_argument('-rtype', '--rType',default='Kings',choices=['Kings','Guttman','primitive','strong'],help='Type of the shortest distance determnination algorithm to use to determine the rings for a given structure, as implemented in the R.I.N.G.S code. Def: Kings/Franzblau algorithm.')

parser.add_argument('-vol','--clVol', default=False,action='store_true',  help='Use Voronoi Tesselation to find the real volume of cluster (rather than using the box volume), requires Ovito package.')

parser.add_argument('-co','--cutoff', type=float,default=1.85,help='Cutoff distance for coordiantion analysis and Voronoi volume caclulations')


parser.add_argument('-du','-del','--delUnPhysical', default=False,action='store_true',  help='Delte the structures with unphysical connections')

parser.add_argument('-clus','-cluster','--cluster', default=False,action='store_true',  help='Cluster the input structures and map them onto 2D space.')

#parser.add_argument('-oa','-order','--orderAnalysis', default=False,action='store_true',  help='Perform order paramter analysis (anisotropy vs coordination no ratio)')

parser.add_argument('-cl', '--clist',default='all',nargs='*',help='List of cluster size to consider; e.g. -nl 0,5,6:128  [Def: all]')

parser.add_argument('-sl', '-slist','--symmlist',default=['all'],nargs='*',help='List of symmetry (point/space) group  to consider; e.g. -sl C1 Ih  [Def: all]')

parser.add_argument('-se', '-sexcl','--symm_exclude',default=[],nargs='*',help='List of symmetry (point/space) group  to exclude; e.g. -se C1 Ih  [Def: all]')

parser.add_argument('-ecut','--ecut', type=float, help='Energy cutoff for structures [eV/atom]. DEf: None, all structures included')

parser.add_argument('-pres','--presslim', type=float, default=10**10,help='Pressure upper limit to discard the structure. Def: None, all structures included')

parser.add_argument('-gc','-gcmap','--globalCM', default=False,action='store_true',  help='Use the global energies for the coloring of the structures in the plots.')

parser.add_argument('-u','-unite','--unite', default=False,action='store_true',  help="Returns a unique list of structures (stored in 'UNIQUE' folder) by uniting the 'same' structures based on similarity analysis. (Def: No)")

parser.add_argument('-us','-usoap','-uniteSOAP','--uniteSOAP', default=False,action='store_true',  help="Returns a unique list of structures (stored in 'UNIQUE' folder) by uniting the 'same' structures based on similarity analysis. (Def: No)")

parser.add_argument("-mu","--chempots", type=str, nargs='*', help="Chemical potentials of the constitutent atoms for computing the formation energy per atom for comparing different stoichometries.Energies should be per fu. \nUsage: -mu Li=value Si=value ...")

parser.add_argument('-temp','--temp', type=float, default=298.15,help='Temperature [K] to be used in Boltzmann weighting of the structures in analyses. Def: 298.15K')


parser.add_argument('-save','--save', default=False,action='store_true',  help='Save the set of structures used in the current analysis in the args.sfold dir')

parser.add_argument('-sf','-savefig','--savefig', default=False,action='store_true',  help='Save the visualisation (snapshot) of the given system coloured wrt to the coord no in the current dir.')


parser.add_argument('-sfold','--sfold', default='SAVED',  help='Directory to save structures')

parser.add_argument('-np','--nprocs', type=int, default=4,help='No of processes/threads to use in parallel computations. Def: 4') #TODO: update the entries throughout the code!!

parser.add_argument('-nopl','-noplot','--noplot', default=False,action='store_true',  help='Surpress the plotting only generate the data files). Def: No')

parser.add_argument('-noE','--noE', default=False,action='store_true',  help='Surpress filtering out structures based on (positive) energy during the --args.del process . Def: Off ')


#Features from the old code:
parser.add_argument('-uo','-old','--useOld', default=False,action='store_true',  help='Activate the old version of the code to use the features listed below.')

parser.add_argument('-po','-pores','--pores', default=False,action='store_true',  help='Plot porosity against density (for GS structures only), requires Zeo++ network code.')

parser.add_argument('-el','-elastic','--elasticity', default=False,action='store_true',  help='Plot bulk modulus (elasticity) against density (for GS structures only), requires LAMMPS code.')

parser.add_argument('-sh','-shape','--shape', default=False,action='store_true',  help='Plot radius of gyration (Gr) and asphericity (b) against cluster-size/density (for GS structures only), requires LAMMPS code.')


args = parser.parse_args()

if args.shape or args.elasticity or args.pores or args.useOld:
    print ('Using old code');useOld=1
args.type=set(args.type)
#TODO: Integrate the features from the older code into the new one (instead of calling the useOld flag.

initT=time.time()

lammps_exe="mpirun -np 8 ~/APPS/lammps-August19/lmp_openmpi" #take as user input.
plt.ion() #interactive plotting on


#Determine which cluster sizes to be considered.
if args.clist:
    clist=[]
    args.clist=''.join(args.clist)

    yy=args.clist.split(',')

    for y in yy:
        if y=='all':clist='all';break
        x=y.split(':')
        if len(x)==1:
            clist.append(int(x[0]))
        elif len(x)==2:
            x0=int(x[0])
            x1=int(x[1])
            clist.extend(list(range(x0,x1+1)))
        else:continue
            

cpots={};
if args.chempots:
    for pot in args.chempots:
        x=pot.split("=")
        try:
            if x[0] in cpots: cpots[x[0]].append(float(x[1]))
            else:cpots[x[0]]=[float(x[1])]
        except: print ("Error in parsing "),pot;exit()

    for pot in cpots.keys():  cpots[pot]=np.mean(cpots[pot])
    

cutoff=args.cutoff #Angstroem for determining the coord no

if args.save:
    system('rm -f %s/*.res; mkdir -p %s'%(args.sfold,args.sfold))
    with open('./%s/README'%args.sfold,'w') as f:
        if args.ecut: f.write('#This dataset contains structures with Ecut=%.5f eV/atom and cluster list: %s\n'%(args.ecut,args.clist))
        else:f.write('#This dataset contains structures with no Ecut and cluster list: %s\n'%(args.clist))
if args.type.intersection(['r','ring','rings','all']):  system('mkdir -p data')


#Determine the label based on the CWD
f=os.getcwd()
if 'DFT'in f:label='DFT/PBE'
elif 'EDIP'in f:label='EDIP'
elif 'LCBOP'in f:label='LCBOP'
elif 'AIREBO-M'in f:label='AIREBO-M'
elif 'AIREBO'in f:label='AIREBO'
elif 'REBO'in f :label='REBO'
elif 'ReaxFF'in f:label='ReaxFF'
elif 'GAP-ML-2020'in f:label='GAP20'
elif 'GAP-ML-2017'in f:label='GAP17'
elif 'TERSOFF'in f:label='TERSOFF'
else: label=''
print('Type of potential: ',label)
cwd_label=label

if cwd_label=='DFT/PBE':
    if args.chempots==None:args.chempots=1;cpots={'C':-1.3816}#TODO: Delte this in the final version
#elif cwd_label=='GAP20': args.chempots=1;cpots={'C':-0.9444}

if cwd_label in ['ReaxFF','GAP17',]: Escale=23.0 #kcal/mol to eV (i.e. units real to metal)
else:Escale=1.0  

if args.delUnPhysical: #Filter the structures with more than one clusters in the beginning, using the OVITO tools.
    #from ovito.data import *
    #system('rm -rf NOT-USED; mkdir -p NOT-USED')
    system('mkdir -p NOT-USED')
    
    
    #cutoff=1.70;
    cutoff=args.cutoff   #1.6 is too low!!
    modifier=ClusterAnalysisModifier(cutoff=cutoff,sort_by_size=0, compute_com=0,unwrap_particles=1, compute_gyration=0) #unwrap_particles for making clusters contiguous through PBC; neighbour mode

    modifier2 = CoordinationAnalysisModifier(cutoff = cutoff*0.9, number_of_bins = 200,partial=0)
    
    cnt=0
    #inpfs=glob.glob("./*.res",recursive=0)
    #for inpf in glob.iglob("./**/partial_charges_bader.xyz",recursive=1):
    #for inpf in tqdm(inpfs):  #TODO: use Parallel u_pqtm
    #for inpf in glob.iglob("./*.res",recursive=0):
    def filterData(inpf):
        cnt=0
        try:atoms=ase.io.read(inpf)
        except:system('mv %s NOT-USED'%inpf);return None; #continue
        odata=ase_to_ovito(atoms)
        
        odata.apply(modifier)
        odata.apply(modifier2)
        coord=odata.particles.coordination
        cluster_sizes = odata.tables['clusters'].y
        noClus=len(cluster_sizes)
        #print(inpf,noClus)
        if noClus>1:
            #print(inpf,noClus)
            system('mv %s NOT-USED'%inpf)
            cnt+=1;return cnt; #continue

        #Filter out the ones with very high coordination (unphysical)
        if np.count_nonzero(coord>5)>0:
            system('mv %s NOT-USED'%inpf)
            cnt+=1;return cnt; #continue
            
        try:E=atoms.get_potential_energy()
        except:E=1000
        if not args.noE and E>0: #filter out positive energy ones (non-physical)
            None
            print(inpf,E)
            system('mv %s NOT-USED'%inpf)
            cnt+=1;return cnt; #continue

        return cnt

    #exit()

    inpfs=glob.glob("./*.res",recursive=0)
    for res in p_imap(filterData,inpfs,num_cpus=args.nprocs):
        if res != None:
            cnt+=res
            
    if cnt>0:print ('%d unphyiscal structures have been filtered out...'%cnt)
    print("Elapsed time: %.2f sec."%( time.time()-initT))

    #exit()

if args.unite:
    import glob
    #from matador.hull import QueryConvexHull
    #from matador.query import DBQuery
    from matador.scrapers.castep_scrapers import res2dict
    from matador.fingerprints.similarity import get_uniq_cursor
    from matador.export import doc2res
    from matador.utils.cursor_utils import filter_cursor, get_array_from_cursor
    from matador.utils.chem_utils import get_root_source
    
    cursor = [res2dict(inpf)[0] for inpf in glob.iglob("./*.res",recursive=0)] #args.inpf]
    #print(cursor)

    print('Filtering unique structures... {}'.format(len(cursor)))
    uniq_list,dupe_dict,_,_ = list(get_uniq_cursor(cursor, debug=0,sim_tol=0.1,energy_tol=5e-3)) #def:energy=1e-2,sym=0.1; working sets: 0.1 and E=1e-3 (skips more structures); 0.01 and 1e-2 / 1e-3
    #ucursor = [cursor[ind] for ind in uniq_list] #Needed here?
    #print('Final cursor length... {}'.format(len(ucursor)))
    #print (uniq_list)
    #uniq_list=np.array(uniq_list) #Needed for the min-energy reorder part below
    dupe_list=[]
    for key in dupe_dict.keys(): #keys and values are based on the indices from the original list
        #Note: one can either use 'get_root_source(cursor[key]' or cursor[key]['source'], where cursor is dict with the keys: 
        #vals=list(dupe_dict[key])
        vals=dupe_dict[key]
        if len(vals)>0:
            try:
                x=[get_root_source(cursor[key])]
                x.extend([get_root_source(cursor[val]) for val in vals if cursor[val]['space_group']==cursor[key]['space_group']]) #we have some fail-safe filters here (i.e. symmetry checks) #if cursor[val]['space_group'] not in ['C1','P1'] and
               
                if len(x)>1:
                    """
                    if 0: # make sure to save the minimum-energy from a given duplicate set (not a random one). #TODO: fix this, not working .
                        x=np.array(x)
                        Es=[cursor[key]['enthalpy']]
                        Es.extend([cursor[val]['enthalpy'] for val in vals])
                        mInd=np.argmin(Es)
                        print(mInd)
                        if mInd>0:
                            uniq_list=np.where(uniq_list==key,vals[mInd-1],uniq_list)
                            x[[0,mInd-1]]=x[[mInd-1,0]] #Place the lowest energy str from the set to first place
                    """
                    dupe_list.append(x)

            except:raise Exception; continue

    #print(['.res '.join(x)+'.res' for x in dupe_list])

    print('Unique structure set is saved in UNIQUE folder...')
    ucursor = [cursor[ind] for ind in uniq_list]
    print('Final cursor length... {}'.format(len(ucursor)))
    dname = 'UNIQUE' ; os.system('rm -rf %s; mkdir -p %s'%(dname,dname))
    for i,doc in enumerate(ucursor):
        #print(doc)
        doc2res(doc, dname + '/' + get_root_source(doc),spoof_titl=0) #when spoof_titl=0, the original header is kept.

    for i,dupe in enumerate(dupe_list):
        repNo=len(dupe)
        fname=dupe[0]
        if 0: print(repNo, '.res '.join(dupe)+'.res') #verbose
        
        # update the number of hits in the res header line for a correct cryan analysis. 
        system("sed -i 's/n - 1/n - %d/g' %s/%s.res"%(repNo,dname,fname))
        system("sed -i 's/n -1/n - %d/g' %s/%s.res"%(repNo,dname,fname))
               
        
    print("Elapsed time: %.2f sec."%( time.time()-initT))

    #system('cd %s; ca -r -cl  | sort -n  -k 5 -k 4 | head -n 50;cd ..'%dname)


    #exit()

if args.uniteSOAP:
    import ase.io
    from ase.io import read
    from ase.build import molecule
    from ase import Atoms

    import glob

    from dscribe.descriptors import SOAP
    from dscribe.kernels import REMatchKernel

    from ase.build import molecule

    from sklearn.preprocessing import normalize
    from sys import exit
    import time

    initT=time.time()
    inpfs = glob.glob('./*.res') 

    print('Reading %d input files...'%len(inpfs))

    """
    from joblib import parallel_backend,Parallel, delayed
    with parallel_backend('loky', n_jobs=8): #'threading'
        structures=Parallel()(delyaed=([ase.io.read(x) for x in inpfs])
    """
    
    
    def readASE(inpf):
        #structures=[ase.io.read(x) for x in inpfs]
        return ase.io.read(inpf) #structures

    iterator = p_uimap(readASE,inpfs,num_cpus=args.nprocs) #def: all CPUs used. #slightly faster than direct mapping (p_map or p_umap)!! unordered is slightly fater than ordered.
    structures=[x for x in iterator]

    #structures=p_umap(readASE,inpfs,num_cpus=4) #def: all CPUs used.
    #structures=p_map(readASE,inpfs,num_cpus=4) #def: all CPUs used.
    
    #structures=[ase.io.read(x) for x in inpfs] #or use serial reading.
    #"""
    
    print("Elapsed time: %.2f sec."%( time.time()-initT))

    print('\nLength of initial cursor: ', len(structures))


    # Let's create a list of structures and gather the chemical elements that are
    # in all the structures.
    species = set()
    energies=[];symmdata=[];nHits=[]
    for x in structures:#Atoms objects
        #dict_keys(['name', 'pressure', 'energy', 'spacegroup', 'times_found']) for x.info
        species.update(x.get_chemical_symbols())
        try:energies.append(x.get_potential_energy()/len(x)/Escale) #E/atom
        except:energies.append(0.0)
        try:symmdata.append(x.info['spacegroup'])
        except:symmdata.append('C1')
        try:nHits.append(x.info['times_found'])
        except:nHits.append(1)

    energies=np.array(energies)
    symmdata=np.array(symmdata)
    nHits=np.array(nHits)
    usymm=np.unique(symmdata)
        

    #Setup the descriptor object. The exact setup depends on the used descriptor and your use case. Notice that typically you will want to instantiate only one descriptor object which will handle all structures in your dataset. You should read the original articles or the specific tutorials to understand the meaning of different settings. For machine learning purposes you may also want to cross-validate the different settings to find the best-performing ones.
    # Let's configure the SOAP descriptor.
    from dscribe.descriptors import SOAP

    sparse=1
    soap = SOAP(
        species=species,
        periodic=False,
        rcut=3.7, #3.7A used for C by Gabor, Deringer, 2.0,3.5 and 7.0A in Partok,Gabor (SketchMap)
        nmax=8, #both 16 in Gabor et al. PRB, 2020, max=9 is allowed for Gaussian type basis fncs
        lmax=8,
        average="inner",
        sparse=sparse
    )

    #The averaging mode over the centers of interest. Valid options are:
    #”off”: No averaging.
    #”inner”: Averaging over sites before summing up the magnetic quantum numbers: pZ1,Z2nn′l∼∑m(1n∑ici,Z1nlm)∗(1n∑ici,Z2n′lm)
    #”outer”: Averaging over the power spectrum of different sites: pZ1,Z2nn′l∼1n∑i∑m(ci,Z1nlm)∗(ci,Z2n′lm)
    #One way of turning a local descriptor into a global descriptor is simply by taking the average over all sites. DScribe supports two averaging modes: inner and outer. The inner average is taken over the sites before summing up the magnetic quantum number. Outer averaging instead averages over the power spectrum of individual sites. In general, the inner averaging will preserve the configurational information better but you can experiment with both versions.  

    # Let's create SOAP feature vectors for each structure
    print('Creating SOAP feature vectors for all structures...')
    feature_vectors = soap.create(structures, n_jobs=cpu_count())
    print("Elapsed time: %.2f sec."%( time.time()-initT))
    if sparse: feature_vectors=np.array(feature_vectors.toarray())
    #As SOAP is a local descriptor, it also takes as input a list of atomic indices or positions. If no such positions are defined, SOAP will be created for each atom in the system. 

    #The result will be a feature vector and not a matrix, so it no longer depends 
    #on the system size. This is necessary to compare two or more structures with
    #different number of elements. We can do so by e.g. applying the distance 
    #metric of our choice.
    from scipy.spatial.distance import pdist, squareform  #This is serial and reuires lots of mem
    from sklearn.metrics.pairwise import pairwise_distances,pairwise_distances_chunked

    import numpy as np

    print('Calculating similarity matrix...')
    #distance = squareform(pdist(feature_vectors)) #from scipy, serial; def metric is euclidean
    #os.environ['OMP_NUM_THREADS']='4'
    distance = pairwise_distances(X = feature_vectors, metric = 'euclidean', n_jobs = 4) #-1 use all CPUs. #from sklearn, parallel  ; def metric is euclidean; #CHECK: pairwise_distances_chunked for a low-memory version
    #distance = next(pairwise_distances_chunked(X = feature_vectors, metric = 'euclidean', n_jobs = 4),working_memory=)

    #print('Similar structures',)
    doneList=[];dupe_dict={}
    thresh=5.0 #similarity threshold [in A] #7.0 A kadar ciktim

    #TODO: Only compare structures of the same size!! Do a cSize loop!!
    for i,dist in enumerate(distance): #iterate over structures
        if i in doneList: continue #do not redo same structures.
        slist=np.where(dist<thresh)[0]
        doneList.extend([x for x in slist])
        
        if len(slist)>1:
            #Now add the other structures similar to other structures in the similarity list (slist)
            new_slist=dc(slist)
            for sim in slist:
                if sim==i:continue
                cdist=distance[sim]
                clist=np.where(cdist<thresh)[0]
                new_slist=np.append(new_slist,clist)

            slist=np.unique(new_slist)#[0]
            doneList.extend(slist)
            
            Es=energies[slist]
            if abs(np.min(Es)-np.max(Es))>0.01: #catch if any outlier in the list in terms of energy per atom.
                if len(slist)==2:
                    #print('Not really a similar set');
                    continue #this is not a real similar set.

                else: #Locate the outlier(s) in a str list of 3 or more items.

                    slist_new=slist.tolist()
                    #for j,E in enumerate(Es):
                    for j in range(1,len(Es)):
                        #if j==i:continue
                        
                        if not np.isclose(Es[0],Es[j],atol=0.01) :
                            print ('Diff energy:',i,inpfs[i],symmdata[0],Es[0],symmdata[j],Es[j])
                            doneList.remove(slist[j])
                            slist_new.remove(slist[j])
                            continue
                        if not symmdata[0]!=symmdata[j]:
                            print ('Diff symm group:',i,inpfs[i],symmdata[0],symmdata[j])
                            doneList.remove(slist[j])
                            slist_new.remove(slist[j])
                        
                            #print('after deletion: ',i, slist_new);
                    slist=dc(slist_new)
                    
        #Use the min-energy structure as the representative one.
        mind=np.argmin(energies[slist])
        dupe_dict[slist[mind]]=slist

    print('Length of final cursor: ', len(dupe_dict.keys()))
    
    print("Elapsed time: %.2f sec."%( time.time()-initT))


    dname = 'UNIQUE-SOAP' ; os.system('rm -f %s/*.res; mkdir -p %s'%(dname,dname))
    print('\nWriting structure files in %s...'%dname)
    for key in dupe_dict:
        slist=dupe_dict[key]
        fname=inpfs[key]
        nHit=nHits[key]-1 #Not to lose the already-existing degenracy info
        nSim=len(slist) #Get the new degeneracy info based on the current similarity analysis.
        repNo=nHit+nSim
        if 0: print(key,[inpfs[sl] for sl in slist],repNo)

        # update the number of hits in the res header line for a correct cryan analysis.
        system("cp %s %s/%s"%(fname,dname,fname))
        if nSim>1:
            system("sed -i 's/n - 1/n - %d/g' %s/%s"%(repNo,dname,fname))
            system("sed -i 's/n -1/n - %d/g' %s/%s"%(repNo,dname,fname))


    print("Elapsed time: %.2f sec."%( time.time()-initT))    



#All PG numbers in increasing order (from https://en.wikipedia.org/wiki/Point_group)
allPGDict={'C1':1, 'Ci':1, 'S2':1, 'Cs':2,'C1v':2,'C1h':2,'C12':2,'C2':3,'C3':4,'C3i':4,'C4':5,'C5':6, 
           'C6':7,'C7':7,'Cn':8,'C2v':9,'C3v':10,'C4v':11,'C5v':12,'C6v':13,'Cnv':14,
          'C2h':15,'C3h':16,'C4h':17,'C5h':18,'C10':18,'C6h':19,'Cnh':20,'S4':21,'S6':22,
           'S8':23,'S10':24,'S12':25,'S2n':26,'D2':27,'D3':28,'D4':29,'D5':30,'D10':30,
           'D6':31,'Dn':32,'D2h':33,'D3h':34,'D4h':35,'D5h':36,'D6h':37,'Dnh':38,'D(inf)h':38,
           'D2d':39,'D3d':40,'D4d':41,'D5d':42,'D6d':43,'Dnd':44,'T':45,'T23':45,'Th':46,
           'Td':47,'O':48,'Oh':49,'I':50,'Ih':51   }

class CalcSymm(AtomsProperty): #Added by BK

    """
    CalcSymm

    Property representing the energy calculated by an ASE calulator

    """

    default_name = 'calc_symm'
    default_params = {}

    @staticmethod
    def extract(s):
        try:
            sg=s.info['spacegroup']
            #return symmdict[sg]
            return allPGDict[sg]
        except:
            return None

def parsegene_symm(c): #added by BK
    return np.array(CalcSymm.get(c))

#gene_symm= Gene(name='symm',parser=parsegene_symm) #added by BK
#print(gene_symm)

from dscribe.descriptors import SOAP
from dscribe.kernels import REMatchKernel

#from ase.build import molecule

from sklearn.preprocessing import normalize
from sys import exit
import time


class CalcSoap(AtomsProperty): #Added by BK

    """
    CalcSoap

    Property representing the SOAP parameter for a structure (given as ASE Atoms) as calculated by Dscribe.

    """

    default_name = 'calc_soap'
    default_params = {} #{'periodic':False,'rcut':3.7,'nmax':8,'lmax':8,'avarage':'inner','sparse':0}  
    
    @staticmethod
    def extract(s):
        try:
            #s (self) stands for Atoms object
            #sg=s.info['spacegroup']
            #return symmdict[sg]
            species=s.get_chemical_symbols()
            sparse=1
            soap = SOAP(
                species=species,
                periodic=0, #params['periodic'],
                rcut=3.7, #rcut, #3.7A used for C by Gabor, Deringer, 2.0,3.5 and 7.0A in Partok,Gabor (SketchMap)
                nmax=8, #nmax, #both 16 in Gabor et al. PRB, 2020, max=9 is allowed for Gaussian type basis fncs
                lmax=8, #lmax,
                average='inner', #average,
                sparse=sparse, #sparse
            )

#The averaging mode over the centers of interest. Valid options are:
#”off”: No averaging.
#”inner”: Averaging over sites before summing up the magnetic quantum numbers: pZ1,Z2nn′l∼∑m(1n∑ici,Z1nlm)∗(1n∑ici,Z2n′lm)
#”outer”: Averaging over the power spectrum of different sites: pZ1,Z2nn′l∼1n∑i∑m(ci,Z1nlm)∗(ci,Z2n′lm)
#One way of turning a local descriptor into a global descriptor is simply by taking the average over all sites. DScribe supports two averaging modes: inner and outer. The inner average is taken over the sites before summing up the magnetic quantum number. Outer averaging instead averages over the power spectrum of individual sites. In general, the inner averaging will preserve the configurational information better but you can experiment with both versions.  

            # Let's create SOAP feature vectors for each structure
            feature = soap.create(s, n_jobs=1)#cpu_count())
            #if sparse: feature=feature.tolist()
            #As SOAP is a local descriptor, it also takes as input a list of atomic indices or positions. If no such positions are defined, SOAP will be created for each atom in the system. 

            return feature #[0] #return a 1D array instead of 2D array
        except:
            return None

def parsegene_soap(c): #added by BK
    #return np.array(CalcSoap.get(c))
    return np.array([x.toarray()[0] for x in CalcSoap.get(c)])

#gene_soap= Gene(name='soap',parser=parsegene_soap,weight=1.0, params={}) #added by BK
#print(gene_soap)

class CalcSymm(AtomsProperty): #Added by BK

    """
    CalcSymm

    Property representing the energy calculated by an ASE calulator

    """

    default_name = 'calc_symm'
    default_params = {}

    @staticmethod
    def extract(s):
        try:
            sg=s.info['spacegroup']
            #return symmdict[sg]
            return allPGDict[sg]
        except:
            return None

def parsegene_symm(c): #added by BK
    return np.array(CalcSymm.get(c))

if args.cluster:
    from soprano.analyse.phylogen import Gene, PhylogenCluster

    # To carry out the analysis we need to define a PhylogenCluster object. This will need as input some Gene objects.
    # The phylogenetic nomenclature is just an analogy to the way phylogenetic analysis is carried out in biology.
    # Ideally, we're trying to give each structure a "DNA" of sorts, then compare them amongst themselves to find
    # which ones are more closely related.
    # Finding the right properties to use to distinguish between structures is key here. In this examples it's pretty
    # simple but we'll still illustrate a couple different ways to get there
    #######
    #Prepare the input data
    ########
    print('Preparing the input data (collection of structures)')
    inpfs = glob.glob('./*.res') #Use RES files, FASTER than CIFs, and RES files store the energy and symm info!!
    aColl = AtomsCollection(inpfs, progress=True) # "progress" means we will visualize a loading bar

    try:energies=np.array([x.get_potential_energy()/len(x)/Escale for x in aColl]) #E/atom
    except:energies=np.array([0.0 for x in aColl])
    try:symmdata=np.array([x.info['spacegroup'] for x in aColl])
    except:symmdata=np.array(['C1' for x in aColl])
    usymm=np.unique(symmdata)
    
    # This gene represents the length of the three lattice parameters
    gene_abc = Gene(name='latt_abc_len', weight=1.0, params={})

    # This gene represents the linkage list property as seen in tutorial 2
    gene_lnk = Gene(name='linkage_list', weight=1.0, params={})#{'size':10})#def:all nearest neighbours/all pair distances e.g. C(60,2) pairs

    gene_coord = Gene(name='coord_histogram', weight=1.0, params={'s1':'C','s2':'C','max_coord':6})
    gene_ener= Gene(name='energy')

    gene_symm= Gene(name='symm',parser=parsegene_symm) #added by BK
    gene_soap= Gene(name='soap',parser=parsegene_soap)#,weight=1.0, params={'periodic':False,'rcut':3.7,'nmax':8,'lmax':8,'avarage':'inner','sparse':0}) #added by BK



    #from ./soprano/properties/linkage/linkage.py
    #see GeneDictionary in ./soprano/analyse/phylogen/genes.py for a full list
    #name='bonds' ; 'linkage_list' ; 'molecules';'molecule_n';'molecule_mass';
    #'molecule_com';'hydrogen_bonds_n';'dihedral_angle_list'; 'bond_graph';
    # 'energy'

    # We can try these separately or together
    #phClust1 = PhylogenCluster(aColl, [gene_abc])
    #phClust2 = PhylogenCluster(aColl, [gene_lnk])
    #phClust3 = PhylogenCluster(aColl, [gene_abc, gene_lnk]) # ORIGINALversion: In this case they get chained together,
                                                            # and the relative weights are used
    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_ener,gene_symm])
    #phClust3 = PhylogenCluster(aColl, [gene_coord,gene_ener,gene_symm])

    #phClust3 = PhylogenCluster(aColl, [gene_coord,gene_lnk,gene_ener])
    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_symm,gene_ener])

    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_coord,gene_ener,gene_symm])
    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_coord,gene_symm,gene_ener])
    #phClust3 = PhylogenCluster(aColl, [gene_symm,gene_ener])

    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_coord,gene_ener])
    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_ener])
    #phClust3 = PhylogenCluster(aColl, [gene_lnk,gene_coord])
    #phClust3 = PhylogenCluster(aColl, [gene_coord])
    #Energy and symm Genes can't be used alone.

    phClust3 = PhylogenCluster(aColl, [gene_soap])
    #phClust3 = PhylogenCluster(aColl, [gene_soap,gene_lnk])
    #phClust3 = PhylogenCluster(aColl, [gene_soap,gene_ener])
    #phClust3 = PhylogenCluster(aColl, [gene_soap,gene_ener,gene_symm])
    #phClust3 = PhylogenCluster(aColl, [gene_soap,gene_symm])

    # Here's a summary of the generated "genes"
    genes, genes_info = phClust3.get_genome_vectors_norm()

    print("---- Genetic strings for each structure (normalised) ----\n")
    #print('\n'.join([str(g) for g in genes]), '\n')
    print(genes,np.shape(genes))
    print('Info:\t', genes_info, '\n\n') # This tells us which genes are present and how long the respective fields are

    ############
    # Now Generate a 2D mapping
    ########
    print('Now reducing the dimensionality, i.e. mapping the data points onto 2D...')
    
    # Here we're using the total principal component method,
    #print(help(phClust3.create_mapping))
    method='total-principal' #default
    #method='clafic'  #Diagonal plot
    #method='fukunaga-koontz' #better? very slow
    #method='optimal-discriminant' #not working, singular matrix
    mapx, mapy = phClust3.create_mapping(method=method)
    np.savetxt('data-cluster.dat',np.array([[mapx[i].real,mapy[i].real] for i in range(len(mapx))]),header='Descriptors: %s'%(', '.join([x[0] for x in genes_info])),fmt='%.5f')
    
    #print (mapx,mapy)
    #TODO: use the scikit mapping algorithms (e.g. kernel PCA and alike) instead of those from SOPRANO
    #https://scikit-learn.org/stable/modules/unsupervised_reduction.html 
    #from sklearn.cluster import KMeans
    #from sklearn.decomposition import PCA
    print(np.shape(mapx),np.shape(mapy))

    minE=min(energies);maxE=max(energies)
    sizes=np.array([1-(E-minE)/(maxE-minE)*0.95 for E in energies])#0.95 scaling to keep the largest Energy one too

    fig, ax = plt.subplots(1, 1)
    #fig = plt.figure(figsize=(9.,6.))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    ax2=ax1.twinx() 

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    #print(colors)
    markers=['s','o', '+', 'x', '*', 'v', '^','<','>']
    plt.gca().set_prop_cycle(marker=markers) 

    c=[allPGDict[x] for x in symmdata];c=energies;c=['C%d'%allPGDict[x] for x in symmdata]
    ax1.scatter(mapx, mapy, c=c,cmap='inferno', s=[sizes*50], alpha=0.8); 
    #ax1.scatter(mapx,mapy);quit()
    for i,sy in enumerate(usymm): ax2.plot(sy,label=sy,c='C%d'%allPGDict[sy],ms=3,lw=0,marker='o')
    handles, labels = ax2.get_legend_handles_labels()
    #fig.legend(handles,labels,ncol=1,loc='right',fontsize='x-small',bbox_to_anchor=(1.0, 0.5))#, borderaxespad=0.) 
#fig.legend(handles,labels,fontsize='x-small',ncol=int(len(usymm)/3),loc='upper center', borderaxespad=0) 
    ax2.set_axis_off();ax2.set_visible(0)
    #plt.tight_layout()
    ax1.set_title('Descriptors: %s'%(', '.join([x[0] for x in genes_info])))
    print("Elapsed time: %.2f sec."%( time.time()-initT))

    plt.title(cwd_label)
    if not args.noplot: plt.show(block=1)

#Only do the following analysis when chosen.
if args.delUnPhysical or args.unite or args.uniteSOAP or args.cluster: exit()

    
if not args.useOld:
    cluster_mod=ClusterAnalysisModifier(cutoff=args.cutoff,sort_by_size=0, compute_com=1,unwrap_particles=0, compute_gyration=1) #unwrap_particles for making clusters contiguous through PBC; neighbour mode
    #This modifier provides the option to "unwrap" the coordinates of particles belonging to the same cluster, mapping them from inside the periodic simulation domain to outside in order for the cluster to appear as a contiguous object. Note that the unwrapping will yield meaningful results only for finite clusters. It will lead to undefined results for infinite clusters that are connected to themselves through a periodic cell boundary
    

    system('mkdir -p NOT-USED')
    if args.savefig: system('mkdir -p FIGURES')
    conv=1.6605391 #amu/A^3 to g/cm^3


    if args.type.intersection(['r','ring','rings','all']):
        odir='rings-out'
        system('mkdir -p %s'%odir)
        if 1:
            print('Converting res files into xyz format for the RINGS calculations...')
            system("ase-file-converter-bk.py -i *.res -ot xyz -odir data/ &>/dev/null; #mv *.xyz data/");
            print("Elapsed time: %.2f sec."%( time.time()-initT))    


    
    coord_mod = CoordinationAnalysisModifier(cutoff = args.cutoff, number_of_bins = 200,partial=0)
    
    #print(args.symmlist)
    ##################3
    # 1st input structure loop to get the energies
    ################
    from p_tqdm import * #p_imap, p_map,p_umap, p_uimap #u:unordered, i:iterator #pip inst
    from collections import defaultdict
    allEs=[] #all Energy values in the input structure order
    allStr=[] #keeps the valid str info
    cnt=0;pcnt=0;energy={};energy=defaultdict(list)
    inpfs=glob.glob("./*.res",recursive=0)
    #for inpf in glob.iglob("./*.res",recursive=0):
    #for inpf in tqdm(inpfs):
    def readFiles(inpf):
        allEs=[] #all Energy values in the input structure order
        allStr=[] #keeps the valid str info
        cnt=0;energy={};pcnt=0

        #none=[]
        try: atoms=ase.io.read(inpf)
        except:system('mv %s NOT-USED/'%inpf);return None;#continue
        cSize=len(atoms)
        try:symm=atoms.info['spacegroup']
        except:symm='C1'
        try:pres=atoms.info['pressure']
        except:pres=0.0
        #times_found=atoms.info['times_found']
        if not(clist=='all' or cSize in clist): return None;#continue #No need
        if not args.symmlist[0]=='all' and symm not in args.symmlist:return None;#continue
        if symm in args.symm_exclude:return None;#continue
        if abs(pres)>args.presslim: pcnt+=1;return [cnt,pcnt];#continue
        key=str(cSize)
        try:
        #if 1:
            if not args.chempots:E=atoms.get_potential_energy()/cSize/Escale #energy per atom
            else:
                syms=atoms.get_chemical_symbols()
                atypes,x=np.unique(syms,return_counts=1)
                acounts=dict(zip(atypes,x))
                E=(atoms.get_potential_energy()-np.sum([cpots[x]*acounts[x] for x in atypes]))/cSize/Escale #Formation energy using the reference chemical potentials
        except:E=1000
        if not args.noE and E>0: #filter out positive energy ones (non-physical)
            None
            print(inpf,E)
            cnt+=1;
            system('mv %s NOT-USED'%inpf);
            return [cnt,pcnt] #continue

        if key in energy: energy[key].append(E)
        else: energy[key]=[E]
        allEs.append(E)
        allStr.append(inpf)

        return allEs,allStr,energy,cnt,pcnt

    # End of 1st input structure loop #
    ##########

    inpfs=glob.glob("./*.res",recursive=0)
    #allEs,allStr,energy,cnt,pcnt=p_umap(readFiles,inpfs,num_cpus=4)
    for res in p_imap(readFiles,inpfs,num_cpus=args.nprocs):
        if res != None:
            if len(res)==2:  cnt+=res[0];   pcnt+=res[1]
            else:
                allEs.extend(res[0])
                allStr.extend(res[1])
                for key in res[2].keys():energy[key].extend(res[2][key]) #dict   #energy.update
                cnt+=res[3];  pcnt+=res[4]
            

    #print(allEs,allStr,energy,cnt,pcnt)
    
    if cnt>0:print ('%d structures with no energy info or unrealistic structure (i.e postive total E) have been filtered out...'%cnt)
    if pcnt>0:print ('%d structures with large cell stress/pressure have been filtered out...'%pcnt)
    print("Elapsed time: %.2f sec."%( time.time()-initT))    

    #print(len(allEs),len(allStr))

    
    #######
    # Now do some preprations for plotting and filtering based on relative energies
    #######
    keys=energy.keys()
    #keys=sorted(energy.keys())
    mins={};maxs={}
    for key in keys:
        mins[key]=np.min(energy[key])
        if args.ecut:   maxs[key]=np.max([x for x in energy[key] if np.abs(x-mins[key])<=args.ecut]) #maxs[key]=dc(mins[key])
        else: maxs[key]=np.max(energy[key])
        #print (key,energy[key],mins[key],maxs[key])
        
    if args.ecut!=None and args.ecut>=0.0:
        if args.ecut==0.0: args.ecut=0.0001
        
    #print('mins,maxs:')
    #for key in keys: print(key,energy[key],mins[key],maxs[key])
    
    ###################
    #2nd structure loop  (Now do analysis)
    #################
    #btype={};shape={};energy={};cnt=0
    btypes=[];btdict={};shapes=[];rdfs={};rings={};chains={};cnt=0;blens={};bangs={}
    btdict=defaultdict(list);rdfs=defaultdict(list);chains=defaultdict(list)
    toDel=[]; saveList=[]; #energy={};
    #inpfs=glob.glob("./*.res",recursive=0)
    inpfs=allStr
    #for strind,inpf in enumerate(tqdm(inpfs)):
    ##for strind,inpf in enumerate(inpfs): #glob.iglob("./*.res",recursive=0)):
    def sortData(inp):
        #print(inp)
        inpf=inp[0];strind=inp[1];E=inp[2];#mins=inp[3];maxs=inp[4]#allEs=inp[2]
        btypes=[];btdict={};shapes=[];rdfs={};rings={};chains={};cnt=0;blens={};bangs={}
        toDel=[]; saveList=[]; #energy={};
        #try:
        atoms=ase.io.read(inpf)
        #except:system('mv %s NOT-USED'%inpf);continue
        odata=ase_to_ovito(atoms);odata_cp=odata.clone()
        cSize=len(atoms)
        if not(clist=='all' or cSize in clist): return None; #continue #No need
        key=str(cSize)
        try:seed=inpf.split('/')[-1].split('.')[0]
        except:print ('Error with file name: %s'%inpf);return None; #continue
        
        #try:E=allEs[strind]
        #except:print('Error in index %d %s'%(strind,inpf))
        
        #Determine/filter-out the high-energy structures (i.e. the ones above the cutoff energy).
        if args.ecut!=None and args.ecut>=0.0:
            #if args.ecut==0.0: args.ecut=0.0001
            if abs(E-mins[key])>args.ecut: #then pass this high energy one
                toDel.append(strind);return [toDel,saveList]; #continue
                
        #Save the current structure set in SAVED folder.
        if args.save:
            system('cp %s.res %s-final.res %s/ 2> /dev/null '%(seed,seed,args.sfold))
            saveList.append(inpf);return [toDel,saveList]; #continue

        if args.clVol: vol=getVol_OVITO(odata,1.85)#args.cutoff)#Get the cluster volume (rather than using the volume of the box) via alpha-shape Tesselation (i.e. volume covered by the Voronoi surface) implemented in OVITO visualiser.
        else:vol=atoms.get_volume() #in A^3
        
        try:
            if vol!=0.0: dens=np.sum([atoms[at].mass for at in range(len(atoms))])/vol*conv
            else: dens=0.0
        except:dens=0.0

        #print('my vol: ',vol,'my density: ',dens,atoms)
        

        #Coord analysis is switched on by default.
        odata.apply(coord_mod)#this is needed in each case.
        #if args.type.intersection(['c','coord','all']):odata.apply(coord_mod)
        if args.type.intersection(['s','shape','all']): odata.apply(cluster_mod)
        
        coord=odata.particles.coordination #for each atom
        sp=np.count_nonzero(coord==2)
        sp2=np.count_nonzero(coord==3)
        sp3=np.count_nonzero(coord==4)
        oth=np.count_nonzero((coord<=1) | (coord>4 ))
        #if args.delUnPhysical and oth>0: continue # just skip the structures with unphysical coordinations.
        tot=float(sp+sp2+sp3+oth)
        #if tot==0: btypes.append([cSize,0.,0.,0.,0.,dens])
        #else:btypes.append([cSize,sp3/tot,sp2/tot,sp/tot,oth/tot,dens])
        btypes.append([cSize,sp3/tot,sp2/tot,sp/tot,oth/tot,E,dens])
        if key in btdict:btdict[key].append([sp3/tot,sp2/tot,sp/tot,oth/tot,dens,np.mean(coord)])
        else:btdict[key]=[[sp3/tot,sp2/tot,sp/tot,oth/tot,dens,np.mean(coord)]]

        #SAVED-FIGURES
        if args.savefig:
            os.chdir('./FIGURES')
            sg=atoms.info['spacegroup'].replace('(','').replace(')','')
            fname="%s/C%d_%.2feV_%s-%s"%('.',len(atoms),E,sg,seed)
            colors=np.array([cr for cr in coord])
            if not os.path.exists(fname+'_top.png'):
                ase_save_fig(atoms,fname+'_top',colors,radius=0.9,view='top') #1.21*sum(cov_radii)=1.85A
                ase_save_fig(atoms,fname+'_xside',colors,radius=0.9,view='xside')
                ase_save_fig(atoms,fname+'_yside',colors,radius=0.9,view='yside')
                ase_save_fig(atoms,fname+'_ortho',colors,radius=0.9,view='ortho')
                
            os.chdir('../')

       
        if args.type.intersection(['rdf','all']):
            coord_mod2 = CoordinationAnalysisModifier(cutoff = 8, number_of_bins = 1000,partial=0)
            odata_cp.apply(coord_mod2) #We need to apply this seperately as it is weirdly not applied to the first entry/structure.
            if key in rdfs: rdfs[key].append(odata_cp.tables['coordination-rdf'].xy())
            else: rdfs[key]=[odata_cp.tables['coordination-rdf'].xy()]

        #Compute the shape related parameters (anisotrpy, radius of gyration etc.
        #if args.orderAnalysis:
        if args.type.intersection(['s','shape','all']):

            ctable = odata.tables['clusters'] #we assume to be only one cluster
            #print('ctable: ',ctable)
            Rg=ctable['Radius of Gyration'][0]
            COM=ctable['Center of Mass'][0]
            gtensor=ctable['Gyration Tensor'][0]
            #cluster_sizes = odata.tables['clusters'].y; noClus=len(cluster_sizes)

            #print (inpf, atoms,dens, 'Rg: ',Rg,'COM: ', COM,'COM from ASE: ',atoms.get_center_of_mass(),'GTensor: ',gtensor)

            gtensMat=np.zeros((3,3))
            gtensMat[0,0]=gtensor[0] #XX
            gtensMat[1,1]=gtensor[1] #YY
            gtensMat[2,2]=gtensor[2] #ZZ
            gtensMat[0,1]=gtensor[3] #XY
            gtensMat[1,0]=gtensor[3] #XY
            gtensMat[2,0]=gtensor[4] #XZ
            gtensMat[0,2]=gtensor[4] #XZ
            gtensMat[1,2]=gtensor[5] #YZ
            gtensMat[2,1]=gtensor[5] #YZ

            #Get egienvalues of gyration tensor
            evals=sorted(np.linalg.eigvals(gtensMat))

            lx=evals[0];ly=evals[1];lz=evals[2]
            #print (gtensMat,lx,ly,lz)

            #Following LAMMPS definitions
            c=lz-0.5*(ly+lx) #acylindiricty
            b=ly-lx #asphericity
            k=1.5*(lx**2+ly**2+lz**2)/(lx+ly+lz)**2-0.5 #relative shape aniostropy

            #print ('eigvals: ', evals)
            #print (inpf,'Rg: ',Rg,'b,c,k :',b,c,k)

            shapes.append([cSize,dens,Rg, b,c,k])#asp,acyl,aniso
            #if key in shape: shape[key].append([Rg,b,c,k])
            #else: shape[key]=[Rg,b,c,k]

        if args.type.intersection(['r','ring','rings','all']):  #TODO: move this part to after filtering based on structure energy below (to avoid do costly RIGNS calcs for unncessary geoms).
            fname=inpf
            odir='rings-out'
            wdir='%s/%s/'%(odir,seed)
            system('mkdir -p %s/'%(wdir))
            os.chdir('%s'%(wdir))
            #print('cwd: ',os.getcwd())
            if not (os.path.exists('./data')):  system('ln -s ../../data .')
            if not os.path.exists('data/%s.xyz'%seed):  system("ase-file-converter-bk.py -i ../../%s.res -ot xyz -o data/%s  &> /dev/null "%(seed,seed)) # &>/dev/null  #Now done altogether in the beginning!!
            
            syms=atoms.get_chemical_symbols()
            atypes,x=np.unique(syms,return_counts=1)
            acounts=dict(zip(atypes,x))
            
            #Prepare the R.I.N.G.S code input for the given input structure file. #TODO:use a spearate working dir for each RINGS job, and use only one my-input and options files, pnly options is suitable!!  
            open('my-input','w').write("""#######################################
#       R.I.N.G.S. input file         #
#######################################
Title                    #
%d                                   # no of atoms
1                                    # no of elements
%s                                   #  atom types
1                                    # no of structures
0                                     # #use lattice paramter type
15.000000 15.00000 15.00000           # Cell params won't matter as we are interested in non-PBC
90.000000 90.00000 90.00000           #
2.5                                   # MD time step
ANI                                   # xyz format
%s.xyz               # input file
200                                   # Real space discretization for the g(r) calculations.
500                                   #Reciprocal space discretization for the S(q) calculations.
25                                    # Maximum modulus of the reciprocal space vectors for the S(q) calculations.
0.125                                 #Smoothing factor for the S(q) calculated using the sampling of vectors from the reciprocal space.
90                                    # Angular discretization between [0-180] degrees used to compute the angle distribution.
20                                    # Real space discretization for the voids and for the ring statistics.
10                                    # (Maximum search depth)/2 for the ring statistics - (min = 2) def: 10
15                                    # Maximum search depth for the chain statistics - (min = 2)
#######################################
C C %.2f                              #Cutoff radius from the partials g(r) *
Grtot %.2f                            #Cutoff radius from the total g(r)
"""%(len(atoms),' '.join(atypes), seed, args.cutoff,args.cutoff))

            
            open('options','w').write("""#######################################
        R.I.N.G.S. options file       #
#######################################
 PBC             .false.               # switch off for clusters
 Frac            .false.              # Fractional (.true.) or real (.false.) atomic coordinates.
 g(r)            .false.              # Evaluation of the radial distribution functions,
 S(q)            .false.              # Evaluation of the neutron and X-rays structure factors, S(q), from FFT of the g(r):
 S(k)            .false.              # Evaluation of the neutron and X-rays structure factors, S(q), from single atomic correlations (Debye equation):
 gfft(r)         .false.              # Evaluation of the radial distribution functions, g(r) from the Fourier transform of the structure factor calculated using the Debye equation
 MSD             .false.              #
 atMSD           .false.              #
 Bonds           .true.               #
 Angles          .true.              #
 Chains          .true.               #
---- ! Chain statistics options ! -----
 Species           0                  # 0 = All atom types are used
 AAAA            .false.              # Only look for AAAA chains.
 ABAB            .false.              # Only look for ABAB chains, homopolar bonds are taken into account when creating the nearest neighbors table.
 1221            .false.              # Only look for 1-(2)n-1 coordinated atoms chains, ie. isolated chains.
---------------------------------------
 Rings           .true.               #
---- ! Ring statistics options ! -----
 Species           0                  #
 ABAB            .false.              # Only look for ABAB rings. Homopolar bonds are ignored when creating the nearest neighbors table. Only relevant for multi-element systems.
 Rings0          .false.              # Look for all closed paths in the box. #Does not work together with others!!
 Rings1          .false.              #  Look for King's shortest path rings - no homopolar(A-A)  bonds
 Rings2          .false.              #  Look for Guttman's shortest path rings - no homoplar
 Rings3          .true.               # Look for King's shortest path rings - with homoplar bonds
 Rings4          .true.               # Look for Guttman's shortest path rings - with homoplar
 Prim_Rings      .true.               #  look for primitive rings (irreduible)
 Str_Rings       .false.              # this analysis is problemetaic for amorphous structures and reuqires large mem.
 BarycRings      .false.              # Positions of the barycenter of the rings.
 Prop-1          .true.              # Analysis of standard properties for atoms involved in rings
 Prop-2          .false.              # Advanced properties of the barycenter of the rings (RDF) and neutron
 Prop-3          .false.              # Neutron structure factor of particles at the origin of rings with n nodes
 Prop-4          .false.              # Advanced properties of the atoms involved in rings 2
 Prop-5          .false.              #  Average neutron structure factors of a ring with n nodes
---------------------------------------
 Voids           .false.              #
#######################################
         Outputting options           #
#######################################
 Evol            .true.              # Evaluation and output of g(r), s(q), ... for each configuration
 Dxout           .false.              #
! OpenDX visualization options !  --
 RadOut          .false.              #
 RingsOut        .false.              #
 DRngOut         .false.              #
 VoidsOut        .false.              #
 TetraOut        .false.              #
 TrajOut         .false.              #
------------------------------------
 Output        %s.ringsout          #
#######################################
# 0 => All closed paths (no rules)
#1 => King's shortest path rings without homopolar bonds
#2 => Guttman's shortest path rings without homopolar bonds
#3 => King's shortest path rings with homopolar bonds
#4 => Guttman's shortest path rings with homopolar bonds
#5 => Primitive rings
#6 => Strong rings
"""%(seed))

            #cmd='rings my-input &> %s/%s.ringsout'%(odir,seed)
            cmd='rings my-input &> ../%s.ringsout'%(seed)
            if not os.path.exists('../%s-rings-3.dat'%(seed)) : #& os.path.exists('%s/%s-rings-3.dat'%(odir,seed)) & os.path.exists('%s/%s-rings-4.dat'%(odir,seed)) & os.path.exists('%s/%s-rings-5.dat'%(odir,seed))): #(os.path.exists('%s.ringsout'%(seed)))
            #if not os.path.exists('%s/%s{.ringsout,-rings-3.dat,-rings-4.dat,-rings-5.dat}'%(odir,seed)):
                #system(cmd);time.sleep(1.0)           

                """
                cnt2=0;flag=0
                while cnt2<3 and not flag: #try three times, each rigns calc shouldn't take more than 1 min.
                    try:
                        x=run('rings my-input',stdout=DEVNULL,stderr=DEVNULL,check=0,shell=1,timeout=40).returncode #capture_output=0  #Waits until the process is finished and returns 0 when finished w/out errors.  #WORKS BEST!!
                        if x==0:flag=1;break
                    except: cnt2+=1;continue
                """
                try:x=run('rings my-input',stdout=DEVNULL,stderr=DEVNULL,check=0,shell=1,timeout=60).returncode
                except:
                    try:run('rings my-input',stdout=DEVNULL,stderr=DEVNULL,check=0,shell=1,timeout=60).returncode
                    except:
                        try:run('rings my-input',stdout=DEVNULL,stderr=DEVNULL,check=0,shell=1,timeout=60).returncode
                        except:None
                
                #if args.rType=='Kings':system('mv rstat/evol-RINGS-3.dat %s/%s-rings-3.dat'%(odir,seed))#Kings
                #elif args.rType=='Guttman':system('mv rstat/evol-RINGS-4.dat %s/%s-rings-4.dat'%(odir,seed))#Guttman
                #elif args.rType=='primitive':system('mv rstat/evol-RINGS-5.dat %s/%s-rings-5.dat'%(odir,seed)) #Primitve #RINGS-res-5.dat
                #if args.rType=='strong': system('mv rstat/evol-RINGS-6.dat %s/%s-rings-6.dat'%(odir,seed)) #strong
                #if x==0:
                if os.path.exists('rstat/evol-RINGS-3.dat'):
                    #if os.path.exists('%s-rings-3.dat'%(seed)) :
                    system('mv rstat/evol-RINGS-3.dat ../%s-rings-3.dat'%(seed))#Kings
                    system('mv rstat/evol-RINGS-4.dat ../%s-rings-4.dat'%(seed))#Guttman
                    system('mv rstat/evol-RINGS-5.dat ../%s-rings-5.dat'%(seed)) #Primitve #RINGS-res-5.dat
                    #system('mv rstat/evol-RINGS-6.dat %s-rings-6.dat'%(seed)) #strong
                    system('mv cstat/chains.dat ../%s-chains.dat'%(seed))
                else: 
                    print('Problem with RINGS calculation of %s.res, skipping...'%seed);os.chdir('../../'); return None;
                        

            #Read the cstat and rstat outputs.
            #if x==0:#if the rings run finished witout errors
            try:
                chain=np.loadtxt('../%s-chains.dat'%(seed),unpack=1)
                if key in chains:chains[key].append(chain)
                else:chains[key]=[chain]
            except: None


            try:
            #if 1:
                if args.rType in ['Kings']:ring=np.loadtxt('../%s-rings-3.dat'%(seed),unpack=1,skiprows=5) #King's shortest path algo
                elif args.rType in ['Guttman']:ring=np.loadtxt('../%s-rings-4.dat'%(seed),unpack=1,skiprows=5) #Guttman's shortest path algo
                elif args.rType in ['primitive']:ring=np.loadtxt('../%s-rings-5.dat'%(seed),unpack=1,skiprows=5) #Primitve ring algo
                #elif args.rType in ['strong']:ring=np.loadtxt('../%s-rings-6.dat'%(seed),unpack=1,skiprows=5) #Primitve ring algo
            except: print('Problem with RINGS calculation of %s.res, skipping...'%seed);os.chdir('../../'); return None; #continue
                
            if key in rings:rings[key].append(ring)
            else:rings[key]=[ring]

            if 1:
                #system('rm -rf angles/ bonds/ rstat/ cstat/ tmp/ lirr*.dat lred*.dat') #clean up the RINGS output folders.  #; mv my-input options %s'%odir)
                os.chdir('../')
                system('rm -rf %s'%seed)
                os.chdir('../')
            else:os.chdir('../../')
            

        if args.type.intersection(['ba','bondangle','all']):
            bond_mod1=CreateBondsModifier(cutoff = args.cutoff)
            bond_mod2=BondAnalysisModifier(bins=360,cosine_mode=0,length_cutoff=args.cutoff*1.2)
            odata.apply(bond_mod1);   odata.apply(bond_mod2)

            ba_hist = odata.tables['bond-angle-distr'].xy().T  # orig size: #bins x 2
            bl_hist = odata.tables['bond-length-distr'].xy().T
            #print(ba_hist,bl_hist)
            if key in blens: blens[key].append(bl_hist) # blens size: #str x 2 x #bins
            else: blens[key]=[bl_hist]
            if key in bangs: bangs[key].append(ba_hist)
            else: bangs[key]=[ba_hist]
            

        return btypes,btdict,shapes,rdfs,rings,chains,toDel, saveList,blens,bangs#maxs,#cnt,
  

    strind=0
    for res in p_imap(sortData,[[inpf,i,allEs[i]] for i,inpf in enumerate(inpfs)],num_cpus=args.nprocs):  #,mins,maxs
        #btypes,btdict,shapes,rdfs,rings,chains,toDel, saveList
        if res != None:
            if len(res)==2:
                toDel.extend(res[0])
                saveList.extend(res[1])
            else:
                btypes.extend(res[0])
                shapes.extend(res[2])
                toDel.extend(res[6])
                saveList.extend(res[7])
                for key in res[1].keys():btdict[key].extend(res[1][key]) #dict  
                for key in res[3].keys():rdfs[key].extend(res[3][key]) #dict   
                for key in res[4].keys():
                    if key in rings: rings[key].extend(res[4][key]) #dict
                    else:rings[key]=res[4][key]
                for key in res[5].keys():chains[key].extend(res[5][key]) #dict
                #for key in res[8].keys():maxs[key]=res[8]
                for key in res[8].keys(): #blens[key]=res[8]
                    if key in blens:blens[key].extend(res[8][key])
                    else:blens[key]=res[8][key]
                for key in res[9].keys(): #bangs[key]=res[9]
                    if key in bangs:bangs[key].extend(res[9][key])
                    else:bangs[key]=res[9][key]

        
        #strind+=1
        
    #print(btypes,btdict,shapes,rdfs,rings,chains,toDel, saveList)
    print("Elapsed time: %.2f sec."%( time.time()-initT))    
          
    # End of input structure loop #
    ##########


    ################
    #Do some preprations for plotting
                
    gmin=np.min(list(mins.values()))
    gmax=np.max(list(maxs.values()))
    print('Global E_min= %.5f eV/atom E_max=%.5f eV/atom'%(gmin,gmax))
    
    btypes=np.array(btypes)
    shapes=np.array(shapes)

    #NowDelete the structures above the ecut threshold (higher neergy)
    if len(toDel)>0:
        print('%d structures above the cutoff energy (%.2f eV/atom) are ignored.'%(len(toDel),args.ecut))

        
    if args.save:
        print('%d structures are saved in %s/'%(len(saveList),args.sfold));exit()
        #print('%d structures are saved in %s/'%(len(saveList),args.sfold));system('cp %s %s/ 2> /dev/null '%(' '.join(saveList),args.sfold)); exit()
    
    #if args.type.intersection(['r','ring','rings','all']): system('rm -rf angles/ bonds/ rstat/ cstat/ tmp/ lirr*.dat lred*.dat; mv my-input options %s'%odir) #clean up the RINGS output folders.
    
    alphas=[];#toDel=[]
    for i,bt in enumerate(btypes):
        key=str('%d'%bt[0]);E=bt[-2]
        if args.globalCM:
            mins[key]=gmin;maxs[key]=gmax #set this to use global values for the colormap; otherwise the structures are coloured wrt the relative energy to the local minimum for each cluster size.
            alphas.append(E) #Use the total E/atom directly in colouring.
        else:
            if maxs[key]==mins[key]:x=1.0
            else:x=1-(E-mins[key])/(maxs[key]-mins[key])*0.90 #FIX: max points are not shown, having a alpha=0 value, so we use a scaling of 0.95 to set alpha values to 0.05-1.0
            alphas.append(x) #normalised relative E/atom for the given cluster size

    alphas=np.array(alphas)
    #print('my alphas', alphas)
        
    btypes=btypes.T
    try: shapes=shapes.T
    except:None
 
    # End of preparations part
    ###########################

    
    #########
    # Plotting part
    #########
    print('Plotting the data from the final set of %d structures...'%(len(btypes.T)))
    font = {'family': 'serif', 'size': 18}
    plt.rc('font', **font)

    if args.type.intersection(['c','coord','all']):
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))
        #fig, ax = plt.subplots(1, 1)
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx()    

        #btype[1]:sp3, btype[2]:sp2, btype[3]:sp, btype[4]:others
        sc=ax1.scatter(btypes[0],btypes[1],label='sp3',marker='o',lw=0,c=alphas, cmap=mycmap('C0',args.globalCM))
        sc=ax1.scatter(btypes[0],btypes[2],label='sp2',marker='o',lw=0,c=alphas, cmap=mycmap('C1',args.globalCM))
        #sc=ax1.scatter(btypes[0],btypes[3],label='sp',marker='o',lw=0,c=alphas, cmap=mycmap('C2',args.globalCM))

        ax1.set_ylim(0,1)
        ax1.set_ylabel('Fraction')
        ax2.set_yticks([])
        ax1.set_xlabel('Cluster size')
        #fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.  #Fix the color label (all blue!)

        if args.globalCM:
            #Colorbar is  only meaningful when using the global values (not individual ones for each cluster size)
            step=0.5
            ax2.set_ylim(gmin,gmax)
            cbar=plt.colorbar(sc)#,ticks=[x for x in np.arange(gmin,gmax+step,step)])
            cbar.set_ticks(np.arange(gmin,gmax,step))
            cbar.ax.set_yticklabels(['%.1f'%x for x in np.arange(gmin,gmax,step)])

            cbar.set_label('Energy [eV/atom]', rotation=90)
            cbar.ax.set_ylim(gmin,gmax)
            #plt.ylim(gmin,gmax)

        plt.title(cwd_label)
        if not args.noplot:         plt.show(block=0)

        ####
        # Now plot the Boltzmann-Weighted (average) fraction vs. cluster size.
        ####
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))

        #First sort the data.
        abtypes={}
        keys=sorted([int(x) for x in energy.keys()])
        for key in keys:#loop over clustr sizes
        #for key in tqmd(keys):#loop over clustr sizes
            key=str(key)
            BW=np.array(boltz_dist(energy[key],T=args.temp)) #Boltzmann weight for each structure in a given size
            abtypes[key]={} #average bond types info for all structures of a given size
            noStr=0
            BWs=[]
            for sind,btype in enumerate(btdict[key]):#loop over structures
                btype=np.array(btype)
                #if np.isclose(BW[sind],0.0,atol=1e-4):continue #This was added to speed up a bit, by skipping high-energy structures, but causes errors with keys.          
                BWs.append(BW[sind])
                #btype[0]: cluster size, btype[0]:sp3, btype[1]:sp2, btype[2]:sp, btype[3]:others
                if  len(abtypes[key].keys())==0:#first str
                    abtypes[key]['sp3']=[btype[0]]
                    abtypes[key]['sp2']=[btype[1]]
                    abtypes[key]['sp']=[btype[2]]
                    abtypes[key]['avg']=[btype[-1]]
                    abtypes[key]['dens']=[btype[-2]]
                else:
                    abtypes[key]['sp3'].append(btype[0])
                    abtypes[key]['sp2'].append(btype[1])
                    abtypes[key]['sp'].append(btype[2])
                    abtypes[key]['dens'].append(btype[-2])                    
                    abtypes[key]['avg'].append(btype[-1])
                
            #print ('BWs:', BWs)

            keys2=abtypes[key].keys()
            for key2 in keys2: #sp,sp2,sp3,others
                abtypes[key][key2]=np.average(abtypes[key][key2],weights=BWs,axis=-1) #TODO: Add std. dev.
            #end of cluster size loop
            
        #print('Average btypes',abtypes)
        
        #r2_vals=[];p_vals=[];stderrs=[]
        mc=[]
        for i,hybrid in enumerate(['sp3','sp2']): #,'sp']):
            #try:
            data=np.array([[float(key), float(abtypes[key][hybrid])] for key in abtypes.keys()]).T
            plt.plot(data[0],data[1],label=hybrid,marker='o',ms=4,lw=1,color='C%d'%i)
            #except:   print('Hata: ',i,key,abtypes[key])
            


            #Draw the best-fit line (options: least-squares,  linear regression or polynomials).
            x=data[0][25:];y=data[1][25:]; #?? Skip some of the first values to get a better fit ??
            if len(x)==0 or len(y)==0: x=data[0];y=data[1]
            m,c=linear_fit(x,y)              
            startPt=[x[0],x[-1]]
            endPt=[m*x[0]+c,m*x[-1]+c]
            label='y=%.1e*x+%.2f'%(m,c)
            label=None
            plt.plot(startPt, endPt, lw=3,color='C%d'%i,label=label) #plot the best fit line
            #plt.plot(x, m*x + c, 'r', label='Fitted line')
            mc.append([m,c])
            
        #Determine the intersection point by interpolation/extrapolation
        ix=(mc[1][1]-mc[0][1])/(mc[0][0]-mc[1][0]) #x1=x2=(c2-c1)/(m1-m2)

        print('Slope of the sp2 bestfit line: %.3e'%mc[1][0])
        print('Slope of the sp3 bestfit line: %.3e'%mc[0][0])
        print('sp2 and sp3 best-fit lines intersect at x=%.1f, y=%.1f '%(ix,mc[0][0]*ix+mc[0][1]))
        plt.plot(x[0],y[0],ms=0,lw=0,label=r'$C_{int}$=%d'%(np.around(ix)))
        #fig.legend(ncol=3,loc='upper center') #fontsize='x-small'
        plt.legend(ncol=3,fontsize='small') #loc='upper right', fontsize='x-small'
        plt.xlabel('Cluster Size')
        plt.ylabel('Average Fraction')
        plt.ylim(0,1)
        
        plt.title(cwd_label)
        if not args.noplot:  plt.show(block=0)

        ####
        # Now plot the Boltzmann-Weighted (average) fraction vs. cluster size.
        ####
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))        
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx()    

        for i,hybrid in enumerate(['avg']): #,'sp']):
            data=np.array([[float(key), float(abtypes[key][hybrid])] for key in abtypes.keys()]).T
            ax1.plot(data[0],data[1],label=hybrid,marker='o',ms=4,lw=1,color='C0') #color='C%d'%i)
            np.savetxt('data-mean-coord-no.dat',data.T,header='Cluster Size vs. Average Coord No')

        for i,hybrid in enumerate(['dens']): #,'sp']):
            data=np.array([[float(key), float(abtypes[key][hybrid])] for key in abtypes.keys()]).T
            ax2.plot(data[0],data[1],label=hybrid,marker='o',ms=4,lw=1,color='C1')#color='C%d'%i)

        ax1.set_xlabel('Cluster Size')
        ax1.set_ylabel('Average Coordination Number')
        ax2.set_ylabel('Density [g/cc]')
        #ax2.set_ylim(0,10)
        
        ax2.annotate(r'$R_{0}$ = %.2f A'%args.cutoff,(0.45,0.95),fontsize='small',xycoords='axes fraction')
        
        #ax1.spines['right'].set_color('C0')
        ax1.yaxis.label.set_color('C0')
        ax1.tick_params(axis='y', colors='C0')
        
        #ax2.spines['left'].set_color('C1')
        ax2.yaxis.label.set_color('C1')
        ax2.tick_params(axis='y', colors='C1')



        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)


        #####
        # Now print out the relevant data files
        #####

        #np.savetxt('data-mean-coord-C%s.dat'%key,np.array([[int(key),ab['sp3'],ab['sp2'],ab['sp']] for ab in abtypes[key]]), header='#ClusterSize sp3 sp2 sp ratios')

        with open('data-mean-coord.dat','w') as f:
            f.write('#Slope of the sp2 bestfit line: %.3e\n'%mc[1][0])
            f.write('#Slope of the sp3 bestfit line: %.3e\n'%mc[0][0])
            f.write('#sp2 and sp3 best-fit lines intersect at x=%.1f, y=%.1f \n'%(ix,mc[0][0]*ix+mc[0][1]))
            f.write('#cSize sp3 sp2 sp ratios and average coord number per atom\n')

            for key in abtypes.keys():
                ab=abtypes[key]
                f.write('%6s %5.3f %5.3f %5.3f %5.3f\n'%(key,ab['sp3'],ab['sp2'],ab['sp'],ab['avg']))

    ####
    # Plot the BW bondd length and bond angle distributions
    ####

    if args.type.intersection(['ba','bondangle','all']):
        if not args.noplot:
            f1 = plt.figure(figsize=(9.,6.))
            f2 = plt.figure(figsize=(9.,6.))
            ax1 = f1.add_subplot(111)
            ax2 = f2.add_subplot(111)

        #Sort the data first
        keys=sorted([int(x) for x in energy.keys()])
        for key in keys: #loop over cluster sizes
            key=str(key)
            BW=np.array(boltz_dist(energy[key],T=args.temp)) #Boltzmann weight for each structure in a given size
            blen=np.array(blens[key])#.T #original blens: (#bins x 2 x #structures)
            bang=np.array(bangs[key])#.T

            sh=np.shape(blen[-1])
            ablen=np.zeros((sh))#.T #avg rdf
            abang=np.zeros((sh))
            noStr=len(blen)
            #ablen[0]=blen[-1][0] #x: distances, y: no of counts
            #ablen[1]=np.average(blen[:][1],weights=BW,axis=0)
            sm=np.sum(BW) 
            #"""
            for i,r in enumerate(blen): #loop over structures
                #r=r.T
                b=bang[i]
                if i==0:ablen[0]+=r[0];abang[0]+=b[0]
                if not np.isclose(BW[i],0.0,atol=1e-4):
                    ablen[1]+=r[1]*BW[i]/sm 
                    abang[1]+=b[1]*BW[i]/sm

                if not args.noplot and 0:#Plot individual contribution of each structure
                    ax1.plot(r[0],r[1],label='C%s-%s-%.1f'%(key,i,BW[i]),lw='%.3f'%BW[i])
                    ax2.plot(b[0],b[1],label='C%s-%s-%.1f'%(key,i,BW[i]),lw='%.3f'%BW[i])

            #No need for these as we divide by the sum of weights.
            #ablen[1]/=float(len(blen))
            #abang[1]/=float(len(abang))          
   
            if 0:  ablen[1]/=np.max(ablen[1]) ; abang[1]/=np.max(abang[1])  #normalise
            
            if not args.noplot:
                ax1.plot(ablen[0],ablen[1],label='Mean-C%s'%(key),lw=2) 
                ax1.set_xlabel('Bond Length [A]')
                ax1.set_ylabel('Count')
                ax1.set_xlim(1,args.cutoff*1.2)

                ax2.plot(abang[0],abang[1],label='Mean-C%s'%(key),lw=2) 
                #plt.legend(ncol=3,fontsize='x-small') #loc='upper center'
                ax2.set_xlabel('Bond angle [deg]')
                ax2.set_ylabel('Count')
                ax2.set_xlim(0,180)

            #WARNING there may be too many labels !!!
            if not args.noplot and len(keys)<10:
                f1.legend(ncol=3,fontsize='x-small') #loc='upper center'
                f2.legend(ncol=3,fontsize='x-small') #loc='upper center'

            np.savetxt('data-mean-blengths-C%s.dat'%key,ablen.T)
            np.savetxt('data-mean-bangles-C%s.dat'%key,abang.T)
            
        if not args.noplot:
            plt.title(cwd_label)
            plt.show(block=0)        

    ####
    # Plot the Energy/atom vs. cluster size
    ####
    if args.type.intersection(['e','energy','all']):

        if not args.noplot: fig = plt.figure(figsize=(9.,6.))
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx()      

        sc=ax1.scatter(btypes[0],btypes[-2],label='Energy',marker='o',lw=0,c=alphas, cmap=mycmap('b',args.globalCM))#alpha=alphas) #s=10**2 #ec=None

        ax1.set_ylabel('Energy [eV/atom]')
        ax2.set_yticks([])
        #ax2.axes.get_yaxis().set_visible(False)
        #plt.ylabel('Fraction')
        ax1.set_xlabel('Cluster size')
        #fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.  #Fix the color label (all blue!)

        np.savetxt('data-energy.dat',[[btypes[0][i],btypes[-2][i]] for i in range(len(btypes[0]))],header='cSize, Energy [eV/atom]')
        
        if args.globalCM:
            #Colorbar is  only meaningful when using the global values (not individual ones for each cluster size)
            step=0.5
            ax2.set_ylim(gmin,gmax)
            cbar=plt.colorbar(sc)#,ticks=[x for x in np.arange(gmin,gmax+step,step)])
            cbar.set_ticks(np.arange(gmin,gmax,step))
            cbar.ax.set_yticklabels(['%.1f'%x for x in np.arange(gmin,gmax,step)])

            cbar.set_label('Energy [eV/atom]', rotation=90)
            cbar.ax.set_ylim(gmin,gmax)
            #plt.ylim(gmin,gmax)

        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)

    #exit()
    ###
    #Plot density of energies
    ###
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx()      

        #hist,bins,_=ax1.hist(btypes[-2],bins=200,density=1,facecolor='b',alpha=0.75,histtype='stepfilled',stacked=1); 
        hist,bins=np.histogram(btypes[-2],bins=200,density=1) #normed=1) #density=0,
        density = stats.gaussian_kde(btypes[-2])
        dens=[density(bn) for bn in bins]
        #hist/=np.sum(hist*np.diff(bins))
        #dens=[density(bn)/np.sum(hist*np.diff(bins)) for bn in bins]
        #sm=np.sum(hist*np.diff(bins))
        #print(hist)
        #hist/=sm
        #hist=[hs/np.sum(hist[i]*(bins[i+1]-bins[i]))) for i,hs in enumerate(hist) ]
        #hist/=np.max(hist);  dens/=np.max(hist)
        
        if 0:
            ax2.bar(bins[0:-1],hist,width=0.025)
            ax1.plot(bins,dens,label='C%d'%int(key))
        else:
            ax=sns.distplot(btypes[-2], hist=True, kde=True, 
                            bins=200, color = 'darkblue', 
                            hist_kws={'edgecolor':'black'},
                            kde_kws={'linewidth': 4},
            )
            pickle.dump(ax, open("data-energy-dens.pickle", "wb"))
        
        ax1.set_xlabel('Energy [eV/atom]')
        ax2.set_yticks([])
        ax2.axes.get_yaxis().set_visible(False)
        ax1.set_ylabel('Density')
        
        #plt.legend(ncol=2,loc='upper right',fontsize='medium')

        np.savetxt('data-energies.dat',btypes[-2],header='Energies')
        np.savetxt('data-energy-dens.dat',[[bins[i],hist[i],density(bins[i])] for i in range(len(hist))],header='Energy, #hits, density')
        
        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)

        
    ###
    #Plot structural dens of states
    ###

    if args.type.intersection(['d','dos','all']):
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))

        keys=[];vals=[]
        for key in energy.keys():
            for val in energy[key]:
                keys.append(int(key))
                vals.append(val)


        #plt.hist(vals,bins=50,histtype='stepfilled',density=1,color=['steelblue' for x in vals], edgecolor='none')#dens means normalise
        plt.hist2d(keys,vals,bins=50,density=1,cmap='Blues')#cmap='inferno'#cmap='YlOrBr')#cmap='Wistia',cmap=plt.cm.jet, #norm=mcolors.PowerNorm(gamma))


        if 0: #Apply Gaussian smearing
            # TODO: the x and y range is in terms of number of bins, convert to real ranges !!
            from scipy.ndimage import gaussian_filter
            from matplotlib import cm

            data = np.histogram2d(keys, vals, bins=150)[0]
            data = gaussian_filter(data, sigma=2)#orig:5
            #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
            norm = cm.colors.Normalize(vmax=1, vmin=0)
            norm = cm.colors.Normalize()
            plt.pcolormesh(data.T, cmap='Blues', shading='gouraud',norm=norm) #cmap='inferno'  #'Blues_r' is the reversed order!!
            #plt.hist2d(data.T[0],data.T[1],bins=50,density=0,cmap='Blues')#cmap='inferno'#cmap='YlOrBr')#cmap='Wistia',cmap=plt.cm.jet, #norm=mcolors.PowerNorm(gamma))  #DOES NOT WORK
            fig.canvas.draw()

        plt.xlabel('Cluster size')
        plt.ylabel('Energy [eV/atom]')

        plt.grid(True)
        cb=plt.colorbar()
        cb.set_label("density", rotation=90)
        #cb.ax.set_ylim(0,1)

        plt.title(cwd_label)
        if not args.noplot:         plt.show(block=0)

        #exit()

    
    ###
    # Plot Radius of gyration and anisoptrpy vs cluster size
    ###
    #if args.orderAnalysis:
    if args.type.intersection(['s','shape','all']):
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))
        #fig, ax = plt.subplots(1, 1)
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx()

        #[cSize,dens,Rg, b,c,k])#asp,acyl,aniso
        sc=ax1.scatter(shapes[0],shapes[2],label='Gyration Radius',marker='o',lw=0,c=alphas, cmap=mycmap('C0',args.globalCM))
        sc=ax2.scatter(shapes[0],shapes[5],label='Shape anisotropy',marker='o',lw=0,c=alphas, cmap=mycmap('C1',args.globalCM))


        ax1.set_ylabel('Radius of Gyration [A]',color='C0')
        ax2.set_ylabel('Anisotropy',color='C1')
        ax2.set_ylim(0,1)
        ax1.set_xlabel('Cluster size')
        
        #ax1.spines['right'].set_color('C0')
        ax1.yaxis.label.set_color('C0')
        ax1.tick_params(axis='y', colors='C0')
        
        #ax2.spines['left'].set_color('C1')
        ax2.yaxis.label.set_color('C1')
        ax2.tick_params(axis='y', colors='C1')
        
        #ax1.legend(ncol=2,loc='upper center',fontsize='small') #this shows all legends from all subplots combined.  #Fix the color label (all blue!)

        #Print the values into files
        np.savetxt('data-mean-shapes.dat',[[x[0], x[2],x[5]] for x in shapes.T],header='Cluster size, Gyration Radius (R_g), Shape anisotropy (k)')#,fmt=['%-5d','%5.3f','%5.3f'])
            

        if args.globalCM:
            #Colorbar is  only meaningful when using the global values (not individual ones for each cluster size)
            step=0.5
            ax2.set_ylim(gmin,gmax)
            cbar=plt.colorbar(sc)#,ticks=[x for x in np.arange(gmin,gmax+step,step)])
            cbar.set_ticks(np.arange(gmin,gmax,step))
            cbar.ax.set_yticklabels(['%.1f'%x for x in np.arange(gmin,gmax,step)])

            cbar.set_label('Energy [eV/atom]', rotation=90)
            cbar.ax.set_ylim(gmin,gmax)
            #plt.ylim(gmin,gmax)

        plt.title(cwd_label)
        if not args.noplot:         plt.show(block=0)
        #exit()


        ###
        #Order paramter plot: sp2 vs. shape anisotropy (k)
        ###
        #if args.type.intersection(['o','order','all']):

        if not args.noplot: fig = plt.figure(figsize=(9.,6.))

        plt.hist2d(btypes[2],shapes[5],150,density=0,cmap='YlOrBr')#cmap='Blues'#cmap='Wistia',cmap=plt.cm.jet, #norm=mcolors.PowerNorm(gamma))
        #btypes[1]:sp3, btypes[2]:sp2, btypes[3]:sp, btypes[4]:others
        plt.grid(True)
        plt.colorbar()

        plt.ylim(0,1)
        plt.xlim(0,1)
        plt.xlabel('sp2 fraction')
        plt.ylabel('Anisotropy')

        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)

        #exit()

    ###
    # Plot Rings statistics #
    ###
    if args.type.intersection(['r','ring','rings','all']):

    #plot Boltzmann-weighted/avg no of rings for each ring size (avg. over cluster size/no of rings??) over all structures.
    
        if len(clist)<=12 and not args.noplot: fig = plt.figure(figsize=(9.,6.))
        
        #Sort the data first
        allRings={}#stores the all ring size info for all structures for all cluster sizes.
        allChains={}
        keys=sorted([int(x) for x in energy.keys()])
        #print(keys)
        for key in keys:#loop over clustr sizes
        #for key in tqmd(keys):#loop over clustr sizes
            key=str(key)
            BW=np.array(boltz_dist(energy[key],T=args.temp)) #Boltzmann weight for each structure in a given size
            crings={} #list of rings of different sizes for all structures of a given size
            #cchains={}
            noStr=0
            BWs=[]
            for sind,ring in enumerate(rings[key]):#loop over structures
                totR=np.sum(ring[:][1],axis=0) #Only consider R3-R12 or all upto R20??
                #print(key,sind,ring[:][1],totR)
                ring=ring.T
                #if np.isclose(BW[sind],0.0,atol=1e-4):continue
                BWs.append(BW[sind])
                #Using the first value (Rc*noAtoms=noRings) from
                # n           Rc            PN           Pmax          Pmin

                for rind,r in enumerate(ring):#loop over ring sizes
                    key2=r[0]#ring size
                    #val=np.around(r[1])/float(key) #no of ring per atom (Rc) #*BW[sind] #now done below
                    if totR!=0:val=r[1]/totR #fraction of ring  #*BW[sind] #now done below
                    else:val=0.
                    if key2 in crings:  crings[key2].append(val)  #+=val
                    else:crings[key2]=[val]#0.0#
                noStr+=1
            #print ('End crings:',crings,'BWs:',BWs)
                
            keys2=crings.keys()
            for key2 in keys2:#loop over ring sizes
                crings[key2]=np.average(crings[key2],weights=BWs,axis=-1) #TODO: Add std. dev and stderr.??

            if len(clist)<=12: plt.plot([x for x in crings.keys()],[y for y in crings.values()],label='C%s'%(key),ms=4,marker='o',lw=1) 
            allRings[key]=crings

            #TODO: Take average of chain no over available structures
            
            if 1:
                #print('Writing the Boltz-Weighted rings count per atom data for each cluster sizes in data-mean-rings-Cx.dat...')
                np.savetxt('data-mean-rings-C%s.dat'%key,[[x, crings[x]] for x in keys2],header='Ring size vs average ring fraction and average chain no',fmt=['%-5d','%5.3f'])
                #np.savetxt('data-mean-chains-C%s.dat'%key,[,header='Ring size vs average ring count per atom',fmt=['%-5d','%5.3f'])
            
            #End of cluster size loop

        #plt.bar(aring[0],aring[1],width=0.05)
        if len(clist)<=12: 
            plt.legend(ncol=3,fontsize='x-small') #loc='upper center'
            plt.xlabel('Ring Size (R_x)')
            #plt.ylabel('Average Ring Count per Atom') #  rings per atom (Rc)
            plt.ylabel('Average Ring Fraction') #  rings per atom (Rc)
            plt.xlim(3,12)
        
            plt.title(cwd_label)
            if not args.noplot: plt.show(block=0)

        ###
        # Plot the BW-average R3,R4,R5,R6 vs. cluster size
        ###
        if not args.noplot: fig = plt.figure(figsize=(9.,6.))
        #fig, ax = plt.subplots(1, 1)
        gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0])  #Top subplot
        ax2=ax1.twinx() 
        
        #"""
        handles=[]
        #datas=[]
        for rSize in reversed([5,6,7,8]):#[3,4,5,6]
            rkey=float('%.1f'%(rSize))
            data=np.array([[float(key), float(allRings[key][rkey])] for key in allRings.keys()]).T
            ax1.plot(data[0],data[1],label='R%d'%rSize,marker='o',ms=2) # ring per atom  #ms=4, ,fillstyle='bottom'
            #ax2.fill_between(data[0],data[1],y2=0,alpha=1.0,label='R%d'%rSize) #does not look good with many cluster sizes.  #alpha=0.8
            #handles.append(h)
            #datas.append(data)
            


        handles, labels = ax1.get_legend_handles_labels()
        plt.legend(handles=reversed(handles),labels=reversed(labels),ncol=4,loc='upper center') #fontsize='x-small'
        #fig.legend
        
        ax1.set_xlabel('Cluster Size')
        #ax1.set_ylabel('Average Ring Count per Atom')
        ax1.set_ylabel('Average Ring Fraction') #  rings per atom (Rc)                   
        ax2.set_yticks([])
        ax1.set_ylim(0,1)

        #Write the  BW-rings vs cluster size data
        print('Writing the Boltzmann-Weighted rings count per atom data for all cluster sizes in data-mean-rings.dat...')
        with open('data-mean-rings.dat','w') as f:
            f.write('#cSize, R3-R10 counts per atom \n')

            keys=sorted([int(x) for x in allRings.keys()])#Re to R20 data available.
            for key in keys: #loop over cSizes
                key=str(key)
                cr=allRings[key]
                f.write('%6s %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n'%(key,cr[3.0],cr[4.0],cr[5.0],cr[6.0],cr[7.0],cr[8.0],cr[9.0],cr[10.0]))
        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)
        


    ###
    # Plot Boltzmann-weighted (and averaged) RDF vs. distance  for each cluster size.
    ###
    if args.type.intersection(['rdf','all']):
        if not args.noplot:         fig = plt.figure(figsize=(9.,6.))

        #Sort the data first
        keys=sorted([int(x) for x in energy.keys()])
        for key in keys:
            key=str(key)
            BW=np.array(boltz_dist(energy[key],T=args.temp)) #Boltzmann weight for each structure in a given size
            #rdfs[key] has the size of (#str for a given size x #bins x 2)
            rdf=np.array(rdfs[key])
            sh=np.shape(rdf[-1])  
            ardf=np.zeros((sh)).T #avg rdf
            noStr=len(rdf)
            for i,r in enumerate(rdf): #loop over structures
                r=r.T
                if i==0:ardf[0]+=r[0]
                if not np.isclose(BW[i],0.0,atol=1e-4):
                    ardf[1]+=r[1]*BW[i] #CHECK THIS (use the one similar to RINGS version). No need!
                    #Plot individual contribution of each structure
                    if 0: plt.plot(r[0],r[1],label='C%s-%s-%.1f'%(key,i,BW[i]),lw='%.3f'%BW[i])
            ardf[1]/=float(len(rdf)) 
            ardf[1]/=np.max(ardf[1]) #normalise
            plt.plot(ardf[0],ardf[1],label='C%s'%(key),lw=2) #WARNING there may be too many labels !!!
            plt.legend(ncol=3,fontsize='x-small') #loc='upper center'
            plt.xlabel('Distance (r) [A]')
            plt.ylabel('Mean g(r)')
            plt.xlim(1,8)

            np.savetxt('data-mean-rdf-C%s.dat'%key,ardf.T)
            
        plt.title(cwd_label)
        if not args.noplot: plt.show(block=0)        

    
    print("Elapsed time: %.2f sec."%( time.time()-initT))
    if not args.noplot: plt.show(block=1)        









    
###########################
#      OLDER VERSION      #
###########################

if not args.useOld: exit()

#########################################
# The only-ground-state-structures loop #
#########################################
#Iterate over minimum-energy structures for each composition
data=Popen4('ca -s -cl | sort -n -k 3')[0]

myDict={};btype=[];pores=[];mods=[];shapes=[]
for st,ln in enumerate(data):
    try:
        words=ln.split()
        fname=words[0].split('/')[-1].split('.')[0]
        totalE=float(words[1])
        cSize=sum([int(x) for x in words[2:-5]]) #account for multiple atom types!!! Additional Seperate columns for the number of atoms from each type
        if not(clist=='all' or cSize in clist): continue

        cType=str(words[-5])
        symm=str(words[-4])
        NoT_found=int(words[-3])
        NoT_searched=int(words[-2].split('-')[-1])#tilda at the end

        myDict[str(cSize)]=[totalE/float(cSize),symm,fname]

    except: continue; #raise Error; continue

    if args.save: system('cp %s.res %s-final.res %s/ 2> /dev/null '%(fname,fname,args.sfold))
    #Get the coord number for atoms in the minimum energy structures for each cluster size.
    try:atoms=ase.io.read(fname+'.res')
    except:atoms=ase.io.read(fname+'-final.res')
    
    #system('ase-file-converter-bk.py -i %s.res -ot xyz  &> /dev/null'%fname)
    #pipeline = import_file(fname+'.xyz') #can import multiple files at once
    #odata = pipeline.compute()
    odata=ase_to_ovito(atoms)
    #print(odata.particles.count);stdout.flush()

    
    vol=atoms.get_volume() #in A^3
    conv=1.6605391 #amu/A^3 to g/cm^3

    if args.clVol:#Get the cluster volume (rather than using the volume of the box) via Voronoi Tesselation (i.e. volume covered by the Voronoi surface) implemented in OVITO visualiser.
        #system('ase-file-converter-bk.py -i %s.res -ot xyz  &> /dev/null'%fname)
        if clflag:  vol=getVol_OVITO(odata,args.cutoff)#'%s.xyz'%fname)
        else: vol=float(Popen4('ovito-voro.py %s.xyz'%fname)[0][1].split(':')[-1])
        #print (fname,len(atoms),vol)

    if vol==0.0:continue
    dens=sum([atoms[at].mass for at in range(len(atoms))])/vol*conv

    modifier = CoordinationAnalysisModifier(cutoff = args.cutoff, number_of_bins = 200,partial=0)
    #pipeline.modifiers.append(modifier)
    odata.apply(modifier)
    #print(odata.tables['coordination-rdf'].xy())

    #i = neighbor_list('i', atoms, cutoff);  coord = np.bincount(i) #this does not count the no of neighbours correctly. use ovito instead!!!
    coord=odata.particles.coordination #for each atom
    #print(coord);stdout.flush()

    #TODO:filter out the Oxygen ones!!
        
    sp=np.count_nonzero(coord==2)
    sp2=np.count_nonzero(coord==3)
    sp3=np.count_nonzero(coord==4)
    oth=np.count_nonzero((coord<=1) | (coord>4 ))
    #if args.delUnPhysical and oth>0: continue # just skip the structures with unphysical coordinations.
    tot=float(sp+sp2+sp3+oth)
    if tot==0: btype.append([cSize,0.,0.,0.,0.,dens])
    else:btype.append([cSize,sp3/tot,sp2/tot,sp/tot,oth/tot,dens])

    if args.pores:
        #Get the porosity info using Zeo++
        if not os.path.exists('%s.cif'%fname):system("ase-file-converter-bk.py -i %s.res -ot cif  -ow &>/dev/null"%fname)
        if not os.path.exists('%s.zeoout'%fname):  system("network-zeo++ -ha -sa 1.82 1.82 1000 -chan 1.82  -vol 1.82 1.82 1000 %s &> %s.zeoout"%(fname+'.cif',fname))
        dens=0;SAVol=0;NSAVol=0;SAVolFrac=0;NSAVolFrac=0
        with open('%s.vol'%fname, 'r') as f:
            ln=f.readlines()[0]
            words=ln.split()
            #print(words)
            for i in range(2,len(words),2):
                word=words[i]
                if word=='Density:':
                    dens=float(words[i+1])
                elif word=='AV_Volume_fraction:' :
                    SAVolFrac=float(words[i+1])
                elif word=='AV_cm^3/g:' :
                    SAVol=float(words[i+1])
                elif word=='NAV_Volume_fraction:':
                    NSAVolFrac=float(words[i+1])
                elif word=='NAV_cm^3/g:':
                    NSAVol=float(words[i+1])
                #print(dens,SAVol,NSAVol,SAVolFrac,NSAVolFrac)
        pores.append([dens,SAVol,NSAVol,SAVolFrac,NSAVolFrac,cSize])

    if args.shape:
        
        if not os.path.exists('%s.data'%fname): system("ase-file-converter-bk.py -i %s.res -ot lammps-data -ow &> /dev/null "%fname)

        #prepare LAMMPS input (TODO: automate the mass input!!)

        str1="""
units           metal
dimension       3
boundary        p p p
#atom_style      atomic
atom_style charge

#read in structure data file
read_data %s.data


#This does not matter as we don't do any simulation step, only sturcture analysis.
pair_style cedip
pair_coeff * * C

mass 1 12.0107

thermo 1

#set atom all image 0 0 0
#compute myCoord all coord/atom cutoff %.2f
compute Rg      all  gyration
compute 1 all gyration/shape Rg
#dump 1 all custom 1 file.trj id type x y z c_myCoord

thermo_style custom time step vol density c_Rg c_1[1] c_1[2] c_1[3] c_1[4] c_1[5] c_1[6]

run 0

#This compute (gyration/shape) calculates a global vector of length 6, which can be accessed by indices 1-6. The first three values are the eigenvalues of the gyration tensor followed by the asphericity, the acylindricity and the relative shape anisotropy.

"""%(fname,args.cutoff)
        lammps_exe='lmp_serial'
        
        if not os.path.exists('./%s.shapeout'%fname):
            open('in.shape','w').write(str1)
            cmd="%s <in.shape>%s.shapeout"%(lammps_exe,fname)
            system(cmd);# print(cmd)
            
        x=[float(xx) for xx in grep('           0        0','%s.shapeout'%(fname))[0].split()]
        Rg=x[4];asp=x[-3];acyl=x[-2];aniso=x[-1]
        shapes.append([cSize,dens,Rg,asp,acyl,aniso])
        
        
    if args.elasticity:
        bmod=0.
        if not os.path.exists('%s.data'%fname): system("ase-file-converter-bk.py -i %s.res -ot lammps-data -ow &> /dev/null "%fname)
        
        cmd="%s <in.elastic>%s.eout"%(lammps_exe,fname)
        
        if not os.path.exists('./%s.eout'%fname):
            open('fname.mod','w').write('read_data %s.data'%fname)
            print(cmd)
            system(cmd)

        try:
            x=grep('Bulk Modulus = ','%s.eout'%fname)[0]
            bmod=float(x.split()[-2])
        except: #try once again
            system(cmd)
            x=grep('Bulk Modulus = ','%s.eout'%fname)[0]
            #if x!='':
            try:bmod=float(x.split()[-2])
            except: print('Error with the elasticity calculation for %s.data'%fname);bmod=0.

        try:smod=float(grep('Shear Modulus 2 = ','%s.eout'%fname)[0].split()[-2])
        except:print('Error with the elasticity calculation for %s.data'%fname);smod=0.
        mods.append([cSize,dens,bmod,smod])

#print (mods)
        

with open('hybrid-data', 'w') as f:
    f.write('#cSize, density, sp3, sp2, sp, others\n')
    for b in btype:
        f.write('%d %8.3f %.3f %.3f %.3f %.3f\n'%(b[0],b[-1],b[1],b[2],b[3],b[4]))

print("Elapsed time: %.2f sec."%( time.time()-initT))

############################################
# PLOTTING PART for Ground-State-Only Loop #
############################################

font = {'family': 'serif', 'size': 18}
plt.rc('font', **font)
#fig = plt.figure(figsize=(11.69, 8.27))
fig = plt.figure(figsize=(9.,6.))
#fig = plt.figure()
gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
ax1 = plt.subplot(gs[0])  #Top subplot
#ax2 = plt.subplot(gs[1]  , sharey=ax1)
ax2=ax1.twiny()


if args.clVol:ax1.set_xlabel('Density [g/cm^3]');ax2.set_xticklabels([])
else:ax1.set_xlabel('Cluster size');ax2.set_xlabel('Density [g/cm^3]')
ax1.set_ylabel('Fraction')

#ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))
#ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,100))
plt.ylim(0,1)


btype=np.array(btype).T
if args.clVol:
    ax1.plot(btype[-1],btype[1],label='sp3',marker='o',lw=0)
    ax1.plot(btype[-1],btype[2],label='sp2',marker='o',lw=0)
    ax1.plot(btype[-1],btype[3],label='sp',marker='o',lw=0)
    ax1.plot(btype[-1],btype[4],label='others',marker='o',lw=0)
else:
    ax1.plot(btype[0],btype[1],label='sp3')#,marker='o',lw=0) marker is needed for multi-atom-type ones
    ax1.plot(btype[0],btype[2],label='sp2')#,marker='o',lw=0)
    ax1.plot(btype[0],btype[3],label='sp')#,marker='o',lw=0)
    ax1.plot(btype[0],btype[4],label='others')#,marker='o',lw=0)

    ax2.plot(btype[-1],btype[1],label='sp3',lw=0)
    ax2.plot(btype[-1],btype[2],label='sp2',lw=0)
    ax2.plot(btype[-1],btype[3],label='sp',lw=0)
    ax2.plot(btype[-1],btype[4],label='others',lw=0)

ax1.legend()
if not args.noplot: plt.show(block=0)
#plt.close(fig)

pores=np.array(pores).T
if args.pores:
    fig = plt.figure(figsize=(9.,6.))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    #ax2 = plt.subplot(gs[1]  , sharey=ax1)
    ax2=ax1.twinx()

    #[dens,SAVol,NSAVol,SAVolFrac,NSAVolFrac,cSize]
    ax1.plot(pores[-1],pores[1]+pores[2],label='Total Void Volume',c='k',marker='o',lw=0)
    ax2.plot(pores[-1],pores[3]+pores[4],label='Void Fraction',c='b',ls='--',marker='o',lw=0)
    if 1:
        ax1.plot(pores[-1],pores[1],label='SA Volume',c='b')#,marker='o',lw=0) marker is needed for multi-atom-type ones
        ax1.plot(pores[-1],pores[2],label='NSA Volume',c='r')#,marker='o',lw=0)
        ax2.plot(pores[-1],pores[3],label='SA Vol Frac',c='b',ls='--')#,lw=0)
        ax2.plot(pores[-1],pores[4],label='NSA Vol Frac',c='r',ls='--')#,lw=0)

    
    #ax1.set_xlabel('Density [g/cc]')
    ax1.set_xlabel('Cluster size')
    ax1.set_ylabel('Volume [cm^3/g]')
    ax2.set_ylabel('Fraction')

    #fig.legend()
    fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.
    if not args.noplot: plt.show(block=0)
    #plt.close(fig)


mods=np.array(mods).T
if args.elasticity:
    fig = plt.figure(figsize=(9.,6.))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    #ax2 = plt.subplot(gs[1]  , sharey=ax1)
    ax2=ax1.twinx()

    #[cSize,dens,bmod,smod]
    ax1.plot(mods[1],mods[2],label='Bulk Modulus',c='k',marker='o',lw=0)
    ax2.plot(mods[1],mods[3],label='Shear Modulus',c='b',ls='--',marker='o',lw=0)
    
    ax1.set_xlabel('Density [g/cc]')
    ax1.set_ylabel('Bulk Modulus [GPa]')
    ax2.set_ylabel('Shear Modulus [GPa]')

    #fig.legend()
    fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.
    if not args.noplot: plt.show(block=0)
    #plt.close(fig)



shapes=np.array(shapes).T
if args.shape:
    fig = plt.figure(figsize=(9.,6.))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    #ax2 = plt.subplot(gs[1]  , sharey=ax1)
    ax2=ax1.twinx()

    #[cSize,dens,bmod,smod]
    ax1.plot(shapes[0],shapes[5],label='Anisotropy',c='k',marker='o',lw=0)
    ax2.plot(shapes[0],shapes[4],label='Acylindricity',c='b',ls='--',marker='o',lw=0)
    
    #ax1.set_xlabel('Density [g/cc]')
    ax1.set_xlabel('Cluster Size')
    ax1.set_ylabel('Anisotropy')
    ax2.set_ylabel('Acylindricity [A^2]')

    #fig.legend()
    fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.
    if not args.noplot: plt.show(block=0)
    #plt.close(fig)


    fig = plt.figure(figsize=(9.,6.))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    #ax2 = plt.subplot(gs[1]  , sharey=ax1)
    ax2=ax1.twinx()

    #[cSize,dens,bmod,smod]
    ax1.plot(shapes[0],shapes[2],label='Gyration',c='k',marker='o',lw=0)
    ax2.plot(shapes[0],shapes[3],label='Asphericity',c='b',ls='--',marker='o',lw=0)
    
    #ax1.set_xlabel('Density [g/cc]')
    ax1.set_xlabel('Cluster Size')
    ax1.set_ylabel('Radius of Gyration [A]')
    ax2.set_ylabel('Asphericity [A^2]')

    #fig.legend()
    fig.legend(ncol=2,loc='upper center') #this shows all legends from all subplots combined.
    if not args.noplot: plt.show(block=0)
    
fig = plt.figure(figsize=(9.,6.))

##########################
# The all-structure Loop #
##########################
#Iterate over all structures for each composition.
myDict2={};minE=1**10;minKey=''
data=Popen4('ca -r -cl ')[0] #-t 100000 option is not needed as by default all structures are shown
#data=Popen4('ca -r -cl | sort -n  -k 5 -k 4') #sort by composition and then relative energy. Not needed
for i,ln in enumerate(data):
    words=ln.split()
    fname=words[0]
    P=float(words[1])
    totalE=float(words[2])
    E_atom=float(words[3])
    
    #cSize=sum([int(x) for x in words[2:-5]]) #account for multiple atom types!!! Additional Seperate columns for the number of atoms from each type
    cSize=0
    for wd in words[3:-3]:
        try: cSize+=int(wd)
        except:
            try:cSize+=int(wd[-4:])
            except:continue
    cType=str(words[-3])
    symm=str(words[-2])
    NoT_found=int(words[-1])
    #NoT_searched=int(words[-2].split('-')[-1])#tilda at the end
    """
    cSize=int(words[4]) #TODO:account for multiple atom types!!!
    cType=str(words[5])
    symm=str(words[6])
    NoT_found=int(words[7])#.split('-')[0])
    #NoT_searched=int(words[6])
    """

    #Apply filters:
    #if args.delUnPhysical:
    #    if oth>0: continue # just skip the structures with unphysical coordinations.
    
    if not(clist=='all' or cSize in clist): continue
    key=str(cSize)
    #val=[totalE/float(cSize),symm,fname]
    #if key!=minKey:val=[E_atom+myDict[key][0],symm,fname]
    if E_atom>=0.: val=[E_atom+myDict[key][0],symm,fname]#add the relative energy to the minimum energy for the given composition.
    else:val=[E_atom,symm,fname]
    if key in myDict2: myDict2[key].append(val)
    else: myDict2[key]= [val]
    
    if val[0]<minE:minE=val[0];minKey=key
     

############################################
# Plotting Part for the All-Structure-Loop #
############################################
keys=myDict2.keys()
keys=sorted(keys)
for key in keys:
    vals=myDict2[str(key)]
    #print(key,vals)
    nosymm=[[],[]]
    symm=[[],[]]
    #nosymm=[];symm=[]
    for val in vals: #Energy,symm,fname
        #print(val)
        if 0: sc=23.0 #kcal/mol to eV conversion #needed for ReaxFF(i.e. real units)
        else:sc=1.
        if val[1] in ['C1','P1']: nosymm[0].append(float(key));nosymm[1].append(val[0]/sc)
        else : symm[0].append(float(key));symm[1].append(val[0]/sc) #symm.append([key,val])

    #plt.plot(nosymm[:][0],nosymm[:][1],marker='x',color='k')
    plt.plot(nosymm[0],nosymm[1],marker='x',color='k')
    plt.plot(symm[0],symm[1],marker='o',color='r')
    #plt.plot(nosymm[:][0],nosymm[:][1],marker='o',color='k')


if args.clVol: plt.xlabel('Density [g/cm^3]')
else:plt.xlabel('Cluster size')
plt.ylabel('Energy [eV/atom]')
#print("Elapsed time: %.2f sec."%( time.time()-initT))

if not args.noplot: plt.show(block=1)




#DELETED PARTS:

"""
        #TODO:Take this out of loop for efficiency
        pipeline = import_file(seed+'.xyz') #needs to be coverted from res format
        pipeline.add_to_scene()
        odata2=pipeline.compute()
        #odata.vis.DataVis(1)
        vp = Viewport()
        vp.type = Viewport.Type.Perspective
        vp.camera_pos = (-100, -150, 150)
        vp.camera_dir = (2, 3, -3)
        #vp.fov = math.radians(60.0)
        vp.render_image(size=(800,600), filename="%s.png"%seed, background=(0,0,0), frame=-1)
        #vp.render_anim(size=(800,600), filename="animation.avi", fps=20)
        #Save the visuals of the system
        odata.apply(ColorCodingModifier(
            property = 'Coordination',
            gradient = ColorCodingModifier.Jet()
))

        """
 
