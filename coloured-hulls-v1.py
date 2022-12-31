#!/usr/bin/env python3
#Adapted from Angela Harpers script.

import matplotlib,time
#matplotlib.rcParams['figure.figsize'] = (8, 6)
matplotlib.rcParams['savefig.dpi'] = 300
# %matplotlib inline 
#from matador.query import DBQuery
#from matador.hull import QueryConvexHull
#from matador.similarity.similarity import get_uniq_cursor
#from matador.scrapers.castep_scrapers import res2dict
#from matador.utils.cursor_utils import filter_cursor, get_array_from_cursor
#from matador.export import doc2res
#from ilustrado.ilustrado import ArtificialSelector
#import seaborn as sns

import numpy as np
import matplotlib.pyplot as plt

from os import chdir, makedirs
from glob import glob

from sys import argv,exit
from matador.utils.chem_utils import get_root_source

from matador.utils.cursor_utils import get_array_from_cursor
from matador.query import DBQuery
from matador.hull import QueryConvexHull
from matador.hull import EnsembleHull
from matador.scrapers import castep2dict, res2dict, bands2dict
from matador.utils.cursor_utils import filter_unique_structures, display_results
from matador.crystal import Crystal
from matador.utils.cell_utils import standardize_doc_cell, get_spacegroup_spg, get_space_group_label_latex
from matador.utils.chem_utils import get_formula_from_stoich
from matador.hull.hull_temperature import TemperatureDependentHull
from matador.orm.spectral import VibrationalDOS, ElectronicDispersion
import matador.plotting
from matador import __version__
print(f"This notebook run was performed with matador {__version__}.")

import argparse

#from matador.utils.cursor_utils import get_guess_doc_provenance
def get_guess_doc_provenance(sources, icsd=None):
    """ Returns a guess at the provenance of a structure
    from its source list.

    Return possiblities are 'ICSD', 'SWAP', 'OQMD' or
    'AIRSS'.
    """
    prov = 'AIRSS'
    if sources is str:
        sources = [sources]
    for fname in sources:
        #if not isinstance(fname,str): fname=fname[0]
        if (fname.endswith('.castep') or fname.endswith('.res') or
                fname.endswith('.history') or 'OQMD' in fname):
            if 'collcode' in fname.lower() or 'colcode' in fname.lower() or 'collo' in fname.lower():
                if fname.split('/')[-1].count('-') == 2 + fname.lower().count('oqmd'):
                    prov = 'SWAPS'
                else:
                    prov = 'ICSD'
            elif 'swap' in fname.lower():
                prov = 'SWAPS'
            elif '-ga-' in fname.lower():
                prov = 'GA'
            elif icsd is not None:
                prov = 'ICSD'
            elif 'oqmd' in fname.lower():
                prov = 'OQMD'
            elif '-icsd' in fname.lower():
                prov = 'ICSD'
            elif 'mp-' in fname.lower() or ('mvc-' in fname.lower() and '_crystal' in fname.lower()):
                prov = 'MP'
            #elif 'config-enum' in fname.lower():
            elif 'config_' in fname.lower():
                prov = 'ConfigEnum'
            else:
                prov='AIRSS'
    return prov

def plot_point(doc,marker='o',label=None,side='black',z=10000000000003,color='white'):
    #formation_enthalpy = get_array_from_cursor(doc, 'formation_enthalpy_per_atom')
    if label is not None:
        try:ax.scatter(doc['concentration'][0], doc['formation_enthalpy'], marker=marker,  s=50, c=color, zorder=z, label=label) #edgecolors=side, lw=1.5, #orig size=50
        except:print('Error in %s'%doc['source'][0])
    else:
        try:ax.scatter(doc['concentration'][0], doc['formation_enthalpy'], marker=marker, s=50, c=color, zorder=z) #edgecolors=side, lw=1.5, #handletextpad =??
        except:print('Error in %s'%doc['source'][0])



#from matador.utils.cursor_utils import get_array_from_cursor
## here we grab the first element of the concentration array, which only contains 1 element for a binary system
#concentrations = get_array_from_cursor(hull.hull_cursor, ['concentration', 0])
#formation_enthalpy = get_array_from_cursor(hull.hull_cursor, 'formation_enthalpy_per_atom')
#hull_distances = get_array_from_cursor(hull.hull_cursor, 'hull_distance')
#for conc, eform, hdist in zip(concentrations, formation_enthalpy, hull_distances):
#    print(f"{conc:12.4f} {eform:12.4f} {1000*hdist:12.4f}")

initT=time.time()

parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')

parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')

parser.add_argument("-c","--comp", type=str,required=1,
                    help="Chemical composition to be analysed. No default")
args = parser.parse_args()



#query=DBQuery(**{'composition': 'LiVNbO', 'intersection': True, 'subcmd': 'query', 'biggest': True, 'details': 0, 'cutoff': [300,301], 'kpoints': 0.07, 'kpoint_tolerance': 0.1, 'source': 0,'summary':1,'no_plot':False , 'hull_cutoff': 0.00, 'db': 'bk393-LiPS'})  # 
#cursor=query.cursor

cursor, failures = res2dict(args.inpf, as_model=True)


# filter for uniqueness
filtering = 0
if filtering:
    polished_cursor = filter_unique_structures(cursor, sim_tol=0.1, enforce_same_stoich=True, quiet=True)
else:
    polished_cursor = cursor

# do some pruning: reevaluate symmetries and reduce cells
polished_cursor = [Crystal(standardize_doc_cell(doc)) for doc in polished_cursor]


polished_hull = QueryConvexHull(
    cursor=polished_cursor, 
    #species='Li3NbO4:LiVO2',
    species=args.comp,
    no_plot=1,
    #hull_cutoff=0.05,
    hull_cutoff=1.5,
)




import matplotlib.pyplot as plt
#import seaborn as sns


font = {'family': 'Arial', 'size': 24, 'weight':'normal'}
plt.rc('font', **font)
fs=24
plt.rc('axes', titlesize=fs)     # fontsize of the axes title
plt.rc('axes', labelsize=fs)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=fs)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fs)    # fontsize of the tick labels

fig, ax = plt.subplots(figsize=(8,6))
from matador.plotting import plot_2d_hull#,plot_3d_hull

#polished_hull.set_plot_param()
ax = plot_2d_hull(polished_hull,ax=ax, show=False,plot_points=False)
"""
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(18)
    #item.set_fontsize(32)
"""
counts={};colors={};fnames=[];minEf=1000.0
#for doc in cursor:
for doc in polished_hull.cursor:
    if 'BKP' in doc['source'][0]:continue
    #delete duplicates based on filename: even faster
    fname=doc['source'][0].split('/')[-1]
    if fname in fnames:continue
    fnames.append(fname)
    
    prov = get_guess_doc_provenance(doc['source'])
    if prov == 'ICSD':marker='o';color='blue';z=1
    elif prov == 'OQMD':marker='o';color='black';z=2
    elif prov == 'SWAPS':marker='o';color='red';z=3
    elif prov == 'MP':marker='o';color='magenta';z=4
    elif prov == 'ConfigEnum':marker='o';color='green';z=7 ; #continue
    elif prov == 'GA':marker='o';color='cyan';z=6
    elif prov == 'AIRSS':marker='o';color='orange';z=7 ; continue
    #black; color='xkcd:gold'
    if doc['space_group'] == "P1":marker='x';#continue
    st=(doc['stoichiometry'])
    """ #Normally we have this part.
    if len(st)==3:#Order is Li, P, S
      if st[0][1] == 7 and st[1][1]==1 and st[2][1]==2: None ;#continue #Li29P9S 
      elif st[0][1] == 8 and st[1][1]==2 and st[2][1]==1: None; #continue #Li32P10S
      elif st[0][1] == 11 and st[1][1]==3 and st[2][1]==1: None;#continue#Li37P11S2
      elif st[0][1] == 5 and st[1][1]==1 and st[2][1]==1: None;#continue#Li37P11S2
      else: continue #only take the above ternaries
    """

    #TODO: gather all data first store them and then plot as a whole (to spped up plotting).
    if prov not in counts: counts[prov]=1
    else: counts[prov]+=1
    
    if counts[prov] == 1:
        #print(dir(doc))
        plot_point(doc,label=prov,color=color,marker=marker,z=z)
    else:
        plot_point(doc,color=color,marker=marker,z=z)

    if prov not in colors: colors[prov]=color
    #if doc['formation_enthalpy']<minEf: minEf=doc['formation_enthalpy']
    
print("Number of structures from each source: ")
provs=sorted(counts.keys())
for x in provs: print("%s: %d"%(x,counts[x]))
print ("Total: %d"%(sum(counts.values())))

if 0:#Filter out same structures based on PDF overlap.
    print('Filtering unique structures... {}'.format(len(cursor)))
    uniq_list, _, _, _ = list(get_uniq_cursor(cursor[1:-1], debug=False))
    cursor = [cursor[1:-1][ind] for ind in uniq_list]

#mn=-0.10;mx=0.25
mn=-0.05;mx=0.300
ax.set_ylim(mn, mx)
ax.set_yticks([val for val in np.arange(mn,mx,0.05)]) #minEf*100-15
ax.set_yticklabels(['{:.2f}'.format(val) for val in ax.get_yticks()])

ax.set_xticks([0., 0.25, 0.33, 0.5, 0.66 , 0.75, 1.0]) #minEf*100-15
ax.set_xticklabels(['{:.2f}'.format(val) for val in ax.get_xticks()])

fs=28
plt.rcParams['font.size']=fs
#"""
ax.xaxis.label.set_size(fs)
ax.yaxis.label.set_size(fs)
#ax.xtick.label.set_size(fs)
#ax.ytick.label.set_size(fs)
ax.xaxis.set_tick_params(labelsize=fs)
ax.yaxis.set_tick_params(labelsize=fs)
#plt.xticks
#plt.tick_params(labelsize=fs)
ax.title.set_fontsize(fs)
ax.title.set_fontname('Arial')

for tick in ax.get_xticklabels():
    tick.set_fontname("Arial");tick.set_fontweight("normal")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")
#"""     

#"""
for tick in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    tick.set_fontsize(fs);tick.set_fontname("Arial");tick.set_fontweight("normal")
    #item.set_fontsize(32)
#"""



#Prepare the legend
handles=[]
for x in provs:
    #h, = plt.scatter([0],[0], marker='o',c=colors[x])#,edgecolors='black', lw=1.5, s=50)
    h, = plt.plot(0, marker='o',c=colors[x],linestyle='')#,edgecolors='black', lw=1.5)#, s=50)
    handles.append(h)
    
l = ax.legend(handles,provs,loc='lower center', frameon=True, fancybox=True, shadow=False, framealpha=1, facecolor='white', edgecolor='black', ncol=6, mode="expand", fontsize=24,handletextpad=0.2) #loc='lower center'

#l = ax.legend(handles,provs,loc='right', frameon=True, fancybox=True, shadow=False, framealpha=1, facecolor='white', edgecolor='black', ncol=1, fontsize=24,handletextpad=0.2)

plt.rc('font', **font)
fs=24
plt.rc('axes', titlesize=fs)     # fontsize of the axes title
plt.rc('axes', labelsize=fs)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=fs)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fs)    # fontsize of the tick labels

plt.savefig('myhull.png', bbox_inches='tight',dpi=600)
print("Elapsed time: %.2f sec."%( time.time()-initT))

plt.show()


if 0:
  dir_name = 'SAVED'
  for doc in cursor_1:
    try:doc2res(doc, dir_name + '/' + get_root_source(doc))
    except:continue


exit()



"""
plt.rc('font', size=fs)          # controls default text sizes
plt.rc('axes', titlesize=fs)     # fontsize of the axes title
plt.rc('axes', labelsize=fs)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=fs)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fs)    # fontsize of the tick labels
#plt.rc('legend', fontsize=fs)    # legend fontsize
plt.rc('figure', titlesize=fs)  # fontsize of the figure title
"""



hull_1 = QueryConvexHull(cursor=cursor_1,no_plot=True ,intersection= True, biggest= True, hull_cutoff= 0.25 , summary=True ) #, 'cutoff': 300, 'kpoints': 0.07, 'kpoint_tolerance': 0.01,'summary':True})
#hull_1 = QueryConvexHull(**{'cursor'=cursor_1,'no_plot'=True,'intersection': True, 'biggest': True, 'hull_cutoff': 0.25, 'details': True, 'cutoff': 300, 'kpoints': 0.07, 'kpoint_tolerance': 0.01,'summary':True})

#hull_2 = QueryConvexHull(cursor_2)

#cursor = [res2dict(res)[0] for res in glob('seed/*.res')]
#hull = QueryConvexHull(cursor=cursor, no_plot=True, kpoint_tolerance=0.03,
#                       summary=True, hull_cutoff=1.5e-1)

print('Filtering down to only ternary phases... {}'.format(len(hull_1.hull_cursor)))
hull_1.hull_cursor = [doc for doc in hull_1.hull_cursor if len(doc['stoichiometry']) == 3]

print('Filtering unique structures... {}'.format(len(hull_1.hull_cursor)))
uniq_list, _, _, _ = list(get_uniq_cursor(hull_1.hull_cursor[1:-1], debug=False))
cursor = [hull_1.hull_cursor[1:-1][ind] for ind in uniq_list]

print('Final cursor length... {}'.format(len(cursor)))
print('over {} stoichiometries...'.format(len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in cursor]))))
print([doc['stoichiometry'] for doc in cursor])

dir_name = 'TOPLA'
for doc in cursor_1:
    try:doc2res(doc, dir_name + '/' + get_root_source(doc))
    except:continue
    
    
exit()

#dict_keys(['_id', 'source', 'user', 'atom_types', 'positions_frac', 'num_atoms', 'stoichiometry', 'num_fu', 'external_pressure', 'space_group', 'num_kpoints', 'kpoints_mp_grid', 'species_pot', 'geom_force_tol', 'elec_energy_tol', 'finite_basis_corr', 'cut_off_energy', 'xc_functional', 'task', 'spin_polarized', 'optimised', 'lattice_cart', 'lattice_abc', 'cell_volume', 'enthalpy', 'enthalpy_per_atom', 'max_force_on_atom', 'stress', 'pressure', 'mulliken_charges', 'mulliken_spins', 'kpoints_mp_spacing', 'castep_version', 'date', 'estimated_mem_MB', 'total_time_hrs', 'peak_mem_MB', 'tags', 'text_id', 'quality', 'elems', 'num_chempots', 'concentration', 'formation_enthalpy', 'hull_distance'])

#print(cursor[0].keys())
#print(polished_hull.cursor[0].keys())
#cursor=[doc for doc in cursor if doc['hull_distance'] <0.20]

"""
query = DBQuery(subcmd='query', composition=['CuLi'],db=['afh41-CuP'],intersection=True,kpoint_mp_spacing=0.05,cutoff=[300,301],summary=True)
cursor=query.cursor
polished_hull = QueryConvexHull(cursor=cursor, subcmd='hull', no_plot=1, composition=['CuLi'], elements=['Cu', 'Li'],summary=True, spin=0,details=False,labels=False,source=False,loose=False,intersection=True,cutoff=[300,301],hull_cutoff=0.15)
cursor=polished_hull.cursor #this one has the concentration fields in the cursor dict
"""

"""
def plot_point(doc,label=None,marker='o',side='black',z=10000000000003,color='white'):
    if label is not None:
        ax.scatter(doc['concentration'][0], doc['formation_enthalpy'], marker=marker, edgecolors=side, lw=1.5, s=50, c=color, zorder=z, label=label)
    else:
        ax.scatter(doc['concentration'][0], doc['formation_enthalpy'], marker=marker, edgecolors=side, lw=1.5, s=50, c=color, zorder=z)
"""

"""
    if  prov == 'ICSD': 
        #count += 1
        if counts[prov] == 1:
            plot_point(doc,label='ICSD',color='blue',marker='o',z=9900000000)
        else:
            plot_point(doc,color='blue',marker='o',z=9900000000)
    elif prov == 'OQMD':
            #_count += 1
            if counts[prov] == 1:
                plot_point(doc,color='xkcd:gold',label='OQMD')
            else:
                plot_point(doc,color='xkcd:gold')
    elif prov == 'SWAPS':
            #__count += 1
            if counts[prov] == 1:
                plot_point(doc,color='red',label='SWAPS')
            else:
                plot_point(doc,color='red')
    elif prov == 'AIRSS':
        #___count += 1
        if counts[prov] == 1:
            plot_point(doc,color='cyan',z=100000000,label='AIRSS')
        else:
            plot_point(doc,color='cyan',z=100000000)
            
    elif prov == 'GA':
        #____count += 1
        if counts[prov] == 1:
            plot_point(doc,color='violet',marker='o',z=10000000,label='GA')
        else:
            plot_point(doc,color='violet',marker='o',z=10000000) 
       """     
