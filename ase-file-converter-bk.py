#!/usr/bin/env python3
#Written by Bora Karasulu, requires ASE.
from __future__ import print_function

import ase.io,time
import argparse,os.path,re
from ase import Atoms
from ase.build import sort,surface,add_adsorbate,minimize_tilt
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

import ase.io.res
from ase.visualize import view

import random
from ase.build import molecule
from ase.geometry import wrap_positions


from ase.geometry import find_mic  #min dist through PBC using the minimum-image covnenti
def minDist_ASE(arr1,arr2,latt): #needs np.array inputs
    if len(arr1)==4: arr1=arr1[1:4]
    if len(arr2)==4: arr2=arr2[1:4]
    D,x=find_mic(np.array([arr1-arr2]),latt,pbc=True)
    return x[0]

def get_random_ab(cell, ntimes, tol):
    ab = []
    ab.append(np.random.rand(2))
    for i in range(1,ntimes+1):
        not_accepted = True
        cnt=0
        while not_accepted:
            if cnt>100:print('not possible to add the adatoms with the given parameters....');break
            cnt+=1
            # generate_random_coordinate
            atmp, btmp = np.random.rand(2)
            #xtmp, ytmp, dummy = np.dot(cell.transpose(), np.array([atmp, btmp, 0.0]))
            tmp = np.dot(cell.transpose(), np.array([atmp, btmp, 0.0]))
            xtmp, ytmp, dummy = tmp
            #print('xtmp ytmp dummy', xtmp, ytmp, dummy)
            # compute distance with all previous coordinates
            r_prev = []
            for ir in ab:
                r1x, r1y, r1z = np.dot(cell.transpose(), np.array([ir[0], ir[1], 0.0]))
                r1 = np.array([r1x, r1y,0.])
                r2 = np.array([xtmp, ytmp,0.])
                #print('r2', r2)
                #r12 = np.linalg.norm(r2-r1)
                r12=minDist_ASE(r2,r1,cell)
                r_prev.append(r12)

            within_tol = True
            for ir in r_prev:
                if ir < tol:
                    within_tol = False
                    break

            if within_tol:
                not_accepted = False
                ab.append(np.array([atmp, btmp]))

    return np.array(ab)

def add_molecule_to_surface(surface, molecule, ntimes=1, nlayers=1,random_ab=False, ab=None, tol=2.0,ori=None):
    # add molecule to the surface
    new_surface = surface.copy()

    for l in range(1,nlayers+1):
        cell = new_surface.cell
        zmax = np.max(new_surface.get_positions()[:,2])
        mol_z = zmax + args.adsorb_height
        mol_c = mol_z/cell[2,2]

        if random_ab:
            ab = get_random_ab(cell, ntimes=ntimes, tol=tol)

        for i in range(ntimes):
            mol_a, mol_b = ab
            mol_pos = np.dot(cell.transpose(), np.array([mol_a, mol_b, mol_c]))
            # translate molecule to new position
            molecule.translate(mol_pos - molecule.get_positions()[0])

            axes=['x','y','z']
            if not ori: #no specific orientation is given so randomly placed
              #Rotate randomly
              axis=axes[np.random.randint(3, size=1)[0]]
              angle=360*np.random.random()-180 #sample btw [-180,180)
              print('Random orientation of the adsorbate on the surface: ', axis,angle)
              molecule.rotate(angle,axis,center='COP')

            else:#take rotation angles from the user
              molecule.rotate(ori[0],'x',center='COP')
              molecule.rotate(ori[1],'y',center='COP')
              molecule.rotate(ori[2],'z',center='COP')

            new_molecule = ase.Atoms(molecule, cell=surface.cell)
            new_molecule.wrap()
            new_surface += new_molecule

    return new_surface




def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def find_interstices(atoms,tol=1e-6,sg=None,verb=None,prim=None):
    from fractions import Fraction
    from ase.spacegroup import Spacegroup,get_spacegroup

    #get spacegroup
    if sg:        SG=Spacegroup(sg)
    else:        SG=get_spacegroup(atoms, symprec=tol)
    if verb: print ("Original space group: ",SG.no,SG.symbol)

    scpos=atoms.get_scaled_positions()
    if verb: print('Original scaled/fractional positions: \n',scpos)   

    if verb:print ('\nReading wyckoff.dat...') #This can be obtained from the VESTA source code package.
    data={}
    with open(environ['HOME']+"/wyckoff.dat",'r') as file:
        flag=0
        i=0
        tmp={}
        for ln in file:
            if '#' in ln[0:1]: continue
            i+=1
            x=ln.split()
            if ln.strip()=='':
                if int(sg.split('-')[0]) != SG.no: tmp={};continue #only store the relevant space group.
                data[sg]=tmp;tmp={};flag=0; continue; #print (sg,data[sg]);continue
            elif len(x)==3 and ln[0] !='(' :wyck='%s%s'%(x[0],x[1]);tmp[wyck]=[]
            elif len(x)>5 and ln[0] !='(' : sg='%s-%s'%(x[0],x[1]) ; 
            elif ln[0]=='(': 
                for crd in x: #TODO: Double check if the x,y,z and addition/substraction of fractions are readand handled correctly (when readin from wyckoff.dat).
                    crd_bkp=deepcopy(crd)
                    crd=crd.strip('(');crd=crd.strip(')')
                    crd=crd.replace('-x','');crd=crd.replace('-y','');crd=crd.replace('-z','')
                    crd=crd.replace('+x','');crd=crd.replace('+y','');crd=crd.replace('+z','')
                    crd=crd.replace('x','');crd=crd.replace('y','');crd=crd.replace('z','')
                    #tmp[wyck].append([float(a) if ['x','y','z'] in a  else 1 for a in crd.split(',')  ])
                    try:tmp[wyck].append([float(Fraction(a)) if a != ''  else 0. for a in crd.split(',')  ])
                    except: print ('Error! :  ', sg,wyck,crd_bkp,crd)


    data2={}
    for key in data.keys(): 
        data2[key]={}
        for w in data[key].keys():
            newcrd=data[key][w][0]
            equiSites=SG.equivalent_sites(newcrd,onduplicates='replace')[0] #Equi site of the first site should cover the other sites as well.
            #print (w,equiSites)
            data2[key][w]=[]
            for crd in equiSites: #SG.equivalent_sites([data[key][w][0]])[0]:
                #if crd not in scpos: #This one is supposed to work but somehow misses some obvious matches (due to the difference being larger than the default value).. Therefore switched to allclose comparison for each array element.
                flag=0
                for sc in scpos:
                    if np.allclose(sc, crd):flag=1;break 
                if not flag:
                    #print('%s,%s: %s not in atoms... %s '%(key,w,crd, equiSites))
                    data2[key][w].append(crd)

    a,b,c= atoms.get_cell()
    for i,key in enumerate(sorted_nicely(data2.keys())):
        if i==1:break #Take only the first setting
        if verb: print(key)
        sites=[];    Va_sites=[]
        for w in sorted_nicely(data2[key].keys()):
            if len(data2[key][w])!=0: 
                if verb: print (w,[" ".join(['%.3f' %y for y in x])  for x in data2[key][w]])
                #for val in data2[key].values():
                sites.extend(data2[key][w])   

        ulist=np.unique(sites,axis=0)
           
        flag=0 #Sometimes np.unique fails to detect equal arrays/sites, so this is a fail-safe approach to delete repeating sites. 
        for crd in ulist:
            flag=0
            for u in Va_sites: 
              if np.allclose(crd,u):flag=1;break 
            if flag:continue
            r=crd[0]*a+crd[1]*b+crd[2]*c #fractional to Cartesian conversion, works!!
            atoms.append(Atom('X',r)) 
            Va_sites.append(crd)

    if verb: print ('\n%d unique vacancy sites detected.'%(len(Va_sites)))#, Va_sites))
    if verb:  
      for x in Va_sites:print (x)
    if verb: print (atoms)
    SG=get_spacegroup(atoms, symprec=tol)
    print ("\nNew space group after adding all interstices: ",SG.no,SG.symbol)

    if verb: print('New scaled/fractional positions: \n',atoms.get_scaled_positions())

    return atoms


parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                    help='an integer for the accumulator')


parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')
parser.add_argument('-ot','--otype', type=str,required=True, help='Output file type')
parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: input name is used as root.')
parser.add_argument('-od','-odir','--odir', type=str,required=False, default='./',help='Output file directory. Def: ./ ')
parser.add_argument('-it','--itype', type=str,required=False, help='Input file type. Def: determined automatically from the extension.')
parser.add_argument('-prim','--ifPrim', default=False,action='store_true', help='print out the primitive cell instead of supercell (requires SPGlib)')
parser.add_argument('-kpath','--kpath', default=False,action='store_true', help='prints out the high-symmetry k-point path for the given symm (requires seekpath package).')
parser.add_argument('-t', '--tol',type=float, default=1e-4,help="The symmetry tolerance.")

parser.add_argument('-all','--ifAll', default=False,action='store_true', help='Consider all the geometries in the input file (e.g. trajectories). Def: only last geometry is considered.')

parser.add_argument('-sep','--separate', default=False,action='store_true', help='Write the multiple steps into separate files. Def: No')

parser.add_argument('-ind','--index', type=int,required=False, help='Index of the geometry to be used  in a multistep trajectory. Def: Last one taken')

parser.add_argument('-skip','--skip', type=int,default=1,required=False, help='To skip steps while reading a trajectory. Def: All steps')

parser.add_argument('-uw','--unwrap',default=False,action='store_true', help='Generate unwrapped coordinates, as needed for MSD and VFAC, in case input coordinates are wrapped into cell. Def: No unwrapping done.')

parser.add_argument('-wr','--wrap',default=False,action='store_true', help='Enforce wrapping of coordinates into the cell. Def: No wrapping done.')


parser.add_argument('-r','--rep', nargs=3, type=int,help="Repeat the resulting structure in a, b, c directions, e.g. -r 2 2 1")

parser.add_argument('-sw','--swaps', nargs='*', type=str,default=None,help='To swap of a given atom type with another. Sample use: -sw Li:K Se,O:S V,N,As:P (and so on)')

parser.add_argument('-int','--interstices', default=False,action='store_true', help='Add the misssing interstices to the structure (Uses Wykoff information, from VESTA - wyckoff.dat). Def: No')

parser.add_argument('-cr','--coord', default=False,action='store_true', help='Coordination analysis')

parser.add_argument('-ns','--noSubs', default=False,action='store_true', help='Do not include the substrate (atoms in the inital frame) in the coord analysis.')


parser.add_argument('-f','--fast', default=False,action='store_true', help='Faster and low-memory conversion using the iterative iread, rather than keeping all timesteps in memory.')

parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')

parser.add_argument('-ovito','--useOvito', default=False,action='store_true', help='Use OVITO io module to read the file (useful e.g.when ASE fails for LAMMPS trajectory/dump files. ef:No0')

parser.add_argument('-airss','--airss', default=False,action='store_true', help='Format the geometry to be used in AIRSS input')

parser.add_argument('-molecule','--molecule', default=False,action='store_true', help='Create the LAMMPS molecule input to be used in deposition, pour etc. simulations.')

parser.add_argument('-optim_xtb','--optim_xtb', default=False,action='store_true', help='Create the LAMMPS molecule input to be used in deposition, pour etc. simulations.')

parser.add_argument('-group','--group', type=str, choices=['no'], help="Group the resulting output files in different directories based on e.g. atom numbers, etc. Options: 'no': total number of atoms, 'xx': XX,")

parser.add_argument('-zlim', '--zlim',type=float, nargs=2, help="Only take atoms within a range of z position, e.g. -zlim 10.0 25.0")

parser.add_argument('-sc','--setcell', nargs=3, type=float,help="Modify cell parameters without distorting the structure")

parser.add_argument('-mint','--min_tilt', default=False,action='store_true', help='Minimise the tilt angle of the given cell.Def: No')

parser.add_argument('-surf','--surface', nargs=3, default=None,type=int,help="Create a surface slab with the given Miller indices from the input structure  using ase.build.surface function, e.g. -surf 1 1 1")
parser.add_argument('-vac','--vacuum', default=15.,type=float,help="Vaccum to apply in z/-z directions for the surface model, only relevant when using the -surf keyword")
parser.add_argument('-lay','--layers',  default=1,type=int,help="No of layers to include in the surface model, only relevant when using the -surf keyword")

parser.add_argument('-add','--adsorb', default=None,type=str,help="Add adsorbate, which can be file name (use ASE supported file types) or a atom/molecule name (again supported by ASE), e.g. -add CO or -add CO.xyz")
parser.add_argument('-add_height','--adsorb_height', default=2.0,type=float,help="Distance to surface for the adsorbate")
parser.add_argument('-add_target','--adsorb_target', default=None,nargs="*",type=str,help="Target surface atom to add adsorbate on. Def: added randomly")
parser.add_argument('-add_count','--adsorb_count', default=1,type=int,help="No of adsorbate to add")
parser.add_argument('-add_nlayers','--adsorb_nlayers', type=int, default=1, help='Number of layers of molecules to add')
parser.add_argument('-add_slayers','--adsorb_slayers', type=int, nargs="*",default=[-1,-2], help='Number of layers to be treated as the surface layers to determine the target atoms. Def: first two are treated as surface')
parser.add_argument('-add_tol', '--adsorb_tol', type=float, default=2.0, help='minimum distance between each of the added atoms/molecules')
parser.add_argument('-add_ab', '--adsorb_ab', nargs=2, type=float, default=None, help='position where the molecule should be positioned in the ab plane. Def: placed randomly')
parser.add_argument('-add_ori', '--adsorb_ori', nargs=3, type=float, default=None, help='Orientation of the adsorbate molecule wrt to the surface. Give rotation angles about the x, y and z axes respectively, e.g. -add_ori 0. 0. 0. to keep the original orientation. Def: randomly oriented')



parser.add_argument('-v','-view','--view', default=False,action='store_true', help='View the final structure')


args = parser.parse_args()

initT=time.time()

#if args.group and args.group=='': args.group='no'

if args.ifPrim or args.kpath: from spglib import *; print ("Symmetry is determined using tol = %f"%args.tol)

if args.otype=='pdb':args.otype='proteindatabank'
elif args.otype=='cell':args.otype='castep-cell'
ext=args.otype.split("-")[-1]
if ext=='proteindatabank':ext='pdb'
if args.ifAll: index="::%d"%args.skip
else:
    if args.index: index=args.index
    else:index=None
    
cnt=0
for inpf in args.inpf:
    print (inpf)
    #if inpf in ['CONTCAR','POSCAR','vasprun.xml','OUTCAR','XDATCAR']: args.itype='vasp'

#    try:
    if 1:
        if args.outf: outf=args.odir+args.outf+"."+ext
        else: outf=args.odir+inpf.split(".")[0]+"."+ext

        if os.path.isfile(outf):
            if not args.overwrite: print ("%s already exists, skipping... "%outf);  continue
            else:  print ("%s already exists, overwriting... "%outf)

        #print (index)
        if args.fast:
            if args.itype: ats=ase.io.iread(inpf,format=args.itype,index=index)
            else: ats=ase.io.iread(inpf,index=index) #ASE automatically determines from the extension.
            system('rm %s'%outf)
            cn=0
            for atoms in ats:
                cn+=1
                if args.swaps:
                    for sw in args.swaps:
                        frm=sw.split(":")[0].split(",")
                        to=sw.split(":")[1]

                        for at in atoms:
                            if at.symbol in frm:  atoms[at.index].symbol=to
                        
                if args.otype=="vasp":ase.io.write(outf,atoms,format=args.otype,vasp5=True,append=1)
                elif args.otype=="lammps-data":ase.io.write(outf,atoms,format=args.otype,atom_style='charge',append=1)
                else: ase.io.write(outf,atoms,format=args.otype,append=1)
                
            print('%d steps were read from %s, skipping %d steps.'%(cn,inpf,args.skip));stdout.flush()
            print("Elapsed time: %.2f sec."%( time.time()-initT))


            cnt+=1
            continue
        
        if args.useOvito:
            from ovito.io import import_file
            from ovito.io.ase import ovito_to_ase

            # Create an OVITO data pipeline from an external file:
            pipeline = import_file(inpf)

            # Evaluate pipeline to obtain a DataCollection:
            if args.index:
                #print(args.index)
                data = pipeline.compute(args.index)
                # Convert it to an ASE Atoms object:
                atoms = ovito_to_ase(data)
            elif args.ifAll:
                atoms=[]
                for frame in range(0,pipeline.source.num_frames,args.skip):
                    data = pipeline.compute(frame)
                    atoms.append(ovito_to_ase(data))
            else:
                data = pipeline.compute(pipeline.source.num_frames)
                # Convert it to an ASE Atoms object:
                atoms = ovito_to_ase(data)

                    
        elif args.itype: 
            if args.itype=='vasp':
                from ase.calculators.vasp import Vasp,Vasp2 #requires ase 3.17.0
                from ase.spacegroup import Spacegroup,get_spacegroup
                
                cwd=os.getcwd()
                dr='/'.join(inpf.split('/')[0:-1])
                if dr=='': dr='./'
                os.chdir(dr)
                try:
                    calc = Vasp(restart=True)
                    atoms = calc.get_atoms()
                    print ("VASP run was read succesfully from OUTCAR.")
                except:
                    print("Error in reading VASP run from OUTCAR")
                    os.chdir(cwd)

                try:atoms.info['energy']=atoms.get_total_energy()
                except:None
                atoms.info['name']=inpf
                P=atoms.info.get('pressure')
                if P is not None: atoms.info['pressure']=P
                else: atoms.info['pressure']=0.
                SG=atoms.info.get('spacegroup')
                if SG is None:   SG=str(get_spacegroup(atoms, symprec=args.tol).symbol.replace(' ',''))
                atoms.info['spacegroup']=SG
                atoms.info['times_found']=0

                os.chdir(cwd)

            else:atoms=ase.io.read(inpf,format=args.itype,index=index)
        else: atoms=ase.io.read(inpf,index=index) #ASE automatically determines from the extension.

        if not args.ifAll and args.setcell: 
            atoms.set_cell(args.setcell,scale_atoms=0) 
            atoms.center()

        if args.ifAll: print('%d steps were read from %s, skipping %d steps.'%(len(atoms),inpf,args.skip));stdout.flush()
        
        if args.ifAll and args.unwrap:
            init_pos=atoms[0].get_scaled_positions()
            #subs=[ids.index for ids in atoms[0]]
            ppos=dc(init_pos)
            noAtoms=len(atoms[0])
            shifts=np.zeros((noAtoms,3))
            print ("Unwrapping the coordinates as requested.")
            for st in range(len(atoms)): #over time/geom step.
                #if args.noSubs and st  in subs: continue #skip the substrate atoms
                ats=atoms[st]
                if  st!=0:
                    cpos=ats.get_scaled_positions(wrap=False)
                    for at in range(noAtoms):
                        for j in range(3):
                            diff=cpos[at][j]-ppos[at][j]
                            if np.fabs(diff) >0.5: #also covers the >1.0 differences.
                                shifts[at][j]+=(diff - np.sign(diff))
                            else: shifts[at][j]+=diff 
                    atoms[st].set_scaled_positions(init_pos+shifts)
                    ppos=dc(cpos)
         
        if args.unwrap and not args.ifAll: print ("Unwrapping of the atomic coordinates requires a multi-step/trajectory input"); exit()
        
        if args.wrap: #wrap coordinates back into cell
            if args.ifAll: 
               atoms_new=[]
               for i in range(len(atoms)):   x=dc(atoms[i]); x.wrap(pbc=1);atoms_new.append(x); 
               atoms=dc(atoms_new)
               #atoms=[wrap_positions(atom)
            else:atoms.wrap()

        if args.ifAll and (args.ifPrim or args.rep or args.interstices or args.airss or args.molecule or args.optim_xtb):  print ("Primitive cell, cell repetition/extension, adding interstices and elemental swapping, AIRSS input and LAMMPS molecule input and OPTIM+XTB input are not compatible with a multi-step/trajectory input"); exit()  #args.swaps or   


        if args.swaps:
            for sw in args.swaps:
                frm=sw.split(":")[0].split(",")
                to=sw.split(":")[1]
                if args.ifAll:
                    for i in range(len(atoms)):
                        for at in atoms[i]:
                            if at.symbol in frm:  atoms[i][at.index].symbol=to
                else:
                    for at in atoms:
                        if at.symbol in frm:  atoms[at.index].symbol=to



        if args.surface:
            atoms=surface(lattice=atoms,indices=args.surface,layers=args.layers,vacuum=args.vacuum, periodic=1)
            atoms=sort(atoms)

        if args.adsorb:
            #info = atoms.info.get('adsorbate_info', {})
            #print(info)
            max_zind = atoms.positions[:, 2].argmax()
            height = atoms.positions[max_zind, 2] + args.adsorb_height
            #print(atoms.get_tags())

            #Get layers from the structure based on z coordinates. #Use 0.7 A bin_width
            zcoords=[at.position[2] for at in atoms]

            try: 
                adatoms=molecule(args.adsorb)
                adatoms.rotate(180,'x',center='COP')
            except:
                try:adatoms=Atoms(args.adsorb)
                except:adatoms=ase.io.read(args.adsorb)

            if args.adsorb_target:
                """
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
                print(counts,bin_edges)
                """

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

                for key in sorted(layers.keys()):  print (key,layers[key])
                layers=[layers[key] for key in sorted(layers.keys())] #holds the atom IDs for each layer in the molecule
                print (layers)
    
                surf_atom_list=[]
                for lay in [layers[i] for i in args.adsorb_slayers]: #iterate over pre-defined surface layers (i.e. first two layer by default)
                    lay=atoms[lay]
                    surf_atom_list.extend([at for at in lay if at.symbol in args.adsorb_target])

                #print(surf_atom_list)
    

                if args.adsorb_count>len(surf_atom_list) and args.adsorb_nlayers==1:
                    print('Cannot add more adatoms than the surface atoms, will only add %d adsorbate atoms/molecules'%len(surf_atom_list))
                    args.adsorb_count=len(surf_atom_list)

                choices=random.sample(surf_atom_list,args.adsorb_count)
                for at in choices:
                    print(at)
                    ab=[at.position[0],at.position[1]]
                    print(ab)
                    add_adsorbate(atoms,adatoms,args.adsorb_height,ab,offset=None,mol_index=0)

            else: #if no specific target on the surface is given for adsorption
                if args.adsorb_ab:
                    atoms = add_molecule_to_surface(atoms, adatoms, ntimes=args.adsorb_count,nlayers=args.adsorb_nlayers, ab=args.adsorb_ab, tol=args.adsorb_tol,ori=args.adsorb_ori)
                else:
                    atoms = add_molecule_to_surface(atoms, adatoms, ntimes=args.adsorb_count, nlayers=args.adsorb_nlayers,random_ab=True, tol=args.adsorb_tol)

            atoms=sort(atoms)
            
        if args.min_tilt:
            minimize_tilt(atoms, order=range(0, 3), fold_atoms=True)

        if args.rep: 
            atoms=atoms.repeat(args.rep)
            atoms=sort(atoms)

        if args.zlim:
            print('Length before deletion: %d'%len(atoms))
            zlim=args.zlim
            sel=[at.index for at in atoms if (zlim[0]>at.position[2] or at.position[2]>zlim[1]) ]
            print('%d atoms to be deleted'%len(sel))
            del(atoms[sel])
            print('Length after deletion: %d'%len(atoms))
        
        #if args.coord and not args.ifAll:
        if args.coord: #TODO: This one is not working well, use the one in cryan_analyse.py
            cutoff=1.85 #Angstroem for determining the coord no: for C-C:1.85
            font = {'family': 'serif', 'size': 18}
            plt.rc('font', **font)
            fig = plt.figure(figsize=(11.69, 8.27))
            if args.ifAll:
                btype=[]
                if args.noSubs:subs=[at.index for at in atoms[0]]
                else:subs=[]
                
                for st in range(len(atoms)): #over time/geom step.
                    #if args.noSubs and st  in subs: continue #skip the substrate atoms
                    ats=atoms[st]
                    if st%50==0: print(st,end=', ');stdout.flush()
                    i = neighbor_list('i', ats, cutoff);  coord = np.bincount(i) #This does not work well!! Does not account for the MIC (or wrongly).
                    #Find alternative in ase.geometry or use ovito
                    
                    #"""
                    sp=0;sp2=0;sp3=0;oth=0
                    for at in range(len(coord)):
                        if at in subs: continue
                        #try:
                        if coord[at] == 2: sp+=1
                        elif coord[at] == 3: sp2+=1
                        elif coord[at] == 4: sp3+=1
                        elif coord[at] <2 or coord[at]>4: oth+=1
                        #except: print ('hata:', at)
                            
                    """
                    sp=np.count_nonzero(coord==2)
                    sp2=np.count_nonzero(coord==3)
                    sp3=np.count_nonzero(coord==4)
                    oth=np.count_nonzero((coord<=1) | (coord>4 ))
                    """
                    
                    tot=float(sp+sp2+sp3+oth)
                    #print(tot)
                    if tot==0: btype.append([st,0.,0.,0.,0.]) #st*args.skip
                    else:btype.append([st,sp3/tot,sp2/tot,sp/tot,oth/tot])
                #print (btype)
                
                with open('%s.data'%(inpf.split('.')[0]), 'w') as f:
                #with open('hybrid.data', 'w') as f:
                    f.write('#timeStep, sp3, sp2, sp, others\n')
                    for b in btype:
                        f.write('%8d %.3f %.3f %.3f %.3f \n'%(b[0],b[1],b[2],b[3],b[4]))
        
                btype=np.array(btype)
                plt.plot(btype.T[0],btype.T[1],label='sp3')
                plt.plot(btype.T[0],btype.T[2],label='sp2')
                plt.plot(btype.T[0],btype.T[3],label='sp')
                plt.plot(btype.T[0],btype.T[4],label='others')
                plt.legend()
                plt.ylim(0,1)
                plt.show()

            cnt+=1
            continue
                
        if args.interstices: #Adding the interstices first, before all other operations (prim, etc.)  is important.
            print('Adding the missing interstices as requested')
            atoms=find_interstices(atoms,tol=1e-6,sg=None,verb=None)#,prim=args.ifPrim)

        if args.ifPrim or args.kpath: 
            scaled_positions= atoms.get_scaled_positions()
            cell=(atoms.cell, scaled_positions, atoms.numbers)
            print ("%s: SG=%s"%(inpf,get_spacegroup(cell,symprec=args.tol)))
            if args.ifPrim:
                lattice, scaled_positions, numbers = find_primitive(cell, symprec=args.tol)
                atoms=Atoms(symbols=numbers,cell=lattice,scaled_positions=scaled_positions,pbc=True)
            if args.kpath:
                #from seekpath import *
                import seekpath
                path_data=seekpath.get_path(cell, with_time_reversal=False, symprec=args.tol)
                #print (path_data.keys())
                path=path_data['path']
                pc=path_data['point_coords']

                #compare SG from SPGlib and pathseek.

                #print (path,pc)
                #print
                print ("%BLOCK PHONON_FINE_KPOINT_PATH")
                for i,p in enumerate(path):
                    if i==0 or i==len(path):
                        if "'" not in p[0]:print ("%s %s %s #%s"%(pc[p[0]][0],pc[p[0]][1],pc[p[0]][2],p[0]))
                        if "'" not in p[1]: print ("%s %s %s #%s"%(pc[p[1]][0],pc[p[1]][1],pc[p[0]][2],p[1]))
                    else:
                        if "'" not in p[1]: print ("%s %s %s #%s"%(pc[p[1]][0],pc[p[1]][1],pc[p[1]][2],p[1]))

                print ("%ENDBLOCK PHONON_FINE_KPOINT_PATH")
                
                try:
                    from ase.dft.kpoints import get_special_points
                    path = get_special_points(atoms.cell,eps=args.tol)
                    print (path)
                except: raise Error; continue


        if args.airss:
            #Create AIRSS input
            seed=inpf.split('.')[0]
            with open(seed+'.airss','w') as f:
                scpos=atoms.get_scaled_positions()
                for ai,at in enumerate(atoms):
                    if ai==0:x='NUM=1 POSAMP=0 ANGAMP=0'
                    else:x=''
                    #f.write("%-2s %12.8f %12.8f %12.8f #%s %% %s \n"%(at.symbol,at.position[0],at.position[1],at.position[2],seed+"_"+str(len(atoms)),x))
                    f.write("%-2s %12.8f %12.8f %12.8f #%s %% %s \n"%(at.symbol,at.position[0],at.position[1],at.position[2],"group"+str(len(atoms)),x))
                    #f.write("%-2s %12.8f %12.8f %12.8f #%s %% %s \n"%(at.symbol,scpos[ai][0],scpos[ai][1],scpos[ai][2],"group"+str(len(atoms)),x))
            continue

        if args.optim_xtb:
            with open('mass','w') as f:
                for at in atoms: f.write("%s %.2f\n"%(at.symbol,at.mass))

            with open('coords','w') as f:
                for at in atoms: f.write('%15.8f %15.8f %15.8f\n'%(at.position[0],at.position[1],at.position[2]))
            
            atoms.write('input.turbo',format='turbomole')
            #TODO:replace $coord with $coord angs and add periodic info in the beginning.
            a,b,c,alpha,beta,gamma=atoms.get_cell_lengths_and_angles()
            str1="""$periodic 3
$cell angs
    %.5f %.5f %.5f %.1f %.1f %.1f
"""%(a,b,c,alpha,beta,gamma)
            
            with open('input.turbo','r') as f:
                lines=f.readlines()
            with open('input.turbo','w') as f:
                f.write(str1)
                for i in range(len(lines)):
                    if i==0: f.write('$coord angs\n')
                    else:f.write(lines[i])
            continue
                
        if args.molecule:
            seed=inpf.split('.')[0]
            #with open(seed+'.data','w') as f:
            with open('C%d.data'%len(atoms),'w') as f: #TODO: fix the key error?? Output is ocrrect.
                f.write('# C%d\n\n'%len(atoms))
                f.write('%d atoms\n\n'%len(atoms))
                types={};tcnt=0;x='';y=''
                for ai,at in enumerate(atoms):
                    if at.symbol not in types: tcnt+=1;types[at.symbol]=tcnt;
                    x+='%5d %5s\n'%(ai+1, types[at.symbol])
                    y+='%5d %12.8f %12.8f %12.8f\n'%(ai+1,at.position[0],at.position[1],at.position[2])
                f.write('Types\n\n'+x+"\n\nCoords\n\n"+y)
            continue
        
        if isinstance(atoms,list) and len(atoms) >1: 
           print ("%d timesteps/geometries were written to %s"%(len(atoms),outf)) #If 1 geometry in Atoms, then it would be an ase.Atoms instance.
           if args.separate:
                   for i,atom in enumerate(atoms):
                       x=outf.split('/')[-1].split('.')
                       if args.otype=="vasp":ase.io.write("%s-%s.%s"%(x[0],str(i),x[1]),atom,format=args.otype,vasp5=True)
                       else: ase.io.write("%s-%s.%s"%(x[0],str(i),x[1]),atom,format=args.otype)

        if args.otype=="vasp":ase.io.write(outf,atoms,format=args.otype,vasp5=True)
        elif args.otype=="lammps-data":ase.io.write(outf,atoms,format=args.otype,atom_style='charge')
        elif args.otype=='castep-cell':ase.io.write(outf,atoms,format=args.otype,positions_frac=1)
        #elif args.otype=='':
        elif args.otype=='res':  
            if args.itype=='vasp':
                Ep=atoms.get_total_energy()
                atoms.info['energy']=Ep
                atoms.info['name']=inpf
                P=atoms.info.get('pressure')
                if P is not None: atoms.info['pressure']=P
                else: atoms.info['pressure']=0.
                SG=atoms.info.get('spacegroup')
                if SG is None:   SG=str(get_spacegroup(atoms, symprec=args.tol).symbol)
                atoms.info['spacegroup']=SG
                atoms.info['times_found']=1

            #try:
            if 1:
                #ase.io.res.write_res(outf, atoms, write_info=0, write_results=1, significant_figures=6)
                myres=ase.io.res.Res(atoms)
                try: Ep=atoms.info['energy']
                except:
                  try:Ep=atoms.get_total_energy()
                  except:Ep=None
                myres.energy=Ep
                myres.write_file(outf, write_info=0,  significant_figures=6) #works!
            #except:ase.io.write(outf,atoms,format=args.otype)

        else: ase.io.write(outf,atoms,format=args.otype)


        if args.group=='no': system('mkdir -p Group-%d; mv %s Group-%d'%(len(atoms),outf,len(atoms)))


        if args.view: view(atoms)

        cnt += 1


"""
    except Exception as ex:  
        print ("Error with %s: %s"%(inpf, type(ex).__name__)), ex.args
        continue
"""

print ("Total of %d file(s) converted into %s format."%(cnt,args.otype))
print("Elapsed time: %.2f sec."%( time.time()-initT))



exit()


#BACKUPS
