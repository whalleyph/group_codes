#!/bin/env python3
import argparse
import os.path
import re

import ase.io
import pymatgen.io.ase
from pymatgen.ext.matproj import *

parser = argparse.ArgumentParser(description='Script for retrieving structures from the Materials Project database for a given stoichiometry/composition.')


parser.add_argument('-s','--stoich', nargs='*', type=str,required=True, help='Chemical system(s) to query, A chemical system (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234). Multiple entries seperated by a white space can be used.')

parser.add_argument('-ot','--otype', type=str,required=False, help='Output file type',default='res')
parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: input name is used as root.')
parser.add_argument('-od','-odir','--odir', type=str,required=False, default=None,help='Output file directory. Def: Coll_$stoich')

parser.add_argument('-key','--APIkey', type=str,required=0, default='ZwyRz2p5d89sNxSH', help='API key to access MP database (you can access it through the MP website account dashboard ')
parser.add_argument('-dt','--dtype', type=str,required=False, help='Data type to query from the database, def: VASP/computational (also include experimental ones), choices: "vasp","exp")',default='vasp')

args = parser.parse_args()


ext=args.otype.split("-")[-1]
if args.otype=='pdb':args.otype='proteindatabank'
elif args.otype=='extxyz':ext='xyz'
if ext=='proteindatabank':ext='pdb'

if args.odir:odir=args.odir
else:odir="Coll_"+"+".join(args.stoich);os.system('mkdir -p %s'%odir)


with MPRester(args.APIkey) as m:
  #strs=m.get_structures(stoich, final=True)
  data=[]
  for stoich in args.stoich:
    data.extend(m.get_data(stoich,data_type=args.dtype))#,prop='material_id')
  print('%d structure(s) were found for %s'%(len(data),", ".join(args.stoich)))
  for dt in data:
    mid=dt['material_id']; form=dt['pretty_formula'];sg=dt['spacegroup']['symbol'].replace('/','').replace('(','').replace(')','')
    str1=m.get_structure_by_material_id(mid, final=True, conventional_unit_cell=1)
    atoms=pymatgen.io.ase.AseAtomsAdaptor.get_atoms(str1)


    outf=odir+"/"+mid+'_'+form+'_'+sg+'.res'

    if args.otype=="vasp":ase.io.write(outf,atoms,format=args.otype,vasp5=True,append=0)
    elif args.otype=="lammps-data":ase.io.write(outf,atoms,format=args.otype,atom_style='charge',append=0)
    else: ase.io.write(outf,atoms,format=args.otype,append=0)

  
#atoms_pmg=pymatgen.io.ase.AseAtomsAdaptor.get_structure(atoms) 
