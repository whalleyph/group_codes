#!/bin/env python3
from pymatgen.ext.optimade import *
import pymatgen.io.ase
import argparse
import os.path
import re
import ase.io
import ase.spacegroup  # .get_spacegroup(atoms, symprec=1e-05)

parser = argparse.ArgumentParser(
    description="Script for retrieving structures from the Materials Project database for a given stoichiometry/composition."
)


parser.add_argument(
    "-s",
    "--stoich",
    nargs="*",
    type=str,
    required=True,
    help="Chemical system(s) to query, A chemical system (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234). Multiple entries seperated by a white space can be used.",
)

parser.add_argument(
    "-ot", "--otype", type=str, required=False, help="Output file type", default="res"
)
parser.add_argument(
    "-o",
    "--outf",
    type=str,
    required=False,
    help="Output file name. Def: input name is used as root.",
)
parser.add_argument(
    "-od",
    "-odir",
    "--odir",
    type=str,
    required=False,
    default=None,
    help="Output file directory. Def: Coll_$stoich",
)

parser.add_argument(
    "-db",
    "--db_list",
    nargs="*",
    type=str,
    required=0,
    default=["oqmd", "aflow", "mp", "nmd", "mpds", "odbx"],
    help="Database list to get structures from, which support OPTIMADE framework. Check https://pymatgen.org/pymatgen.ext.optimade.html for up-to-date list",
)  # ,'mcloud.li-ion-conductors'


args = parser.parse_args()


ext = args.otype.split("-")[-1]
if args.otype == "pdb":
    args.otype = "proteindatabank"
elif args.otype == "extxyz":
    ext = "xyz"
if ext == "proteindatabank":
    ext = "pdb"

if args.odir:
    odir = args.odir
else:
    odir = "Coll_" + "+".join(args.stoich)
    os.system("mkdir -p %s" % odir)


data = {}
for db in args.db_list:
    print("\nCurrent database to extract structures from: ", db)
    with OptimadeRester(aliases_or_resource_urls=db, timeout=60) as m:
        # strs=m.get_structures(stoich, final=True)
        data = {}
        print(m.describe())
        for stoich in args.stoich:
            # data.extend(m.get_data(stoich,data_type=args.MP_dtype))#,prop='material_id')
            # Editted by Ben to retrieve single element queries
            if "-" in stoich:
                els = stoich.split("-")
                nels = len(els)
                form = None
                # If there is a number, or more than one upper case letter, it much be a formula, rather than an element query
            elif (
                any(ch.isdigit() for ch in stoich)
                or sum(1 for ch in stoich if ch.isupper()) > 1
            ):
                form = stoich
                els = None
                nels = None
            else:
                form = None
                els = [stoich]
                nels = len(els)

            print(stoich, form, els, nels)
            try:
                data.update(
                    m.get_structures(
                        elements=els, nelements=nels, chemical_formula_hill=form
                    )
                )
            except:
                None
            # print(m.get_structures(elements=['Li','S'],nelements=2,chemical_formula_hill=None))

        try:
            key = list(data.keys())[0]
            print(key)
            print(
                "%d structure(s) were found for %s"
                % (len(data[key].values()), ", ".join(args.stoich))
            )
        except:
            print("No structure found on ", db)
            continue
        # print(data.keys())

        for dt in data[key]:
            str1 = data[key][dt]
            atoms = pymatgen.io.ase.AseAtomsAdaptor.get_atoms(str1)
            mid = str(dt)
            form = atoms.get_chemical_formula(mode="hill", empirical=1)
            sg = ase.spacegroup.get_spacegroup(atoms, symprec=1e-03).symbol
            sg = sg.replace(
                "/", "").replace("(", "").replace(")", "").replace(" ", "")

            if db == "mp":
                outf = odir + "/" + form + "_" + sg + "_" + mid + ".res"
            else:
                outf = odir + "/" + form + "_" + sg + "_" + db + "-" + mid + ".res"

            if args.otype == "vasp":
                ase.io.write(outf, atoms, format=args.otype,
                             vasp5=True, append=0)
            elif args.otype == "lammps-data":
                ase.io.write(
                    outf, atoms, format=args.otype, atom_style="charge", append=0
                )
            else:
                ase.io.write(outf, atoms, format=args.otype, append=0)


# atoms_pmg=pymatgen.io.ase.AseAtomsAdaptor.get_structure(atoms)
