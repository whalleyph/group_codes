#!/bin/env python
import numpy as np
import optimade
from optimade.client import OptimadeClient
import argparse, os.path, re
import ase.io
import ase.spacegroup  # .get_spacegroup(atoms, symprec=1e-05)
from ase import Atoms

parser = argparse.ArgumentParser(description='Script for retrieving structures from the optimade databases for a given stoichiometry/composition.')

parser.add_argument(
    "-s",
    "--stoich",
    nargs="*",
    type=str,
    required=True,
    help="Specfic stoichiometry to query e.g. Na2O2",
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

args = parser.parse_args()

if args.odir:
    odir = args.odir
else:
    odir = "Coll_" + "+".join(args.stoich)
    os.system("mkdir -p %s" % odir)

my_base_urls = ["http://oqmd.org/optimade/", "https://optimade.materialsproject.org", "https://nomad-lab.eu/prod/rae/optimade/"] #3 databases that are known for crystal structures
m = OptimadeClient(base_urls = my_base_urls)
sql = []
if "stoich" in args:                                 #specific grammar is required for optimade module and the code below fixes the argument to that       
        sql = f'{args.stoich}'
        disallowed_characters = "[]"
        for character in disallowed_characters:
            sql = sql.replace(character, '')
        disallowed_characters = "'"
        for character in disallowed_characters:
            sql = sql.replace(character, '"')

filter = f'chemical_formula_reduced={sql}'
m.get(filter)
results = m.all_results

for url in my_base_urls:                             #Results are received in a dictionaries within dictionaries. Code below sorts through and outputs in to a file in the .res format
    data = results["structures"][f'chemical_formula_reduced={sql}'][url].data
    if len(data) > 0:
        for struct in data:
            coords = (struct["attributes"]["cartesian_site_positions"])
            cell_params = struct["attributes"]["lattice_vectors"]
            els = (struct["attributes"]["species_at_sites"])
            id = (struct["id"])
            form = (struct["attributes"]["chemical_formula_reduced"])
            atoms = Atoms(els, positions = coords, cell = cell_params)
            frac_coords = atoms.get_scaled_positions()
            frac_atoms = Atoms(els, positions = frac_coords)
            outf = odir + "/" + form + '_' + str(id) + '.res'
            ase.io.write(outf, frac_atoms, format = args.otype, append = 0)
