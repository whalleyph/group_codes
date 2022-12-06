import glob
from matador.scrapers.castep_scrapers import res2dict
import pickle
import os

import argparse as ap

parser = ap.ArgumentParser()

parser.add_argument("-s", "--simtol", default = 0.1, type=float)

args = parser.parse_args()

if os.path.exists("datafile.pickle"):
	print("Found the pickle.")
	with open("datafile.pickle" ,'rb') as fi:
		polished_cursor = pickle.load(fi)
else:
	print("Can't find the pickle.")
	cursor, failures = res2dict("*.res", db=False)
	polished_cursor = filter_unique_structures(cursor, sim_tol=args.simtol, enforce_same_stoich=True, quiet=True)
	polished_cursor = [Crystal(standardize_doc_cell(doc)) for doc in polished_cursor]
	with open('datafile.pickle', 'wb') as fh:
		pickle.dump(polished_cursor, fh)

unique_res = [c.source[0] for c in polished_cursor]
os.system("mkdir unique")
for ur in unique_res:
	os.system(f"cp {ur} unique/")


