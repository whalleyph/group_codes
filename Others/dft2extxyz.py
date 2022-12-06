#!/usr/bin/env python3

"""
castep2extxyz.py Benedict Saunders, 2022

A method to extract geometries, energies and forces from the output of
CASTEP runs.

TODO: tqdm for progress bay and parallelisation
"""
import pyfiglet
import argparse

import random
import string

import traceback
import os
from pathlib import Path
from tqdm import tqdm
import numpy as np
from ase.io import read as read_atoms
from ase.io.extxyz import write_extxyz, write_xyz
from matscipy import elasticity
import matplotlib.pyplot as plt


class bcolors:	### ooooo pretty colours
	HEADER = "\033[95m"
	OKBLUE = "\033[94m"
	OKCYAN = "\033[96m"
	OKGREEN = "\033[92m"
	WARNING = "\033[93m"
	FAIL = "\033[91m"
	ENDC = "\033[0m"
	BOLD = "\033[1m"
	UNDERLINE = "\033[4m"


def get_random_string(used):
	chars = f"{string.ascii_uppercase}{string.ascii_lowercase}{string.digits}"
	rs = "".join(random.choice(chars) for i in range(8))
	if rs not in used:
		used.append(rs)
		return rs, used
	else:
		get_random_string(used)


def get_stems(inputs, use_parent=True):
	stems = [Path(x).stem.split(".")[0] for x in inputs]
	if len(stems) != len(set(stems)):
		if use_parent:
			new_stems = [f"{Path(input).parent}_{Path(input).stem}" for input in inputs]
		else:
			used = []
			new_stems = []
			for idx, stem in enumerate(stems):
				rand_string, used = get_random_string(used=used)
				new_stems.append(f"{stem}-{rand_string}")
	else:
		new_stems = stems
	return new_stems


def printc(s, bcolor):
	tqdm.write(bcolor + str(s) + bcolors.ENDC)


def flatten(l):
	return [el for sl in l for el in sl]


def process_atoms(atoms_list, fqhl, fl):
	idc = []
	l = len(atoms_list)
	# Getting the first, quarter-point, halfway point and end point of the trajectories
	if fqhl and not fl:
		if l < 4:
			return atoms_list
		else:
			q = int(l / 4) - 1
			h = int(l / 2) - 1
			processed = [atoms_list[0], atoms_list[q], atoms_list[h], atoms_list[-1]]
			return processed
	# Getting the first and last points in the trajectories
	elif fl and not fqhl:
		if l < 2:
			return atoms_list
		else:
			processed = [atoms_list[0], atoms_list[-1]]
			return processed
	else:
		return atoms_list


def get_virial(atoms):
	stress = elasticity.Voigt_6_to_full_3x3_stress(atoms.get_stress())
	virial_gap = stress * atoms.get_volume() * -1
	virial_flat = [str(a) for a in list(virial_gap.flatten())]
	return virial_gap, " ".join(virial_flat)


def get_stats(nsymbols):
	stat_dict = {}
	symbols = flatten(nsymbols)
	unique = list(set(symbols))
	tot = len(symbols)
	for sym in unique:
		stat_dict[sym] = 100 * symbols.count(sym) / tot
	return stat_dict


def warn_and_quit(s, skip=True):
	if s in ["energy", "forces", "stess tensor"]:
		printc(f"Cannot find the {s}. Quitting...", bcolor=bcolors.FAIL)
	else:
		printc("Failure:", bcolor=bcolors.FAIL)
		tqdm.write(s)
		tqdm.write(traceback.format_exc())
	if not skip:
		exit(1)
	else:
		pass

title = pyfiglet.figlet_format("DFT2XYZ", font="larry3d")
print(title)
print()
parser = argparse.ArgumentParser(
	description="Convert CASTEP/VASP output to extxyz format for GAP-SOAP."
)

parser.add_argument(
	"-i",
	"--inputs",
	help="Input files",
	nargs="+",
	required=True,
)

parser.add_argument(
	"-a",
	"--all",
	help="Use all steps in the input file, not just the final geometry. Files will be saves as <stem>-step<n>.xyz",
	action="store_true",
	required=False,
)

parser.add_argument(
	"-rnd",
	"--random",
	help="use a random 8 charater string as the UUID, rather than the parent directory name ",
	action="store_true",
	required=False,
)


parser.add_argument(
	"-e",
	"--energy",
	help="Include the DFT Energy in the extended xyz file.",
	action="store_true",
	required=False,
)

parser.add_argument(
	"-f",
	"--forces",
	help="Include the DFT Forces in the extended xyz file.",
	action="store_true",
	required=False,
)

parser.add_argument(
	"-v",
	"--virial",
	help="Include the DFT Virial in the extended xyz file.",
	action="store_true",
	required=False,
)

parser.add_argument(
	"-hs",
	"--hessian",
	help="Include the Hessian in the extended xyz file.",
	action="store_true",
	required=False,
)

parser.add_argument(
	"-vvv",
	"--superverbose",
	help="Print ALL the steps!",
	action="store_true",
)

parser.add_argument(
	"-t",
	"--test_size",
	help="Testing set size for GAP. Def.: 200",
	required=False,
	default=200,
)

parser.add_argument(
	"-fqhl",
	"--fqhl",
	help="Get the first, last, middle and (last/4)th ionic steps",
	required=False,
	action="store_true",
)

parser.add_argument(
	"-fl",
	"--firstlast",
	help="Use the first and last geometries in the trajectory.",
	required=False,
	action="store_true",
)

args = parser.parse_args()
fqhl = args.fqhl
fl = args.firstlast
inputs = args.inputs
parent_uuid = not args.random
symbols = []
options = {"E": args.energy, "F": args.forces, "V": args.virial, "H": args.hessian}
if all(value == False for value in options.values()) == True:
	warn_and_quit(
		"You must include at least 1 property: One or more of -e, -f, -v, -hs"
	)

if all(v == True for v in [fl, fqhl]) == True:
	warn_and_quit(
		"You must chose either the first and last, or the first/quarter/half/last steps.\nYou cannot have both!"
	)

# Set the index or slice to retrieve with ase.io.read()
if any([args.all, fqhl, fl]):
	atom_index = ":"
else:
	atom_index = -1

# Determine if the stems are unique, and if not, format them such that they are (with UUIDs or parent directory etc.)
new_stems = get_stems(inputs, use_parent=parent_uuid)
failed = []
cols = ["symbols", "positions"]
for idx, _file in enumerate(tqdm(inputs)):
	tqdm.write(f"Converting {_file}...")
	stem = new_stems[idx]
	try:

		# atoms = read_castep(castep_file)
		raw_atoms_list = read_atoms(_file, index=atom_index)

		# If only the final frame is read, then atoms_list will be an atom object. We
		# need to convert is to something iterable, hence the square braces.
		if atom_index == -1:
			raw_atoms_list = [raw_atoms_list]
			multiple = False
		elif fqhl == False:
			multiple = True
		atoms_list = process_atoms(raw_atoms_list, fqhl, fl)
		tqdm.write(f"  Structures: {len(atoms_list)}")

		energies = []
		updated_atoms = []
		idc = []
		# There's probably a better, more idiomatic way to do this. But hey, it works.
		for idx, atoms in enumerate(atoms_list):
			idc.append(idx)
			if args.superverbose:
				if multiple:
					tqdm.write(f"	 Step {idx+1}")
				else:
					tqdm.write(f"	 Final ionic step")

			### THIS IS THE BIT WE ALL CARE ABOUT!!! ###
			if options["E"]:
				try:
					energy = atoms.get_potential_energy(force_consistent=True)
					atoms.info["energy"] = energy
					energies.append(energy)
				except:
					warn_and_quit("energy")

			if options["F"]:
				try:
					forces = atoms.get_forces()
					# atoms.info["forces"] = forces
				except:
					warn_and_quit("forces")

			if options["V"]:
				try:
					mtx, V_str = get_virial(atoms)
					atoms.info["virial"] = V_str
				except:
					warn_and_quit("stress tensor/virial")

			if options["H"]:
				raise NotImplementedError("Hessian extraction not implemented.")

			updated_atoms.append(atoms)

			# Fitting the energy minimisation curve
			# max_extraction_count = 5
			# steps = len(energies)
			# max_idx = steps - 1
			# to_get = []
			# for x in range(max_extraction_count):
			#	  div = 1
			#	  to_get.append = max_idx/div
			#	  div += 1

			# Writing the file, adding the step index if multiple are used.
			if fqhl:
				name = f"{stem}_sample{idx+1}".split("/")[-1]
			elif multiple:
				name = f"{stem}_step{idx+1}".split("/")[-1]
			else:
				name = stem.split("/")[-1]
			symbols.append(atoms.get_chemical_symbols())
			# Check whether all ionioc steps are selected, or just the last one.

			with open(f"{name}.xyz", "w") as fileobj:
				write_extxyz(
					fileobj,
					atoms,
					write_info=True,
					columns=cols,
					write_results=options["F"],
				)

		# plt.scatter(x=idc, y=energies)
		# plt.show()

	# This prints out the entire traceback if an error if caught (typically with ase.io.read)
	# and will either quit or skip to the next frame.
	except Exception as e:
		warn_and_quit(s=traceback.format_exc(), skip=True)
		failed.append(_file)
	printc(f"Done", bcolor=bcolors.OKGREEN)
printc("=== Complete ===\n", bcolor=bcolors.OKGREEN + bcolors.BOLD)
d = get_stats(symbols)
# print(symbols)
printc(f"Atom occurences", bcolor=bcolors.HEADER + bcolors.BOLD)
for key, val in d.items():
	tqdm.write(f"  {key}\t{val:.3f} %")
print()
printc(f"{len(failed)} Failures", bcolor=bcolors.FAIL + bcolors.BOLD)
for fail in failed:
	printc(f"  {fail}", bcolor=bcolors.WARNING)
# if args.test_size is not None:
#	 printc(
#		 "\n Splitting into training and testing sets for GAP...", bcolor=bcolors.BOLD
#	 )
#	 os.system("mkdir test train")
#	 os.system("mv ISOLATED_* train")
#	 os.system(f"ls *.xyz | shuf -n {int(args.test_size)} | xargs -i mv {{}} test")
#	 os.system("mv *.xyz train")
#	 os.system("cat train/*.xyz > train.xyz")
#	 os.system("cat test/*.xyz > validate.xyz")
#	 os.system("tar -czf all_xyz.tar.gz *.xyz")
"""


gap_fit 
	energy_parameter_name=energy 
	force_parameter_name=forces 
	do_copy_at_file=F 
	sparse_separate_file=T 
	gp_file=GAP_3b.xml 
	at_file=train.xyz 
	default_sigma={
		0.008	#ENERGY
		0.04	#FORCES
		0		#VIRIAL
		0		#HESSIAN
	}
	gap={
		distance_2b
		cutoff=4.0
		covariance_type=ard_se
		delta=0.5
		theta_uniform=1.0
		sparse_method=uniform
		add_species=T
		n_sparse=10
		:
		angle_3b
		cutoff=3.5
		covariance_type=ard_se
		delta=0.5
		theta_fac=0.5
		add_species=T
		n_sparse=30
		sparse_method=uniform
	}


"""
