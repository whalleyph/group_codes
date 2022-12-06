"""
doping4airss.py - Benedict Saunders, 2022

    The script prepares a .cell file that can be used to lithiate a structure
    using the AIRSS program.
"""
import argparse
import re
from pathlib import Path
from pprint import pprint
from typing import *

import numpy as np
from ase.io import read
from chemparse import parse_formula as Formula
from spglib import get_spacegroup

from input4airss import pyairss_input, vprint


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def printc(s, bcolor, v=True):
    vprint(bcolor + str(s) + bcolors.ENDC)


def getSGnum(atoms):
    spg = get_spacegroup(
        (atoms.cell, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-2
    ).split(" ")
    spgd = dict(
        zip(
            [
                "symbol",
                "number",
            ],
            [
                spg[0],
                int(re.sub(r"[()]", "", spg[1])),
            ],
        )
    )
    printc(
        f"Using SpaceGroup {spgd['number']} ({spgd['symbol']})", bcolor=bcolors.OKCYAN
    )
    return spgd


class doping(pyairss_input):
    def __init__(self, struct, enforceSG=False, verbose=True):
        super().__init__(verbose)
        self.verbose = verbose
        self.atoms = read(struct, index=-1)
        self.geomfile = struct
        self.cell_matrix = list(self.atoms.__dict__["_cellobj"])
        self.cell_params = self.atoms.cell.cellpar()
        self.elements = Formula(self.atoms.get_chemical_formula())
        self.SG = enforceSG

    def set_dopants(self, dopants, radius_scaling=1.00):

        ### Dopants must be given in a dictionary: {<key=symbol>:[<lower range>, <upper range>]}
        ### If only one number is specified, it still needs to be a single-element list!

        for dopant in list(dopants.keys()):
            if dopant not in self.elements.keys():
                self.elements[dopant] = 0
        self.get_radii(scaling=radius_scaling)
        self.dopants = dopants

    def make_cell_file(self, seed, vectors=False):
        uid = Path(self.geomfile).stem
        name = f"{uid}-{seed}.cell"
        printc(f"\nWriting input file {name}", bcolor=bcolors.OKCYAN, v=self.verbose)
        self.f = open(name, "w")
        self.f.write(
            f"## Doping {self.geomfile} input for PyAIRSS - Made with doping4airss.py\n"
        )
        self.f.write(f"## Benedict Saunders, 2022\n\n")

        ### Usage of spacegroup
        if self.SG:
            num = str(getSGnum(self.atoms)["number"])
        else:
            printc("Will use random spacegroup...\n", bcolor=bcolors.WARNING)
            num = "10 to 230"
        self.f.write(f"#! define s1 space group number={num}\n\n")

        ### Writing the lattice parameters
        if vectors:
            self.f.write("%BLOCK LATTICE_CART\n")

            for line in self.cell_matrix:
                self.f.write(f"{line[0]}\tl{line[1]}\t{line[2]}\n")
            self.f.write("%ENDBLOCK LATTICE_CART\n")
        else:
            self.f.write("%BLOCK LATTICE_ABC\n")
            lcp = list(np.around(self.cell_params, decimals=4).astype(str))
            lens = "\t".join(lcp[:3])
            angs = "\t".join(lcp[3:])

            self.f.write(f"{lens}\n{angs}\n")
            self.f.write("%ENDBLOCK LATTICE_ABC\n")
        self.f.write("\n")

    def write_position_block(self, fractional=False, precision=6):
        if not fractional:
            positions = self.atoms.get_positions()
            title = "POSITIONS_ABS"
        else:
            positions = self.atoms.get_scaled_positions()
            title = "POSITIONS_FRAC"
        elems = self.atoms.get_chemical_symbols()

        ### This stuff has been taken from the ASE source, as the ASE package seems to
        ### require some instance of a CASTEP calculator to convert an atoms object
        ### into a .cell file.

        fformat = "%{0}.{1}f".format(precision + 3, precision)
        cell_block_format = " ".join([fformat] * 3)
        pos_block_format = "%s " + cell_block_format

        self.f.write(f"%BLOCK {title}\n")
        for idx, elem in enumerate(elems):
            xyz = positions[idx]
            line = pos_block_format % tuple([elem] + list(xyz))
            self.f.write(f"{line}\n")
        self.f.write(f"%ENDBLOCK {title}\n\n")

        ### This concludes the ASE plagarism :)

        ### Defining the separations
        self.write_seperations()

    def randomise_dopants(self):
        ### Defining the randomness of the dopants
        for dopant in list(self.dopants.keys()):
            rng = self.dopants[dopant]
            if len(rng) > 1:
                r = f"{rng[0]} to {rng[1]}"
            else:
                r = rng[0]
            self.f.write(f"#! define n{dopant} int {r}\n")
        self.f.write("\n")

        ### Placing the dopants into the lettuce
        pairs_list = ["".join(c) for c in self.combinations]
        pairs = ", ".join(pairs_list)
        for dopant in list(self.dopants.keys()):
            self.f.write(
                f"#! insert and randomise n{dopant} {dopant} atoms into cell subject to s1, {pairs}\n"
            )
        self.f.write("\n")


def handle_dopant_input(dl):
    d = " ".join(dl)
    dops = {}
    dummy = "___"
    for x in d.split(" "):
        if x.replace(".", "", 1).isdigit():
            try:
                dops[dummy].append(int(x))
            except:
                printc("The number of dopants must be an integer!", bcolor=bcolors.FAIL)
            dops[dummy].sort()
        else:
            dops[x] = []
            dummy = x
    return dops


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Doping4AIRSS")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        nargs="*",
        help="The input files to be doped with things.",
    )
    parser.add_argument(
        "-s",
        "--seed",
        default="doped",
        required=False,
        help="Output file will be named <Input filr>-<seed>.cell. Def.: <Input file>-doped.cell",
    )
    parser.add_argument(
        "-d",
        "--dopants",
        required=True,
        help="The specie(s) to dope the lettuce with. Eg. -d Li 4 Na 3 5 will dope with 4 Li atoms, and between 3 and 5 Na atoms.",
        nargs="*",
    )
    parser.add_argument(
        "-spg",
        "--spg",
        action="store_true",
        default=False,
        help="Enforce the spacegroup of the input structure for the doping procedure.",
    )
    args = parser.parse_args()

    ### Would nice to have these as command line arguments...

    dop_dict = handle_dopant_input(args.dopants)

    for geom in args.input:
        lithiator = doping(geom, args.spg)
        lithiator.set_dopants(dopants=dop_dict, radius_scaling=0.85)
        lithiator.make_cell_file(args.seed)
        lithiator.write_position_block()
        lithiator.randomise_dopants()
        lithiator.finalise(save_as_opt=True)
