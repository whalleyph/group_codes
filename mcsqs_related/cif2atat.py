"""
cif2atat.py by Ben Saunders 2022

DESCRIPTION:
    This utility converts a CIF to the format required by the ATAT toolkit, specifically as an input
    for ATAT's corrdump, which generates the inputs used by the MCSQS (mcsqs) executable.
    For each CIF processed, a directory is made where the output, 'rndstr.in', is saved. Useful when
    it comes to running the tools from ATAT.

    If a CIF is given which does not specify any fractional occupancies, nothing will be done, and it
    will be skipped.

USAGE:
    python cif2atat.py --files <CIF1> <CIF2> etc.

"""

import site
import numpy as np
from ase.io.cif import read_cif
from ase.io.vasp import write_vasp
from contextlib import contextmanager
import os
import argparse
import math
from fractions import Fraction
from ase.build import make_supercell
import primefac as pf


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


def printc(s, bcolor):
    print(bcolor + s + bcolors.ENDC)


@contextmanager
def cd(newdir):
    if not os.path.isdir(newdir):
        os.mkdir(newdir)
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def dec2frac(dec):
    if dec == 0.333:
        frac = (1, 3)
    elif dec == 0.667:
        frac = (2, 3)
    else:
        frac = [int(a) for a in str(Fraction(dec).limit_denominator(100)).split("/")]
    return frac


def approximateEqualFactors(x, n=3):
    prime_factors = list(pf.primefac(x))  # from the primefac module
    if len(prime_factors) == n:
        F = prime_factors
    elif len(prime_factors) < n:
        while len(prime_factors) < n:
            prime_factors.append(1)
        F = prime_factors
    else:
        nth_root = x ** (1.0 / float(n))
        facs = np.ones(n, dtype=int)
        prime_factors.sort(reverse=True)
        for f in prime_factors:
            for idx, i in enumerate(facs):
                test = f * i
                if test < nth_root:
                    facs[idx] = test
                    break
                if idx == n - 1:
                    minimum = np.min(facs)
                    minloc = np.argmin(facs)
                    facs[minloc] = minimum * f
        F = facs
    P = np.zeros((n, n), dtype=int)
    np.fill_diagonal(P, F)
    return F, P


def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)


def flatten(l):
    return [item for sublist in l for item in sublist]


def unpackNestedDict(it):
    if isinstance(it, list):
        for sub_it in it:
            yield from unpackNestedDict(sub_it)
    elif isinstance(it, dict):
        for value in it.values():
            yield from unpackNestedDict(value)
    else:
        yield it


def defineSupercell(atoms):
    """
    Using the occupancies in the CIF, this function calculations the number of sites required in
    the final supercell such that every site will have an integer occupancy.
    """
    sites = len(atoms)
    lst_occs = list(atoms.info["occupancy"].items())
    indexed_partials = []
    unique_occupancies = []
    for idx, occ in enumerate(lst_occs):
        d2l = set(occ[1].items())
        if len(d2l) > 1:
            indexed_partials.append(occ)
            if d2l not in unique_occupancies:
                unique_occupancies.append(d2l)

    non_int_vals = [x[1] for x in flatten(unique_occupancies)]

    i = 1
    if len(non_int_vals) > 0:
        np_nums = np.array(non_int_vals)
        while True:
            new = np_nums * i
            check = all(x.is_integer() for x in new)
            if check:
                break
            else:
                i = i + 1
    F, P = approximateEqualFactors(i)
    return i, P, sites, len(non_int_vals), False


def check_cell(m):
    mx = list(m)
    if len(mx) == 3 and len(mx[1]) == 3:
        return mx
    else:
        new = np.zeros((3, 3), dtype=float)
        np.fill_diagonal(new, mx)
        return new


def convertCIF(
    cif,
    do_corrdump,
    multibodyterms=[0, 0, 0],
    scaled=True,
    write_matrix=True,
    prim=False,
    SC=27,
):
    """
    Reading the atoms object into a dictionary, this function extractions occupancies and information
    about the cell as required by the input specification for ATAT.
    """
    print(f"\nWorking {cif}...")
    st = (
        not prim
    )  ### For whatever reason, subtrans_included=False causes an issue with OMP_NUM_THREADS. Why? Who knows. Solution? Don't use it.
    atoms = list(
        read_cif(
            cif,
            index=slice(None),
            subtrans_included=st,
            primitive_cell=prim,
        )
    )[-1]
    atoms_dict = atoms.__dict__

    if scaled:
        positions = atoms.get_scaled_positions(wrap=False)
    else:
        positions = atoms_dict["arrays"]["positions"]

    SG_kinds = atoms_dict["arrays"]["spacegroup_kinds"]
    occs = atoms_dict["info"]["occupancy"]

    ### Defining th supercell multiplier/tranformation matrix

    mult, multSCmatrix, unit_cell_sites, num_of_int, fail = defineSupercell(atoms)
    if SC != 0:
        mult = SC
        _, multSCmatrix = approximateEqualFactors(SC)

    if num_of_int == 0 or fail == True:
        printc("No fractional occupancies present. Skipping...", bcolor=bcolors.WARNING)
        return 1, None
    tot_atoms = mult * unit_cell_sites

    if tot_atoms >= 100:
        printc(
            f"    >>> WARNING: Supercell will have {tot_atoms} atoms <<<",
            bcolor=bcolors.WARNING,
        )

    # I don't know which one of these it wants! I think the second matrix or second to fourth line is the transformation matrix

    idmtx = np.identity(3)
    # cell_matrix = list(atoms_dict["_cellobj"])
    cell_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    abc = atoms.cell.cellpar()
    cart_cell = atoms.get_cell(complete=True)
    cart_matrix = check_cell(cart_cell)

    ### Making the rndstr.in file ###

    name = cif[:-4]
    with cd(name):
        os.system(f"cp ../{cif} .")
        write_vasp("POSCAR", atoms=atoms)
        with open("rndstr.in", "w") as f:
            f.write(" ".join([str(i) for i in abc]) + "\n")

            # Writing the cell (transformation?) matrix
            if write_matrix:
                for ln in cell_matrix:
                    lns = [str(x) for x in ln]
                    f.write(" ".join(lns) + "\n")

            # Writing the positions and occupancies
            for idx, pos in enumerate(positions):
                coords = " ".join([str(i) for i in pos])
                site_occ_dict = occs[str(SG_kinds[idx])].items()
                if len(site_occ_dict) > 1:
                    site_info = ",".join([f"{i[0]}={str(i[1])}" for i in site_occ_dict])
                else:
                    site_info = list(site_occ_dict)[0][0]
                f.write(coords + " " + site_info + "\n")

        ln = len(positions)
        printc(f'  File "rndstr.in" generated from "{cif}".', bcolors.OKGREEN)

        ### Writing sqscell.out which constrains the output of MCSQS

        super = make_supercell(atoms, P=multSCmatrix)
        super_matrix = check_cell(super.get_cell())
        with open("sqscell.out", "w") as sf:
            sf.write("1\n\n")
            for line in super_matrix:
                strline = [f"{x:.5f}" for x in line]
                joined = "\t".join(strline)
                sf.write(f"{joined}\n")

        ### Running corrdump in preparation for mcsqs
        if do_corrdump:
            corrdump_command = f"corrdump -ro -noe -nop -l=rndstr.in -2={multibodyterms[0]} -3={multibodyterms[1]} -4={multibodyterms[2]} -clus > corr.out"
            try:
                os.system(corrdump_command)
            except:
                printc("  >>> FAILURE: corrdump failed <<<", bcolors.FAIL)

    return 0, name


def checkInputs(args):
    """
    This function checks that, if corrdump is to be run, at least one of the multibody terms are non-zero.
    """

    multibody = [args.two, args.three, args.four]
    if args.corrdump and not any(i > 0 for i in multibody):
        printc(
            ">>> FAILURE: You must define at least one multibody term to be greater than 0. Nothing was done. <<<",
            bcolor=bcolors.FAIL,
        )
        exit(1)
    return multibody


### MAIN ###
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert CIFs into the input required by the ATAT toolkit, specifcally corrdump for mcsqs."
    )
    parser.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="CIFs to be converted. New directory will be created for each.",
    )
    parser.add_argument(
        "-c",
        "--corrdump",
        action="store_true",
        default=False,
        help="Run the corrdump utility for MCSQS preparation",
    )
    parser.add_argument(
        "-2",
        "--two",
        default=0,
        type=int,
        help="Two-body cutoff for generating correlation functions.",
    )
    parser.add_argument(
        "-3",
        "--three",
        default=0,
        type=int,
        help="Three-body cutoff for generating correlation functions.",
    )
    parser.add_argument(
        "-4",
        "--four",
        default=0,
        type=int,
        help="Four-body cutoff for generating correlation functions.",
    )
    parser.add_argument(
        "-s",
        "--fractional",
        action="store_true",
        default=False,
        help="Use fractional coordinate system in MCSQS input file.",
    )
    parser.add_argument(
        "-p",
        "--primitive",
        action="store_true",
        default=False,
        help="I cannot remeber what this does...",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        type=int,
        default=27,
        help="Number of cells to be in supercell. Use 0 to use maximum intger occupancies. Def.: 27, i.e 3x3x3",
    )
    args = parser.parse_args()

    files = args.files

    successful = []
    multibody = checkInputs(args=args)

    for cif in files:
        exit_code, name = convertCIF(
            cif,
            args.corrdump,
            scaled=args.fractional,
            multibodyterms=multibody,
            prim=args.primitive,
            SC=args.matrix,
        )
        if exit_code == 0:
            successful.append(name)
