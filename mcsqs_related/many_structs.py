"""
cif2atat.py by Ben Saunders 2022

DESCRIPTION:
    This utility is for running multiple instances of mcsqs for multiple structures.
    It will determine how many strucutres you are using and assign the maximum number
    of threads to each structure equally, e.g.: with 48 cores and 9 structures, this
    will run 5 instances of mcsqs for each of the structures. 1 core will host this
    script, as the master, and the remaining 2 will remain idle.

    NOTE:
    It is probably best not to use this script for systems where resources are allocated
    on a financial basis, e.g. SULIS or ARCHER.

USAGE:
    python multi_structs.py --dirs <DIR1> <DIR2> etc. --cores

"""

import argparse as ap
import site
import numpy as np
from contextlib import contextmanager
from parallel_mcsqs import paralellise, mcsqs, run
import math
import concurrent.futures as cf
import os


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


def runInDir(dir, coresperdir, min_at, multi=1):
    with cd(dir):
        with open("multiple.of", "r") as f:
            m = int(f.readlines()[0])
        process_ids = [*range(0, coresperdir, 1)]
        sites = 0
        while sites < min_at:
            sites = m * multi
            multi = multi + 1

        o = paralellise(sites * 2, process_ids, coresperdir)
        # print(f"{dir}: {sites}")


def handleDirectories(dirs, totalcores, min_at):
    ndirs = len(dirs)
    coresperdir = int(math.floor(totalcores) / ndirs)
    print(coresperdir)
    arg_list = [[d, coresperdir, min_at] for d in dirs]
    with cf.ProcessPoolExecutor(max_workers=ndirs) as executor2:
        x = executor2.map(runInDir, *zip(*arg_list))
        # print(list(x))
    return x


parser = ap.ArgumentParser(
    description="Set multiple mcsqs sessions off at once, for multiple sctructures."
)

parser.add_argument(
    "-d",
    "--dirs",
    help="Directories to work through.",
    nargs="+",
)

parser.add_argument(
    "-c", "--cores", type=int, help="Number of cores available.", default=1
)
parser.add_argument(
    "-m",
    "--min_atoms",
    type=int,
    help="Multiplier to apply to the number of sites in the supercell, as written in multiple.of",
    default=2,
)
args = parser.parse_args()

print("MCSQS for many structures.")
l = handleDirectories(args.dirs, args.cores, args.min_atoms)
