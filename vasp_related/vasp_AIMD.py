"""     AIMD SETUP FOR VASP - Benedict Saunders, 2022
    This scripts sets up and submits the calculations for AIMD runs over various
    temperatures.

"""

import argparse as ap
import os
from contextlib import contextmanager
from copy import deepcopy as dcp

from pathlib import Path

import numpy as np
from ase.build import make_supercell
from ase.build.tools import sort as sort_atoms
from ase.calculators.vasp import Vasp
from ase.io import read as ase_read
from ase.io.vasp import read_vasp, write_vasp
from chemparse import parse_formula as Formula

from vasp4uspex import handle_hubbard, handle_magmoms


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


def handle_input_combinations(thermostat, ens, inputs_dict):
    inputs = {}
    available_ensembles = {
        "NVE": 1,
        "NVT": 2,
        "NpT": 3,
        "NpH": 4,
    }
    available_thermostats = {
        "Andersen": 1,
        "Nose-Hoover": 2,
        "Langevin": 3,
        "Multiple Andersen": 13,
    }

    try:
        ensemble_idx = available_ensembles[ens]
    except:
        print(f"Ensemble '{ens}' not recognised. Exitting...")
        exit()
    inputs["isif"] = ensemble_idx
    if ensemble_idx == 1:
        thermostat = 1
        inputs["andersen_prob"] = 0

    elif ensemble_idx == 2:
        thermostat = args.thermostat

    elif ensemble_idx == 3:
        if thermostat != 3:
            print(
                "NpT ensemble only available with Langevin thermostat. Continuing with MDALGO=3"
            )
        thermostat = 3

    elif ensemble_idx == 4:
        if thermostat != 3:
            print(
                "NpT ensemble only available with Langevin thermostat. Continuing with MDALGO=3"
            )
        thermostat = 3
        inputs["isif"] = 3
        inputs["langevin_gamma"] = 0

    return {**inputs_dict, **inputs}


def input2poscar(input):
    iloc = input.split("/")
    if iloc[-1] != "POSCAR":
        atoms = ase_read(input, index=-1)
        write_vasp("POSCAR", atoms=atoms)
    else:
        if len(iloc) > 1:
            os.system(f"cp {input} POSCAR")
    return 0


def makesuperposcar(atoms, P):
    os.system("mv POSCAR original_primitive.POSCAR")
    supercell = make_supercell(atoms, P=P)
    sorted_supercell = sort_atoms(atoms=supercell)
    write_vasp("POSCAR", sorted_supercell)

    # Reading the POSCAR so we can get the correct order of atoms for POTCAR and U/MM in INCAR
    with open("POSCAR", "r") as poscar:
        poscar_lines = poscar.readlines()
    return tuple(poscar_lines[5].strip().split())


def makePmatrix(terms):
    M = np.zeros([3, 3])
    for i, t in enumerate(terms):
        M[i][i] = t
    return M


def getKPoints():
    if not Path("KPOINTS").is_file():
        print("POINTS does not exists. Writing new default gamma centered.")
        lines = [
            "1x1x1 Gamma-centered mesh\n",
            " 0 0 0\n",
            "Gamma\n",
            " 1 1 1\n",
            " 0 0 0\n",
        ]
        with open("KPOINTS", "w") as kpoints:
            for line in lines:
                kpoints.write(line)
        return False


def handle_hubbard(sym, luj):
    if luj is None:
        print("Hubbard corrections not set.")
        return None
    labels = ["L", "U", "J"]
    n = 4
    elements = []
    d = {}
    separated = [luj[i : i + n] for i in range(0, len(luj), n)]
    for indiv in separated:
        elements.append(indiv[0])
        d[indiv[0]] = dict(zip(labels, [float(x) for x in indiv[-3:]]))
    for s in sym:
        if s not in elements:
            d[s] = dict(
                zip(labels, [-1, 0, 0])
            )  # Negative 1 for no onsite interation added.
    # pprint(d)
    return d


def make_potcar(ordered, potloc):
    def defs(sym):
        if sym in [
            "Li",
            "K",
            "Ca",
            "Sc",
            "Ti",
            "V",
            "Rb",
            "Sr",
            "Y",
            "Rb",
            "Nb",
            "Mo",
            "Cs",
            "Ba",
            "W",
            "Fr",
            "Ra",
        ]:
            return f"{sym}_sv"
        elif sym in [
            "Na",
            "Cr",
            "Mn",
            "Tc",
            "Ru",
            "Rh",
            "Hf",
            "Ta",
        ]:
            return f"{sym}_pv"
        elif sym in [
            "Ga",
            "Ge",
            "In",
            "Sn",
            "Tl",
            "Pb",
            "Bi",
            "Po",
        ]:
            return f"{sym}_d"
        elif sym in [
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Lu",
        ]:
            return f"{sym}_3"
        else:
            return sym

    cmd = ["cat"]
    for s in ordered:
        cmd.append(f"{potloc.rstrip('/')}/{defs(s)}/POTCAR")
    cmd.append("> POTCAR")
    os.system(" ".join(cmd))


def handle_magmoms(syms, magmoms):
    """
    The magamoms variable should be a list parsed by argparse in the form:
    [...] -mgm Fe 5 Nb 0.6 O 0.6 [...]
    which is then converted to a dictionary:
    d = {
        'Fe': 5.,
        'Nb': 0.6,
        'O': 0.6
        }
    """
    elements = magmoms[::2]
    values = magmoms[1::2]
    d = dict(zip(elements, values))
    init_mgm = []
    for s in syms:
        if s not in elements:
            init_mgm.append(0)
        else:
            init_mgm.append(d[s])
    return d, init_mgm


def process_incar(args, ordered):
    atoms = read_vasp("POSCAR")
    atom_dict = Formula(str(atoms.symbols))

    ### Constructing Hubbard correction input
    if args.hubbard_correctionLUJ is not None:
        L, U, J = [], [], []
        hd = handle_hubbard(ordered, args.hubbard_correctionLUJ)
        print(hd)
        for s in ordered:
            L.append(str(int(hd[s]["L"])))
            U.append(str(hd[s]["U"]))
            J.append(str(hd[s]["J"]))
        LDAU_lines = [
            f"LDAU = .TRUE.\n",
            f"LDAUL = {' '.join(L)}\n",
            f"LDAUU = {' '.join(U)}\n",
            f"LDAUJ = {' '.join(J)}\n",
        ]

    ### Construcintg magnetic moment input
    if args.magmoms is not None:
        print("Using ISPIN = 2")
        mgd, mgl = handle_magmoms(ordered, args.magmoms)
        mg_list = []
        for s in ordered:
            if s not in mgd.keys():
                magmom = 0
            else:
                magmom = mgd[s]
            mg_list.append(f"{int(atom_dict[s])}*{magmom}")
        mg_lines = [f"MAGMOM = {' '.join(mg_list)}\n", "ISPIN = 2\n"]

    with open("INCAR", "r") as incar_r:
        lines = incar_r.readlines()
    new_lines = []
    verboten = ["MAGMOM", "LDAU", "ISPIN", "TEBEG", "POTIM", "NSW", "IBRION"]
    for line in lines:
        if not any(word in line for word in verboten):
            new_lines.append(line)
    new_lines.extend(LDAU_lines)
    new_lines.extend(mg_lines)
    new_lines.append("IBRION = 0\n")
    new_lines.append(f"NSW = {args.steps}\n")
    new_lines.append(f"POTIM = {args.timestep}\n")
    os.system("mv INCAR INCAR.old")
    with open("INCAR", "w") as incar_w:
        incar_w.writelines(new_lines)


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="A straight forward approach to setting up AIMD simulations with VASP.",
        epilog=f"{__file__.split('/')[-1]} Benedict Saunders, 2022",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        default="POSCAR",
        help="The input file, can be anything readable by ASE. Def.: POSCAR",
    )
    parser.add_argument(
        "-T",
        "--temperatures",
        required=True,
        default=None,
        nargs="+",
        help="The temperatures, in Kelvin, for which to perform the AIMD runs at.",
    )
    parser.add_argument(
        "-t",
        "--timestep",
        default=2,
        required=False,
        type=int,
        help="The timestep, in femtoseconds, for the molecular dynamics run.",
    )
    parser.add_argument(
        "-s",
        "--steps",
        default=1000,
        type=int,
        help="Number of steps for the molecular dynamics run.",
    )
    parser.add_argument(
        "-e",
        "--ensemble",
        default="NVT",
        type=str,
        help="MD Ensemble. Choose from: NVT (default), NpT or NVE.",
    )
    parser.add_argument(
        "-a",
        "--thermostat",
        default=2,
        type=int,
        help="Thermostat algorithm. Chose from 1: Andersen, 2: Nose-Hoover (default), 3: Langevin, or 13: Multiple Andersen",
    )
    parser.add_argument(
        "-f",
        "--friction",
        type=float,
        default=None,
        required=False,
        help="Friction coefficient if using Langevin dynamics. Def.: None",
    )
    parser.add_argument(
        "-P",
        "--supercell",
        default="111",
        help="Supercell creation, e.g.: -P 333 for a 3x3x3 cell. Def.: 222",
        required=False,
    )

    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        default=128,
        required=False,
        help="Number of cores to run VASP with. Def.: 128",
    )
    parser.add_argument(
        "-vexe",
        "--vasp_executable",
        type=str,
        default="vasp_std",
        required=False,
        help="VASP executable command. Def.: 'vasp_gam'",
    )
    parser.add_argument(
        "-pp",
        "--vasp_potentials",
        type=str,
        default="~/APPS/vasp.5.4.1/PPs/potpaw_PBE",
        required=False,
        help="Location of potpaw potentials for VASP",
    )

    parser.add_argument(
        "-mpi",
        "--mpi",
        action="store_true",
        help="Overrides the srun command to mpirun, and changes the syntx of the command accordingly.",
    )

    parser.add_argument(
        "-N",
        "--nodes",
        default=1,
        type=int,
        help="Define the number of nodes as required by the srun command. Def.: 1",
    )

    parser.add_argument(
        "-mgm",
        "--magmoms",
        default=None,
        help="Magnetic moments for a colinear calculation. Eg, 'Fe 5.0 Nb 0.6 O 0.6' Def.: None. If a defined element is not present in the POSCAR, no MAGMOM will be set for it.",
        nargs="*",
        required=False,
    )

    parser.add_argument(
        "-LUJ",
        "--hubbard_correctionLUJ",
        default=None,
        nargs="*",
        required=False,
        help="Hubbard corrections. Usage: <element 1> <L1> <U1> <J1> ... <element n> <Ln> <Un> <Jn>. Def.: None",
    )

    parser.add_argument(
        "-x",
        "--submission",
        default="aimd.bash",
        required=False,
        help="The bash file for submission to the sbatch queuing system.",
    )

    parser.add_argument(
        "-dry",
        "--dryrun",
        default=False,
        required=False,
        action="store_true",
        help="Do everything except submit to the queuing system.",
    )

    args = parser.parse_args()

    if args.thermostat == 4:
        args.thermostat = 13
        print(
            "Check your input! We fixed it for you this time, but watch out in future!"
        )

    temps = args.temperatures

    potim = args.timestep
    nsw = args.steps

    md_params = handle_input_combinations(args.thermostat, args.ensemble, {})
    print(md_params)

    _ = input2poscar(args.input)
    terms = list(args.supercell)
    P = makePmatrix(terms=terms)
    required = ["POSCAR", "INCAR", "POTCAR", "KPOINTS", args.submission]

    user_kpoints = getKPoints()
    atoms = read_vasp(f"POSCAR")
    species_ordered = makesuperposcar(atoms=atoms, P=P)
    process_incar(args, species_ordered)
    make_potcar(ordered=species_ordered, potloc=args.vasp_potentials)

    for t in temps:
        with cd(f"T_{t}"):
            req = [f"../{r}" for r in required]
            os.system(f"cp {' '.join(req)} .")
            # Changing the temperatures of each incar file.
            with open("INCAR", "r") as incar_r:
                lines = incar_r.readlines()
            with open("INCAR", "w") as incar_w:
                for line in lines:
                    if not "TEBEG = " in line:
                        incar_w.write(line)
                incar_w.write(f"TEBEG = {t}\n")
            if not args.dryrun:
                os.system(f"sbatch {args.submission}")

    """
    
    What am i doing?
    
    """
