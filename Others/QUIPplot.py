# general imports
import os
from random import triangular
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as cp
import argparse as ap
from pprint import pprint

# ase imports
import ase.io

# from ase import Atoms, Atom
# from ase import units
# from ase.build import molecule

# for MD
# from ase.md.langevin import Langevin
from pathlib import Path

# from ase.io.trajectory import Trajectory


def rms_dict(x_ref, x_pred):
    """Takes two datasets of the same shape and returns a dictionary containing RMS error data"""

    x_ref = np.array(x_ref)
    x_pred = np.array(x_pred)

    if np.shape(x_pred) != np.shape(x_ref):
        raise ValueError("WARNING: not matching shapes in rms")

    error_2 = (x_ref - x_pred) ** 2

    average = np.sqrt(np.average(error_2))
    std_ = np.sqrt(np.var(error_2))

    return {"rmse": average, "std": std_}


def energy_plot(in_file, out_file, ax, title="Plot of energy"):
    """Plots the distribution of energy per atom on the output vs the input"""
    # read files
    in_atoms = ase.io.read(in_file, ":")
    out_atoms = ase.io.read(out_file, ":")
    # list energies
    ener_in = [
        at.get_potential_energy() / len(at.get_chemical_symbols()) for at in in_atoms
    ]
    ener_out = [
        at.get_potential_energy() / len(at.get_chemical_symbols()) for at in out_atoms
    ]
    # scatter plot of the data
    ax.scatter(ener_in, ener_out)
    # get the appropriate limits for the plot
    for_limits = np.array(ener_in + ener_out)
    elim = (for_limits.min() - 0.05, for_limits.max() + 0.05)
    # print(elim)
    # ax.set_xlim(elim)
    # ax.set_ylim(elim)
    # add line of slope 1 for refrence
    ax.plot(elim, elim, c="k")
    # set labels
    ax.set_ylabel("energy by GAP / eV")
    ax.set_xlabel("energy by VASP / eV")
    # set title
    ax.set_title(title)
    # add text about RMSE
    _rms = rms_dict(ener_in, ener_out)
    rmse_text = (
        "RMSE:\n"
        + str(np.round(_rms["rmse"], 3))
        + " +- "
        + str(np.round(_rms["std"], 3))
        + "eV/atom"
    )
    ax.text(
        0.9,
        0.1,
        rmse_text,
        transform=ax.transAxes,
        fontsize="large",
        horizontalalignment="right",
        verticalalignment="bottom",
    )


def force_plot(in_file, out_file, ax, symbol, title="Plot of force"):
    """Plots the distribution of firce components per atom on the output vs the input
    only plots for the given atom type(s)"""

    in_atoms = ase.io.read(in_file, ":")
    out_atoms = ase.io.read(out_file, ":")

    # extract data for only one species
    in_force, out_force = [], []
    for at_in, at_out in zip(in_atoms, out_atoms):
        # get the symbols
        sym_all = at_in.get_chemical_symbols()
        # add force for each atom
        for j, sym in enumerate(sym_all):
            if sym in symbol:
                in_force.append(at_in.get_forces()[j])
                # out_force.append(at_out.get_forces()[j]) \
                out_force.append(
                    at_out.arrays["force"][j]
                )  # because QUIP and ASE use different names
    # convert to np arrays, much easier to work with
    # in_force = np.array(in_force)
    # out_force = np.array(out_force)
    # scatter plot of the data
    ax.scatter(in_force, out_force)
    # get the appropriate limits for the plot
    for_limits = np.array(in_force + out_force)
    flim = (for_limits.min() - 1, for_limits.max() + 1)
    ax.set_xlim(flim)
    ax.set_ylim(flim)
    # add line of
    ax.plot(flim, flim, c="k")
    # set labels
    ax.set_ylabel("force by GAP / (eV/Å)")
    ax.set_xlabel("force by VASP / (eV/Å)")
    # set title
    ax.set_title(title)
    # add text about RMSE
    _rms = rms_dict(in_force, out_force)
    rmse_text = (
        "RMSE:\n"
        + str(np.round(_rms["rmse"], 10))
        + " +- "
        + str(np.round(_rms["std"], 10))
        + "eV/Å"
    )
    ax.text(
        0.9,
        0.1,
        rmse_text,
        transform=ax.transAxes,
        fontsize="large",
        horizontalalignment="right",
        verticalalignment="bottom",
    )


parser = ap.ArgumentParser()
parser.add_argument(
    "-t",
    "--train",
    help="Location of training file, def.: train.xyz",
    default="train.xyz",
)
parser.add_argument(
    "-v",
    "--validate",
    help="Location of training file, def.: validate.xyz",
    default="validate.xyz",
)
parser.add_argument(
    "-qt",
    "--quiptrain",
    help="Location of training file, def.: quip_train.xyz",
    default="quip_train.xyz",
)
parser.add_argument(
    "-qv",
    "--quipvalidate",
    help="Location of training file, def.: quip_validate.xyz",
    default="quip_validate.xyz",
)
parser.add_argument(
    "-m",
    "--model",
    help="Location of GAP XML, def.: GAP.xml",
    default="GAP.xml",
)
parser.add_argument(
    "--save",
    help="Save figure",
    action="store_true",
)

parser.add_argument(
    "--title",
    default="Untitled",
)

parser.add_argument("--fname", help="Image file name", default="fit_all.png")

args = parser.parse_args()

fig, ax_list = plt.subplots(nrows=2, ncols=2, gridspec_kw={"hspace": 0.3})
fig.set_size_inches(20, 15)
ax_list = ax_list.flat[:]

quip_train = args.quiptrain
quip_validate = args.quipvalidate
train = args.train
validate = args.validate
model = args.model

print(f"Plotting {args.title}")

energy_plot(train, quip_train, ax_list[0], "Energy on training data")
energy_plot(validate, quip_validate, ax_list[1], "Energy on validation data")

force_plot(train, quip_train, ax_list[2], "LiPS", "Force on training data")
force_plot(validate, quip_validate, ax_list[3], "LiPS", "Force on validation data")

# if you wanted to have the same limits on the force plots
# for ax in ax_list[2:]:
#    flim = (-20, 20)
#    ax.set_xlim(flim)
#    ax.set_ylim(flim)


fig.tight_layout()
fig.suptitle(args.title, fontsize=36)
if args.save:
    fig.savefig(args.fname)
    print(f"  Saved to '{args.fname}'")
else:
    fig.show()
