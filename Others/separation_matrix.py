from mendeleev import element as elemental_data
import argparse as ap
import itertools
import pandas as pd
import time
import numpy as np
import io

UTIME = int(time.time())

def get_radii(elements, scaling=1.00):

    radii = {}
    r_types = [
        "atomic_radius",
        "atomic_radius_rahm",
        "covalent_radius_bragg",
        "covalent_radius_cordero",
        "covalent_radius_pyykko",
        "covalent_radius_pyykko_double",
        "covalent_radius_pyykko_triple",
    ]

    print(f"Getting covalent radii from Mendeleev")
    for element in elements:
        print(f" > {element}")
        radii.update(
            {element: float(elemental_data(element).covalent_radius) * scaling}
        )
    return radii

def get_seperations(elements, radii):
    seps = {}
    combinations = itertools.combinations_with_replacement(elements, 2)
    for combination in combinations:
        val = (radii[combination[0]] + radii[combination[1]]) / 100
        seps[combination] = val
    return seps

def str_to_io(seps):
    keys = seps.keys()
    iol = []
    for key in keys:
        iol.append(f"{key[0]},{key[1]},{seps[key]:.3f}")
    return "\n".join(iol) 

def show_mx(iostr, elements, h):
    print()
    l = len(elements)
    df = pd.read_csv(io.StringIO(iostr), header=None)
    if h:
        tri = df[2].values
        z = np.zeros((l,l), dtype=float)
        z[np.triu_indices(l, 1)] = tri
        z += z.T
    else:
        z = (
            df.set_index([0, 1])
            .reindex(pd.MultiIndex.from_product([elements, elements]))
            .squeeze()
            .unstack(1)
        )
        z = z.where(~z.isna(), z.T)
    print(z)
    print()

parser = ap.ArgumentParser(description="Get an n x n matrix of the elements you list using accepted atomic/covalent/ionic radii.")
parser.add_argument('elems', metavar='Elems', type=str, nargs='+', help="The elements, in order.")
# parser.add_argument('--hide', action="store_true", help="Hide elemental labels.")
# parser.add_argument("-r", "--radius", default="covalent", required=False, help="The type of radius to use")
args = parser.parse_args()

if len(args.elems) < 2:
    print("What")

radii = get_radii(args.elems, 0.85)
seps = get_seperations(args.elems, radii)
iostr = str_to_io(seps)
show_mx(iostr, args.elems, False)

