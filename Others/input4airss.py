import argparse as ap
import itertools

from chemparse import parse_formula as Formula
from mendeleev import element as elemental_data


def vprint(s, v=True):
    if v:
        print(s)


class pyairss_input:
    def __init__(self, verbose=True):
        self.verbose = verbose

    def set_phases(self, phases):
        if self.verbose:
            vprint(
                "Creation of randomised pyAIRSS input for phases along the tieline:",
                self.verbose,
            )
            vprint(f"{phases[0]} <----------> {phases[1]}\n", self.verbose)

        stoichiometries = []
        combined = {}
        self.phases = phases
        for phase in phases:
            formula = Formula(phase)
            combined.update(formula)
            stoichiometries.append(formula)

        self.elements = set(combined.keys())
        self.stoichs = stoichiometries

    def get_radii(self, scaling=1.00):

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

        vprint(f"Getting covalent radii from Mendeleev", self.verbose)
        for element in self.elements:
            vprint(f" > {element}", self.verbose)
            radii.update(
                {element: float(elemental_data(element).covalent_radius) * scaling}
            )

        self.radii = radii
        return radii

    def make_cell_file(self, name):
        name = f"{name}.cell"
        vprint(f"\nWriting input file {name}", self.verbose)
        self.f = open(name, "w")
        self.f.write("#! define s1 space group number=10 to 230\n\n")

    def write_seperations(self):
        self.combinations = list(
            itertools.combinations_with_replacement(self.elements, 2)
        )
        for combination in self.combinations:
            id = "".join(combination)
            sep = "-".join(combination)
            val = (self.radii[combination[0]] + self.radii[combination[1]]) / 100
            self.f.write(f"#! define {id} minimum separation {sep} {val:.3}\n")
        self.f.write(f"\n")
        return None

    def determine_new_stoichiometries(self, xs, ys, fu=1):
        self.f.write(f"#! define x int {xs[0]} to {xs[1]}\n")
        self.f.write(f"#! define y int {ys[0]} to {ys[1]}\n")
        self.f.write(f"#! define u int {fu}\n\n")

        for stoich in self.stoichs:
            for element in self.elements:
                if element not in stoich.keys():
                    stoich.update({element: 0})

        # Interpolation of phases
        for mol in self.stoichs[0].keys():
            self.f.write(f"#! define n{mol} int x*u*{int(self.stoichs[0][mol])}\n")
        self.f.write(f"\n")
        for mol in self.stoichs[1].keys():
            self.f.write(
                f"#! define n{mol} int n{mol}+(y*u*{int(self.stoichs[1][mol])})\n"
            )
        self.f.write("\n")

    def constrain_and_randomise(self, vol_scale_range=(5, 15)):
        sum_statment = "+n".join(self.elements)
        self.f.write(f"#! define nTotal int n{sum_statment}\n")
        self.f.write(
            f"#! define cv volume constraint nTotal*{vol_scale_range[0]} to nTotal*{vol_scale_range[1]}\n"
        )
        self.f.write(f"#! randomise basis vectors subject to s1, cv\n\n")
        pairs_list = ["".join(c) for c in self.combinations]
        pairs = ",".join(pairs_list)
        for element in self.elements:
            self.f.write(
                f"#! insert and randomise n{element} {element} atoms into cell subject to s1, {pairs}\n"
            )
        self.f.write("\n")

    def finalise(self, save_as_opt=True):
        if save_as_opt:
            self.f.write("\n#! save as optimised")
        self.f.close()
        vprint("Finished.", self.verbose)


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument(
        "-p",
        "--phases",
        help="Stoichiometry of the two phases:\n -p <phase 1 formula> <phase 2 formula>",
        nargs=2,
        type=str,
    )
    parser.add_argument(
        "-s",
        "--seed",
        dest="seed",
        help="Seed name for input file for PyAIRSS. Default: <phase1>-<phase2>",
        default="default-seed",
        type=str,
        required=False,
    )

    args = parser.parse_args()

    if args.seed != "default-seed":
        seed = args.seed
    else:
        seed = f"{args.phases[0]}-{args.phases[1]}"

    new_phases = pyairss_input()
    new_phases.set_phases(args.phases)
    new_phases.get_radii(scaling=0.85)
    new_phases.make_cell_file(seed)
    new_phases.write_seperations()
    new_phases.determine_new_stoichiometries(xs=(1, 15), ys=(1, 15), fu=1)
    new_phases.constrain_and_randomise(vol_scale_range=(5.0, 15.0))
    new_phases.finalise(save_as_opt=True)
