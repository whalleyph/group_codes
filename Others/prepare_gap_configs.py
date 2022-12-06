import argparse as ap
from glob import glob
import os
import sys

def write_file(name, lines):
    with open(name, "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")
    print(f"Written {len(lines)} line(s) to {name}.")
    return 0

def single_config(atoms_filename, gap_file, energy_sigma, forces_sigma, virial_sigma):
    lines = [
        f"atoms_filename={atoms_filename}",
        f"gp_file={gap_file}",
        "gap={",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 zeta=1 covariance_type=dot_product sparse_method=uniform :",
        "soap cutoff=5.0 atom_sigma=0.5 n_max=10 l_max=6 delta=1.0 covariance_type=dot_product zeta=4 n_sparse=2000 sparse_method=cur_points:",
        f"}} default_sigma={{{energy_sigma} {forces_sigma} {virial_sigma} 0.0}} e0={{Li:-3.07586323:P:-5.0837987:S:-4.6637181}}",
        "force_parameter_name=forces",
    ]
    write_file("single_config", lines)
    return 0

def sparsex_config(atoms_filename, gap_file, energy_sigma, forces_sigma, virial_sigma):
    globbed = glob("*.xml.sparseX.GAP_*.input")
    globbed.sort()
    sparsex_inputs = [g.split("/")[-1] for g in globbed]
    lines = [
        f"atoms_filename={atoms_filename}",
        f"gp_file={gap_file}",
        "gap={",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F",
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[0]} Z1=3 Z2=3:",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F",
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[1]} Z1=3 Z2=15:",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F"
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[2]} Z1=3 Z2=16:",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F",
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[3]} Z1=15 Z2=15:",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F",
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[4]} Z1=15 Z2=16:",
        "distance_2b cutoff=6.0 cutoff_transition_width=1.0 delta=0.1 n_sparse=100 n_exponents=14 add_species=F",
        f"zeta=1 covariance_type=dot_product sparse_method=file sparse_file={sparsex_inputs[5]} Z1=16 Z2=16:",

        "soap cutoff=5.0 atom_sigma=0.5 n_max=10 l_max=6 delta=1.0 covariance_type=dot_product add_species=F",
        f"zeta=4 n_sparse=2000 sparse_method=file sparse_file={sparsex_inputs[6]} n_species=3 Z=3 species_Z={{3 15 16}}:",
        "soap cutoff=5.0 atom_sigma=0.5 n_max=10 l_max=6 delta=1.0 covariance_type=dot_product add_species=F",
        f"zeta=4 n_sparse=2000 sparse_method=file sparse_file={sparsex_inputs[7]} n_species=3 Z=15 species_Z={{3 15 16}}:",
        "soap cutoff=5.0 atom_sigma=0.5 n_max=10 l_max=6 delta=1.0 covariance_type=dot_product add_species=F",
        f"zeta=4 n_sparse=2000 sparse_method=file sparse_file={sparsex_inputs[8]} n_species=3 Z=16 species_Z={{3 15 16}}:",
        f"}} default_sigma={{{energy_sigma} {forces_sigma} {virial_sigma} 0.0}} e0={{Li:-3.07586323:P:-5.0837987:S:-4.6637181}}",
    ]
    write_file("config", lines)
    return 0

energy_sigma = 0.001
forces_sigma = 0.005
virial_sigma = 0.01
atoms_filename = sys.argv[2]
gap_file = sys.argv[3]  

if sys.argv[1] == "single":
    single_config(atoms_filename, gap_file, energy_sigma, forces_sigma, virial_sigma)

elif sys.argv[1] == "sparsex":
    sparsex_config(atoms_filename, gap_file, energy_sigma, forces_sigma, virial_sigma)

else:
    print("N O P E")
    exit()

