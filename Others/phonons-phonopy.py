#!/usr/bin/env python3
import phonopy
#from phonopy import Phonopy
from phonopy.units import VaspToTHz

#unitcell, _ = read_crystal_structure("POSCAR-unitcell", interface_mode='vasp')

#Procedure:
#1. Create a supercell from the input geom using 
# phonopy -d --dim="2 2 2" -c POSCAR-unitcell ;  mv SPOSCAR POSCAR
#2. Run the VASP DFPT or Finite-Disp methods (i.e. IBRION=8 or IBRION=6 and NSW=1) seperately and read the vaprun.xml to create FORCE_CONSTANTS file using
# phonopy --fc vasprun.xml


phonon = phonopy.load(
                      #supercell_matrix=[2, 2, 2],
                      #supercell_filename="SPOSCAR",
                      primitive_matrix='auto',
                      unitcell_filename="POSCAR",
                      factor=VaspToTHz,
                      #force_constants_filename="force_constants.hdf5")
                      force_constants_filename="FORCE_CONSTANTS")

"""
path = [[[0, 0, 0], [0.5, 0, 0.5], [0.625, 0.25, 0.625]],
        [[0.375, 0.375, 0.75], [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75]]]
labels = ["$\\Gamma$", "X", "U", "K", "$\\Gamma$", "L", "W"]
qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
phonon.plot_band_structure().show()
"""

phonon.auto_band_structure()

#These commands also work:
#phonon.auto_band_structure(plot=True).show() #also works directly
#phonon.plot_band_structure().show() #this also works
#phonon.auto_total_dos(plot=True).show()
#phonon.auto_projected_dos(plot=True).show()

# To plot DOS next to band structure
phonon.run_mesh([20, 20, 20])
phonon.run_total_dos()
phonon.plot_band_structure_and_dos().show()


# To plot PDOS next to band structure
#phonon.run_mesh([20, 20, 20], with_eigenvectors=True, is_mesh_symmetry=False)
#phonon.run_projected_dos()
#phonon.plot_band_structure_and_dos(pdos_indices=[[0], [1]]).show()


#MESH SAMPLING Calcualtion
#the mesh sampling calclation has to be done in the THz unit, which we do with factor=VaspToTHz
phonon.run_mesh([20, 20, 20])
phonon.run_thermal_properties(t_step=10,
                              t_max=1000,
                              t_min=0)
tp_dict = phonon.get_thermal_properties_dict()
temperatures = tp_dict['temperatures']
free_energy = tp_dict['free_energy']
entropy = tp_dict['entropy']
heat_capacity = tp_dict['heat_capacity']

if 0:
  print ('temperature, free_energy, entropy, heat_capacity:')
  for t, F, S, cv in zip(temperatures, free_energy, entropy, heat_capacity):
    print(("%12.3f " + "%15.7f" * 3) % ( t, F, S, cv ))

phonon.plot_thermal_properties().show()



