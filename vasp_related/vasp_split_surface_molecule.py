#! /usr/bin/env python

import os,sys
import argparse
import ase
import ase.io
import numpy as np
from vasp_shift_origin import shift_origin
import pybel
from myatoms import myAtoms

def extract_molecule(atoms, smiles, molecule_atoms, rmax):

    # atoms: surface + molecule

    # smiles: smiles of the molecule to be extracted
    #         extracted molecule is matched against the connectivity
    #         in this smiles

    # filter_atoms: function that is a filter for atoms in the molecule
    # filter_atoms() should return a list of ase.Atom objects

    print smiles
    molecule_to_extract = pybel.readstring('smi', smiles)
    nbonds_required = molecule_to_extract.OBMol.NumBonds()

    print 'Looking for molecule: %s' % (molecule_to_extract.formula)
    print 'Number of bonds in required molecule: %i' % (nbonds_required)

    molecule = ase.Atoms(molecule_atoms, cell=atoms.get_cell(), pbc=atoms.get_pbc())

    molecule_found = False
    for z in np.arange(-0.25,0.26,0.01):
        for y in np.arange(-0.25,0.26,0.01):
            for x in np.arange(-0.25,0.26,0.01):

                # shift origin of the box so that the molecule lies 
                # completely in the box
                molecule = shift_origin(molecule, x, y, z)

                # write the molecule temporarily so that pybel can read it
                #tmpxyz = 'tmp.xyz'
                #ase.io.write(tmpxyz, molecule, format='xyz')
                #extracted_molecule = pybel.readfile('xyz', tmpxyz).next()

                # calculate the number of bonds in the molecule
                #nbonds = extracted_molecule.OBMol.NumBonds()
                nbonds = myAtoms(molecule).get_number_of_bonds(rmax)

                # check if the number of bonds found in the xyz file matches the
                # number of bonds in the smiles string
                if nbonds == nbonds_required:
                    molecule_found = True
                if molecule_found: break
            if molecule_found: break
        if molecule_found: break

    #os.remove(tmpxyz)

    #pos = molecule.get_scaled_positions()
    #pos = pos+np.array([x, y, z])
    #for i in pos:
    #    print i[0], i[1], i[2]

    if molecule_found:
        return molecule
    else:
        print 'Molecule not found'
        sys.exit(-1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--poscar', nargs='?', default='POSCAR',
                help='POSCAR file to use')

    parser.add_argument('-s', '--shift', nargs='*', default=[0.0, 0.0, 0.0],
                help='POSCAR file to use')

    parser.add_argument('-l', '--smiles', required=True,
                help='SMILES string of the molecule to extract')

    parser.add_argument('-rmax', '--rmax', type=float, default=1.6,
                help='maximum bond length to be looked for')

    args = parser.parse_args()

    print 'Using %s' % (args.poscar)

    atoms = ase.io.read(args.poscar, format='vasp')

    molecule_atoms = [i for i in atoms if i.symbol in ['C', 'H']]

    molecule = extract_molecule(atoms, args.smiles, molecule_atoms, args.rmax)

    ase.io.write('molecule.xyz', molecule, format='xyz')
