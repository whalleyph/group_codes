#!/usr/bin/env python
# encoding: utf-8

import ase.io
import numpy as np
from potcar import potcar
from poscar import poscar
from vtstbader import acf
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Script to get partial charges from a VTST Bader Analysis')

    parser.add_argument('-i1', '--fposcar', default='CONTCAR')
    parser.add_argument('-i2', '--fpotcar', default='POTCAR')
    parser.add_argument('-i3', '--facfdat', default='ACF.dat')
    parser.add_argument('-o', '--output-file', default='partial_charges_bader.xyz')

    args = parser.parse_args()
    fposcar, fpotcar, facfdat, output_file = (
            args.fposcar, args.fpotcar, args.facfdat, args.output_file)

    nions = poscar().read(fposcar).get_nions()
    nelectrons = potcar().read(fpotcar).get_nelectrons()

    valence_electrons = []
    for i in range(len(nions)):
        valence_electrons.extend(nions[i]*[nelectrons[i]])

    bader_charges = acf(facfdat).charges
    atoms = ase.io.read(fposcar)

    partial_charges_bader = [e-q for 
            (e,q) in zip(valence_electrons, bader_charges)]

    #print partial_charges_bader

    atoms.arrays['partial_charges_bader'] = np.array(partial_charges_bader)

    ase.io.write(output_file, atoms, format='extxyz',
            columns=['symbols', 'positions', 'partial_charges_bader'])

    #for atom, e, q in zip(atom_symbols, valence_electrons, bader_charges):
    #    print '%2s %10.2f' % (atom, e-q)

