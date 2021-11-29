#!/usr/bin/env python
# encoding: utf-8

import ase
import ase.io
from incar import incar
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Script to set the DIPOL tag in INCAR')

    parser.add_argument('-i', '--input-file', default='INCAR',
                    help = 'input file')
    parser.add_argument('-p', '--poscar-file', default='POSCAR',
                    help = 'poscar file')
    parser.add_argument('-o', '--output-file', default='INCAR_dipol',
                    help = 'output file')
    parser.add_argument('-a', '--set-all-tags', action='store_true', default=False,
                    help = 'adds the tags IDIPOL and LDIPOL as well')
    
    args = parser.parse_args()
    
    input_file = args.input_file
    poscar_file = args.poscar_file
    output_file = args.output_file
    set_all_tags = args.set_all_tags

    atoms = ase.io.read(poscar_file, format='vasp')
    a, b, c = atoms.get_center_of_mass(scaled=True)

    incar = incar().read(input_file)

    incar.settag('DIPOL = %.4f %.4f %.4f' % (a, b, c))

    if set_all_tags:
        incar.settag('LDIPOL = .TRUE.')
        incar.settag('IDPOL = 3')

    incar.write(output_file)
