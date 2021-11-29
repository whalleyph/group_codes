#! /usr/bin/env python

import argparse
from vasp_poscar_normal_to_selective import normal_to_selective

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="""script to fix an atom in the xy-plane""")
    
    parser.add_argument('-f', '--file', default='POSCAR')
    parser.add_argument('-n', '--atoms-to-be-fixed', required=True)
    parser.add_argument('-o', '--out-file', default='POSCAR.fixed')

    args = parser.parse_args()

    poscar_lines = open(args.file).readlines()

    # determine if the poscar already has selective dynamics
    if poscar_lines[7].strip()[0].lower() == 's':
        selective = True
    else:
        selective = False

    # if not selective convert to selective dynamics file
    if not selective:
        poscar_lines = normal_to_selective(args.file).split('\n')

    # get the lines where all the positions are present
    positions = [l.strip() for l in poscar_lines[9:]]

    # get the atoms to be fixed in the xy-plane
    atoms_to_be_fixed = eval(args.atoms_to_be_fixed)

    # check if it is just one atom
    try:
        atoms_to_be_fixed = [int(atoms_to_be_fixed)]
    except:
        pass

    for i in atoms_to_be_fixed:
        pos = positions[i-1].split()[:-3]
        pos += 'F F T'.split()
        positions[i-1] = '  '.join(pos)

    poscar_lines_before_positions = [l.strip() for l in poscar_lines[:9]]
    poscar_fixed = '\n'.join(poscar_lines_before_positions + positions)

    open(args.out_file, 'w').write(poscar_fixed)
