#! /usr/bin/env python

# valid only for vasp5 poscar files

import argparse

def normal_to_selective(f):
    lines = open(f).readlines()

    # get the total number of atoms in the file
    try:
        natoms = int(lines[6].strip())
    except:
        natoms = sum([int(i) for i in lines[6].strip().split()])

    top = [l.strip() for l in lines[:7]]

    coordinates_type = lines[7].strip()
    positions = [l.strip() for l in lines[8:8+natoms]]

    for i, l in enumerate(positions):
        positions[i] += '  T  T  T'

    selective_poscar = '\n'.join(top + ['Selective Dynamics'] + [coordinates_type] + positions)

    return selective_poscar.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="""script to convert a normal poscar to a
            selective dynamics poscar""")
    
    parser.add_argument('-f', '--file', default='POSCAR')
    parser.add_argument('-o', '--out-file', default='POSCAR.selective')
    
    args = parser.parse_args()

    selective_poscar = normal_to_selective(args.file)

    open(args.out_file, 'w').write(selective_poscar)
