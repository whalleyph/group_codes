#!/usr/bin/env python

import argparse

def fix_coordinates(fposcar, list_of_atoms):
    lines = open(fposcar).readlines()
    lines = [line.strip() for line in lines]
    top_lines = lines[:8]
    coordinates = lines[8:]

    newcoordinates = []
    for i, ic in enumerate(coordinates):
        if i in list_of_atoms:
            line_with_constraints = ic + ' T T F'
        else:
            line_with_constraints = ic + ' T T T'
        newcoordinates.append(line_with_constraints)

    newposcar = top_lines + newcoordinates
    newposcar = newposcar[:7] + ['Selective Dynamics'] + newposcar[7:]
    newposcar = '\n'.join(newposcar)

    return newposcar


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='fix the z coordinates of the given atoms')
    
    parser.add_argument('-n', '--natoms', nargs='*', type=int)
    parser.add_argument('-r', '--range')
    parser.add_argument('-f', '--fposcar', default='POSCAR')

    args = parser.parse_args()

    if args.natoms:
        natoms = [i-1 for i in args.natoms]
    elif args.range:
        natoms = eval('range('+args.range+')')

    newposcar = fix_coordinates(args.fposcar, natoms)

    f = open(args.fposcar, 'w')
    f.write(newposcar)
