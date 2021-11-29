#!/usr/bin/env python
# encoding: utf-8

import argparse
import subprocess as sp

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        'Script to visualize the input and output of a dimer calculation')

    parser.add_argument('-r', '--rotations', nargs='*', default=['0z'])

    args = parser.parse_args()

    rotations = args.rotations

    for r in rotations:
        cmd = "vasp_visualize_input.py -i ./POSCAR ./CONTCAR -rc 1 2 -r {rotation} -o montage_{suffix}.png".format(rotation=r, suffix=r.strip().replace('-', 'm'))
        print cmd
        sp.check_output(cmd, shell=True)
