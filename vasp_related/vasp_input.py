#! /usr/bin/env python

import os
from shutil import copy
import argparse

template_dir = td = '/home/cande/mytools/templates/vasp'


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='')

    parser.add_argument('-t', '--type', required=True)
    parser.add_argument('-o', '--output_directory', default='.')

    args = parser.parse_args()

    calculation_type = args.type
    dst_dir = args.output_directory
    dst_dir = os.path.abspath(dst_dir)

    get_src_files = lambda f: os.path.join(td, f)
    get_dst_files = lambda f: os.path.join(dst_dir, f.split('.')[0])

    if calculation_type == 'molecule':
        src_files = ['INCAR.molecule', 'KPOINTS.gamma']
        src_files_full_path = map(get_src_files, src_files)
        dst_files = map(get_dst_files, src_files)

        map(copy, src_files_full_path, dst_files)
