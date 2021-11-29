#! /usr/bin/env python

import ase
import ase.io
import argparse
import os
import subprocess as sp
from glob import glob
from shutil import rmtree
import sys

def get_fnames_with_labels(fnames, png_files, labels=False, file_tail=False):

    # use file filenames as labels
    if labels and file_tail:
        labels = ['-label %s %s' % (s) for s in zip(fnames, png_files)]

    # remove 'POSCAR' or 'CONTCAR' at the end of the file names
    elif labels and not file_tail:
        labels = ['-label %s %s' % ('_'.join(f.split(os.path.sep)[:-1]), p) 
                                            for f, p in zip(fnames, png_files)]

    # empty labels
    else:
        labels = ['%s' % (s) for s in png_files]

    return ' '.join(labels)

parser = argparse.ArgumentParser(
        description='')

parser.add_argument('-i', '--input-files', nargs='*')
parser.add_argument('-t', '--tmp-directory', default='img')
parser.add_argument('-rc', '--rows-columns', nargs=2, type=int)
parser.add_argument('-r', '--rotation', default='')
parser.add_argument('-o', '--output-file', default='montage.png')
parser.add_argument('-bbox', '--bbox', nargs=4, type=float, default=None)
parser.add_argument('-l', '--labels', action='store_true', default=False)
parser.add_argument('-nft', '--no-file-tail', action='store_false', default=True)
parser.add_argument('-s', '--side', type=int, default=200)

args = parser.parse_args()
rotation = args.rotation.replace('m', '-')

files = args.input_files
tmp_dir = args.tmp_directory

if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

png_files = []
for f in files:
    atoms = ase.io.read(f, format='vasp')
    file_name = f.split(os.path.sep)[-1]
    original_directory = f.split(os.path.sep)[-2]
    png_file = os.path.join(tmp_dir, '%s.%s.png' % (file_name, original_directory))
    ase.io.write(png_file, atoms, format='png', rotation=rotation, bbox=args.bbox)
    png_files.append(png_file)


if args.rows_columns:
    r, c = args.rows_columns
else:
    #c, r = len(png_files), 1
    c = int(len(png_files)**0.5)+1
    r, c = c, c

s = args.side

fnames_with_labels = get_fnames_with_labels(files, png_files, args.labels)
print fnames_with_labels

cmd = ['montage',
    fnames_with_labels,
    '-tile',
    '%ix%i' % (c, r),
    '-frame',
    '5',
    '-geometry',
    '"%ix%i>"' % (s, s),
    args.output_file
    ]


print os.getcwd()
print ' '.join(cmd)

sp.check_call(' '.join(cmd), shell=True)

# remove temporary directory
rmtree(tmp_dir)
