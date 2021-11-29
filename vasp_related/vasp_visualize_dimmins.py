#! /usr/bin/env python

import ase
import ase.io
from glob import glob
import os, shutil
import subprocess as sp

if os.path.exists('min1'):
    dirs = ['min1', '../', 'min2',]
elif os.path.exists('rct'):
    dirs = ['rct', '../', 'prd']

ndirs = len(dirs)
img_dir = 'dim_imgs'

# create all the possible images
if not os.path.exists(img_dir): os.makedirs(img_dir)
for i in xrange(ndirs):
    d = dirs[i]
    fposcar = os.path.join(d, 'POSCAR')
    fcontcar = os.path.join(d, 'CONTCAR')
    for f in [fposcar, fcontcar]:
        if os.path.exists(f):
            fimg = '.'.join(['%2.2i' % (i), f.split(os.path.sep)[-1], 'png'])
            atoms = ase.io.read(f, format='vasp')
            ase.io.write(os.path.join(img_dir, fimg), atoms, format='png')

# make the montage
cmd = ['montage',
    os.path.join(img_dir, '??.POSCAR.png'),
    os.path.join(img_dir, '??.CONTCAR.png'),
    '-geometry',
    '200x200>',
    '-tile',
    'x'.join(['%i' % (ndirs), '2']),
    '-frame',
    '5',
    'dim_montage.png'
    ]

sp.check_call(cmd)
