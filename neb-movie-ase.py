#!/usr/bin/env python

#Also check "/home/bk393/APPS/PYTHON-2.7.13/lib/python2.7/site-packages/ase/io/pov.py"
#Default atom colours for ASE can be foudn in //home/bk393/APPS/PYTHON-2.7.13/lib/python2.7/site-packages/ase/data/colors.py
import numpy as np

from ase import Atoms
from ase.io import write,read
from ase.io import pov
from os import system


images= read('neb.xyz',index=':')
#images = read('neb.traj@-%d:'%args.NoImages)

rotation = '-90x, -90y, 0z' # found using ASE-GUI menu 'view -> rotate'

#Make colors for each atom
#from ase.utils import hsv
#colors = hsv(atoms[0].positions[:, 0])
colors=None


# Textures
#tex = ['jmol',] * 288 + ['glass',] * 288+ ['ase3',] * 288 + ['vmd',] * 288
tex = ['jmol',] * 288 +  ['ase3',] * 288 + ['vmd',] * 288


# keywords
kwargs = { # Keywords that exist for eps, png, and pov
'rotation': rotation,
'colors': colors,
#'radii': None,
'radii'         : .85, # float, or a list with one float per atom
'show_unit_cell': 2,   # 0, 1, or 2 to not show, show, and show all of cell
}

extra_kwargs = { # For povray files only
'display'      : False, # Display while rendering
'pause'        : False, # Pause when done rendering (only if display)
'transparent'  : False, # Transparent background
#'canvas_width' : 1024,  # Width of canvas in pixels #Def: None
'canvas_width' : 640,  # Width of canvas in pixels #Def: None
'canvas_height': None,  # Height of canvas in pixels #Def: None
'camera_dist'  : 50.,   # Distance from camera to front atom
'image_plane'  : None,  # Distance from front atom to image plane
                        # (focal depth for perspective)
'camera_type'  : 'perspective', # perspective, orthographic, ultra_wide_angle
#'camera_type': 'orthographic',  # perspective, ultra_wide_angle
'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
'area_light'   : [(2., 3., 40.) ,# location
                  'White',       # color
                  .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
'background'   : 'White',        # color
'textures'     : tex, # Length of atoms list of texture names #Def: None
'celllinewidth': 0.05, # Radius of the cylinders representing the cell
'transmittances': None,  # transmittance of the atoms

# use with care - in particular adjust the camera_distance to be closer
#'depth_cueing': False,  # fog a.k.a. depth cueing
#'cue_density': 5e-3,  # fog a.k.a. depth cueing
'bondlinewidth': 0.10,  # radius of the cylinders representing bonds
#'bondatoms': pov.get_bondpairs(atoms, radius=1.1) ,  # [[atom1, atom2], ... ] pairs of bonding atoms  #Def: []
'exportconstraints': False,  # honour FixAtoms and mark relevant atoms?
}

kwargs.update(extra_kwargs)


system('rm neb.gif')

# Make the raytraced image
for i,frame in enumerate(images):
  print ('Frame: %d'%i)
  #frame.center() #;frame.wrap()

  #To draw bonds between atom pairs based on (covalent radii*radius).
  if 0:  bp=[ [p[0],p[1]] for p in pov.get_bondpairs(frame, radius=1.1)] #radius =scaling factor for the covalent radii,Def: 1.1
  else: bp=[]
  write('neb-movie-%d.pov'%i, frame, run_povray=True, bondatoms=bp,**kwargs)

system('convert -delay 40 -loop 100 `ls -v neb-movie-*.png` neb.gif')
system('rm neb-movie-*.*')

#povray  +I neb-movie.1.pov +Oneb-movie.1.png +W1024 +H768 +V -D +FN +Q9 +P +UD +UL +UV +A +AM2 +UA

