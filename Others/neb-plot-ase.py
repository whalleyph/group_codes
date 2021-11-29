#!/usr/bin/env python
import matplotlib.pyplot as plt
from ase.neb import NEBTools
from ase.io import read,write
import argparse
from os import system
from os import listdir,system,access,R_OK

parser = argparse.ArgumentParser()
parser.add_argument("-n","--NoImages", type=int,required=0, default=0,
                    help="No of images in the NEB chain (inclding inital and final structures)")
parser.add_argument("-s","--show", action='store_true', default=False,
                    help="Show resulting plot. It is automaticall saved in 'diffusion-barrier.png'")

args = parser.parse_args()

fl = listdir(".")
fl.sort()

if args.NoImages==0:
    
    #for fn in fl:
    #    if fn[0:3]=="IMG" :
    #        ind = int(fn[3:])
    
    args.NoImages=len([fn for fn in fl if fn[0:3]=='IMG'])
    print ("%d images detected..."%args.NoImages)

images = read('neb.traj@-%d:'%args.NoImages)

write('neb.xyz',images,format='xyz')

nebtools = NEBTools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()

print ('Ef,dE',Ef,dE)

# Get the barrier without any interpolation between highest images.
#Ef, dE = nebtools.get_barrier(fit=0)

# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()

# Create a figure like that coming from ASE-GUI.
fig = nebtools.plot_band()
#fig.savefig('diffusion-barrier.png')

# Create a figure with custom parameters.
fig = plt.figure(figsize=(5.5, 4.0))
ax = fig.add_axes((0.15, 0.15, 0.8, 0.75))
nebtools.plot_band(ax)
fig.savefig('diffusion-barrier.png')

system('$HOME/bin/neb-movie-ase.py')

#if args.show:fig.show() #doesn't work
if args.show: system('display diffusion-barrier.png&')


