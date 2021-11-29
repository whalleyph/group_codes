#! /usr/bin/env python

import numpy as np
from cStringIO import StringIO as sio
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from os import system
import argparse

"""Plots the total DoS from the DOSCAR file from VASP"""

parser = argparse.ArgumentParser(description='script to create total DOSplot from a DOSCAR file')

parser.add_argument('-i', '--input_file', default='DOSCAR',
                    help='input file in DOSCAR format')

parser.add_argument('-f', '--fermi', action='store_true', default=False,
                        help='If Fermi level should be shifted to zero.')

parser.add_argument('-x1', '--x1', default=0.0, type=float,
            help='Start point of the x axis')
parser.add_argument('-x2', '--x2', default=0.0, type=float,
            help='End point of the x axis')
parser.add_argument('-y1', '--y1', default=0.0, type=float,
            help='Start point of the y axis')
parser.add_argument('-y2', '--y2', default=0.0, type=float,
            help='End point of the y axis')

parser.add_argument('-s', '--shift', default=0.0, type=float,
            help='Amount of shift in the x direction')

args = parser.parse_args()
fermi=args.fermi

#f = open('DOSCAR')
f = open(args.input_file)
l = f.readlines()

natoms = int(l[0].strip().split()[1])
hold = l[5].strip().split()
ndat, efermi = int(hold[2]), float(hold[3])
#print ndat, efermi

dat = sio(''.join(l[6:6+ndat]).strip())
energy, dos, intdos = np.loadtxt(dat, unpack=True)

#if args.shift != 0.0: energy += args.shift
if fermi: shifted_energy = energy - efermi
elif args.shift != 0.0: shifted_energy = energy + args.shift
else: shifted_energy = energy
dos_per_atom = dos/natoms

#for ienergy, idos, iintdos in zip(energy, dos, intdos):
#    print energy, dos, intdos

fig = plt.figure()
plt.minorticks_on()
ax = fig.add_subplot(111)

ax.plot(shifted_energy, dos_per_atom)


#ax2 = fig.add_subplot(221)
#ax2.plot(shifted_energy, dos_per_atom)

#if fermi and args.shift==0.0 : xlabel = 'E - Efermi [%.2f eV]' % (efermi)
if fermi : xlabel = 'E - Ef [%.2f eV]' % (efermi)
#elif fermi and args.shift!=0.0: xlabel = 'E - Efermi + shift [%.2f and %.2f eV]' % (efermi,args.shift)
elif not fermi and args.shift!=0.0 :xlabel = 'E + shift [%.2f eV]' % (args.shift)
else: xlabel = 'E [eV]'

ax.set_xlabel(xlabel)

x1=args.x1;x2=args.x2;y1=args.y1;y2=args.y2
xx1,xx2,yy1,yy2 = plt.axis()

ax.set_ylabel('DoS/atom')
#if (x1 == 0.0 and x2==0.0 and y1==0.0 and y2==0.0):
#    x1,x2,y1,y2 = plt.axis()
if x1==0.0: x1=xx1
if x2==0.0: x2=xx2
if y1==0.0: y1=yy1
if y2==0.0: y2=yy2


plt.axis((x1,x2,y1,y2))

#plt.autoscale(axis="y")
#ymin,ymax=plt.ylim()
#print ymin,ymax
#plt.axis((x1,x2,0,ymax))

ax = fig.add_subplot(111)
ax.plot(shifted_energy, dos_per_atom)

plt.show()
plt.savefig('dos.png')

system("display dos.png &")

