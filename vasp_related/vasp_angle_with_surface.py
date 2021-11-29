#! /usr/bin/env python

import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-a1', '--a1', type=float, nargs=3, required=True)
parser.add_argument('-a2', '--a2', type=float, nargs=3, required=True)

args = parser.parse_args()
a1 = np.array(args.a1)
a2 = np.array(args.a2)

a1a2 = a2 - a1

# angle between surface normal and a1a2
angle_radians = np.arccos(np.dot(a1a2, np.array([0, 0, 1]))/np.linalg.norm(a1a2))
rad2deg = 180./np.pi
angle_degrees = angle_radians*rad2deg
print 'Angle with surface normal: %.2f' % (angle_degrees)
print 'Angle with surface: %.2f' % (90 - angle_degrees)
