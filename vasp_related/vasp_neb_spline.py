#! /usr/bin/env python

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# x = reaction coordinate
# y = energy

x, y = np.loadtxt('neb.dat', usecols=(1,2), unpack=True)

f = interpolate.interp1d(x,y,kind='cubic')

xnew = np.linspace(min(x), max(x), 100)
ynew = f(xnew)

plt.plot(x,y,'o',xnew,ynew,'-')

plt.show()
