#!/usr/bin/env python
# encoding: utf-8

import numpy as np

class acf():

    """Class to read ACF file produced by Bader Charge Analysis of the VTST code"""

    def __init__(self, fname='ACF.dat'):

        self._fname = fname

        self._data = np.genfromtxt(self._fname, skip_header=2, skip_footer=4)


    @property
    def positions(self):

        return self._data[:,1:4]


    @property
    def charges(self):

        return self._data[:,4]


    @property
    def min_distances(self):

        return self._data[:,5]


    @property
    def atomic_volumes(self):

        return self._data[:,6]
