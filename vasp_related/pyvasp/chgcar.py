#! /bin/env python

import os
from subprocess import Popen
from shutil import move
from sys import stderr
import numpy as np
import cStringIO as sio

class chgcar:
    def __init__(self):

        self._file_ = 'CHGCAR'
        self.tags = {}


    def read(self, file=None):
        if file is not None:
            try:
                file.pos
            except AttributeError:
                self._file_ = open(file)
            else:
                self._file_ = file
                self._file_.seek(0)
        else:
            self._file_ = open(self._file_)

        return self

    def get_magmom_from_locprop(self, file):
        lines = open(file).readlines()
        for line in lines:
            if line.find('MAGMOM') > -1:
                magmoms = [float(i) for i in line.split('=')[-1].split()]
        return magmoms

    def _get_rootdir(self):
        sep = os.path.sep
        rootdir = sep.join(self._file_.name.split(sep)[:-1])
        return rootdir


    def get_magmom(self):
        """Use Marcel's readchgx to do a Voronoi analysis and get the
        magnetic moments on the atoms"""

        readchgx = '/home/cande/bin/readchgx'
        rootdir = self._get_rootdir()
        locprop = os.path.join(rootdir, 'localproperties.txt')
        # Check if there already is the file local_properties.txt 
        # newer than chgcar
        if os.path.exists(locprop) and os.stat(locprop)[-2] > \
                                        os.stat(self._file_.name)[-2]:
            magmoms = self.get_magmom_from_locprop(locprop)
        else:
            cmd = ' '.join([readchgx, self._file_.name])
            print cmd
            Popen(cmd, shell=True).wait()
            magmoms = self.get_magmom_from_locprop('localproperties.txt')
            # Move files
            filelist = 'totalcharge.txt spindensity.txt \
                    localproperties.txt'.split()
            for file in filelist:
                move(file, os.path.join(rootdir, file))

        return magmoms

    def get_nions(self):
        snions = self._file_.readlines(1024)[5]
        iions = [int(i) for i in snions.split()]
        nions = sum(iions)
        return nions

    def get_charges(self):
        """Use Bader analysis to get the charges on each of the atoms in the
        supercell
        """
        if os.environ.has_key('BADER_EXEC'):
            badexexec = os.environ['BADER_EXEC']
        else:
            print 'Set up environment variable BADER_EXEC pointing to the \
                    bader analysis executable'
        sep = os.path.sep
        fchgcar = self._file_.name.split(sep)[-1]
        pwd = os.getcwd()
        dir = self._get_rootdir()
        os.chdir(dir)
        runbader = False
        facfdat = os.path.join(dir, 'ACF.dat')
        if os.path.exists(facfdat):
            if os.stat(fchgcar)[-2] > os.stat(facfdat)[-2]:
                runbader = True
        else:
            runbader = True
        if runbader:
            Popen('bader ' + fchgcar, shell=True).wait()
        os.chdir(pwd)
        facfdat = open(facfdat)
        data = ''.join(facfdat.readlines()[2:-2])
        fdata = sio.StringIO(data)
        bader = np.loadtxt(fdata)
        charges = bader[:,4]

        return  charges


if __name__ == '__main__':

    import argparse, sys

    parser = argparse.ArgumentParser(
                description = '')

    

    print chgcar().read('Al/CHGCAR').get_magmom()
