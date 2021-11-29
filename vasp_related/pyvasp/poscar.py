#! /bin/env python

import os
import numpy as np
from sys import stderr
from copy import copy, deepcopy
from StringIO import StringIO as sio
from tempfile import NamedTemporaryFile as tmpfile

testposcar = sio(
"""\
C Fe Fe Cr oP16 M3C Cementite                      
   1.00000000000000     
     5.1188776123521880    0.0000000000000000   -0.0740365002817798
     0.0000000000000000    6.8414705451276365    0.0000000000000000
    -0.0685192778490072    0.0000000000000000    4.4806853108558959
   4   8   3   1
Direct
  0.7355543135359121  0.3579794705882335  0.5733273659914344
  0.0158074593425641  0.8579793823529442  0.0566685773561339
  0.2408464098378537  0.3579794705882335  0.1952712457326108
  0.4960104305236720  0.8579793823529442  0.7509391029459945
  0.7995046463170539  0.6703278689514439  0.7974745595350671
  0.7995046853795529  0.0456310134014922  0.7974745595350671
  0.9432483498194049  0.1812797350809912  0.3092652589883943
  0.9432484435694082  0.5346790590366557  0.3092652589883943
  0.2844784872042710  0.0453262948261171  0.9650777403434906
  0.2844784872042710  0.6706324728209493  0.9650777927452500
  0.4343502046396406  0.5338348385508587  0.4725855698490415
  0.4343502046396406  0.1821239555667882  0.4725855698490415
  0.0774692983306396  0.3579794705882335  0.8001407909762037
  0.1575615073556639  0.8579793823529442  0.4777199740708724
  0.5817300746639138  0.3579794705882335  0.9767932520953404
  0.6447572183397177  0.8579793823529442  0.2862577805609648
 
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
""")

class poscar:
  def __init__ (self):

    self._file_ = 'POSCAR'
    self.comment = ''
    self.scale   = 1.0
    self.latvec  = np.identity(3, float)
    self.cellvec = self.scale*self.latvec
    self.nions   = np.ones(1, int)
    self.selective = False
    self.coordtype  = 'Direct'
    self.coord   = np.zeros((3), float)
    self.velocities = np.zeros((3), float)

    self.types   = ['Fe']
    self.ids     = ['Fe']

  def read(self, file=None):
      """
      >>> pos = poscar().read(testposcar)

      >>> pos.comment
      'C Fe Fe Cr oP16 M3C Cementite'

      >>> pos.scale
      1.0
      >>> pos.types
      ['C', 'Fe', 'Fe', 'Cr']

      >>> pos.coord
      array([[...]])

      >>> pos.nions
      array([4, 8, 3, 1])

      >>> pos.latvec
      array([[...]])

      >>> pos.cellvec
      array([[...]])

      >>> poscar().read('nonexistent_file')
      Traceback (most recent call last):
      ...
      IOError: [Errno 2] No such file or directory: 'nonexistent_file'

      """

      if file is not None:
          try:
              file.pos
          except AttributeError:
              self._file_ = open(file)
          else:
              self._file_ = file
      else:
          self._file_ = open(self._file_)

      lines = self._file_.readlines()
      lines = [line.strip() for line in lines]

      self.comment = lines[0]
      self.scale   = float(lines[1])

      self.latvec = np.loadtxt(sio('\n'.join(lines[2:5])))

      self.cellvec  = self.scale*self.latvec

      if lines[6][0].lower() == 's' or lines[7][0].lower() == 's':
          self.selective = True
      else:
          self.selective = False

      if self.selective and lines[8][0].lower() in ['c', 'k', 'd']:
          version = '5'
          startcoord = 9
          if lines[8][0].lower() in 'c k'.split():
              coordtype = 'Cartesian'
          elif lines[8][0].lower() == 'd':
              coordtype = 'Direct'
      elif self.selective and lines[7][0].lower() in ['c', 'k', 'd']:
          version = '4'
          startcoord = 8
          if lines[7][0].lower() in 'c k'.split():
              coordtype = 'Cartesian'
          elif lines[7][0].lower() == 'd':
              coordtype = 'Direct'
      elif not self.selective and lines[7][0].lower() in ['c', 'k', 'd']:
          version = '5'
          startcoord = 8
          if lines[7][0].lower() in 'c k'.split():
              coordtype = 'Cartesian'
          elif lines[7][0].lower() == 'd':
              coordtype = 'Direct'
      elif not self.selective and lines[6][0].lower() in ['c', 'k', 'd']:
          version = '4'
          startcoord = 7
          if lines[6][0].lower() in 'c k'.split():
              coordtype = 'Cartesian'
          elif lines[6][0].lower() == 'd':
              coordtype = 'Direct'

      if version == '4':
        self.nions = np.array(lines[5].split()).astype(int)
      elif version == '5':
        self.nions = np.array(lines[6].split()).astype(int)

      try:
        self.comment.split()[len(self.nions)]
      except IndexError:
        self.types = []
      else:
        self.types = self.comment.split()[:len(self.nions)]

      self.ids = []
      if self.types:
        for i in range(len(self.nions)):
          self.ids.extend(((self.types[i]+' ')*self.nions[i]).split())

      #if lines[6].lower()[0] == 's':
      #  self.selective = True
      #  coordtype = lines[7]
      #  startcoord = 8
      #else:
      #  self.selective = False
      #  coordtype = lines[6]
      #  startcoord = 7

      #if coordtype.lower()[0] in 'c k'.split():
      #  self.coordtype = 'Cartesian'
      #elif coordtype.lower()[0] == 'd':
      #  self.coordtype = 'Direct'
      #else:
      #  print >> stderr, 'Unkown coordinate type %s' % (coordtype)
      #  sys.exit(1)

      if len(self.nions) == 1 and self.nions == 1:
          self.coord = np.loadtxt(sio(
              '\n'.join(lines[startcoord:startcoord+self.nions.sum()])),
              dtype=float, usecols=(0,1,2))
          self.coord = np.array(self.coord, ndmin=2)
      else:
          self.coord = np.loadtxt(sio(
              '\n'.join(lines[startcoord:startcoord+self.nions.sum()])),
              dtype=float, usecols=(0,1,2))

      return self

  def _str_ (self):
    comment_str = "%-40s" % (self.comment)

    scale_str = "%19.14f" % float(self.scale)

    latvec_str = ''
    for i in self.latvec:
      latvec_str += " %22.16f%22.16f%22.16f" % tuple(
          [j for j in i]) + '\n'
    latvec_str = latvec_str.rstrip()

    nions_str = ''.join(["%4i" % (i) for i in self.nions])

    if self.selective:
      selective_str = 'Selective Dynamics'

    coordtype_str = self.coordtype

    coord_str = ''
    for i in self.coord:
      coord_str += ''.join(["%20.16f" % (j) for j in i]) + '\n'
    coord_str = coord_str.rstrip()

    if self.selective:
      poscar = '\n'.join([
        comment_str,
        scale_str,
        latvec_str,
        selective_str,
        nions_str,
        coordtype_str,
        coord_str])
    else:
      poscar = '\n'.join([
        comment_str,
        scale_str,
        latvec_str,
        nions_str,
        coordtype_str,
        coord_str])

    return poscar

  def __str__ (self):
    return self._str_()

  def write (self, file='POSCAR'):
    writefile = False

    if os.path.exists(file):
      print 'File exists at %s' % file
      yes = raw_input('Overwrite file? [No]')
      if yes.startswith('y') or yes.startswith('Y'):
        writefile = True
    else:
      writefile = True

    if writefile:
      open(file, 'w').write(self._str_()+'\n')
    return None

  def set_scale(self, scale):
      self.scale = scale
      return self

  def get_nions(self):
    return self.nions

  def set_iontypes(self, iontypes):
    self.types = iontypes
    self.ids = []
    for i in range(len(self.nions)):
      self.ids.extend(((self.types[i]+' ')*self.nions[i]).split())
    return self

  def set_coordinate(self, iion, xyz):
    """Set coordinates of the ith ion with xyz given as a list of [x, y, z]
    """
    self.coord[iion] = xyz
    return self

  def set_comment(self, comment):
    self.comment = comment
    return self
  
  def set_nions(self, nions):
      self.nions = np.array(nions)
      return self

  def del_coordinate(self, iion):
    self.coord = np.vstack([self.coord[:iion], self.coord[iion+1:]])
    print >> stderr, 'Not changing nions!! Change manually'
    return self

  def replicate(self, supercell):

    from utils import replicate
    #from operator import itemgetter

    """Builds a supercell of a POSCAR file
    supercell is a tuple of the form (n1, n2, n3) where n1, n2, n3 are the
    number of cells required in the x, y and z directions"""

    sposcar = deepcopy(self)
    ncells = supercell[0] * supercell[1] * supercell[2]
    sposcar.ids = []
    sposcar.nions = []

    sposcar.nions = [i*ncells for i in self.nions]
    if self.types:
      for i in range(len(sposcar.nions)):
        sposcar.ids.extend(((self.types[i]+' ')*sposcar.nions[i]).split())

    sposcar.scale = 1.0

    print supercell

    sposcar.latvec[0,:] = self.cellvec[0,:]*supercell[0]
    sposcar.latvec[1,:] = self.cellvec[1,:]*supercell[1]
    sposcar.latvec[2,:] = self.cellvec[2,:]*supercell[2]

    sposcar.cellvec[0,:] = self.cellvec[0,:]*supercell[0]
    sposcar.cellvec[1,:] = self.cellvec[1,:]*supercell[1]
    sposcar.cellvec[2,:] = self.cellvec[2,:]*supercell[2]

    if self.coordtype == 'Direct':
      indices = [sum(self.nions[:i+1]) for i in range(len(self.nions))]
      indices[0:0] = [0]

      coord = []
      for i in range(len(indices)-1):
        tmpcoord = replicate(supercell, self.coord[indices[i]:indices[i+1]])

        x = tmpcoord[:,0]/float(supercell[0])
        y = tmpcoord[:,1]/float(supercell[1])
        z = tmpcoord[:,2]/float(supercell[2])

        coord.append((np.array([x,y,z]).transpose()))

      sposcar.coord = np.vstack(coord)
      return sposcar
    else:
      print >> stderr, 'Hmm... needs more work then'
      sys.exit(1)

  def get_spacegroup(self):
    import poscar2findsym as p2f
    return p2f.findsym().read(self._file_).get_spacegroup()


if __name__ == '__main__':
    import doctest, argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--test', action='store_true')

    args = parser.parse_args()

    if args.test:
        doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)
