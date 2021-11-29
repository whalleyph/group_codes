from sys import exit
# import numpy as np

def str2dict(string, pairsep=',', keyvalsep='='):
  """
  Converts string to a dictionary based on a pair separator (pairsep) and a
  key value separator (keyvalsep)
  
  >>> str2dict('eggs are bad, carrots are good, milk are great', ',', 'are')
  {'EGGS': 'bad', 'CARROTS': 'good', 'MILK': 'great'}
  
  >>> print str2dict('isif=2, nsw=0, algo=fast')
  {'ISIF': '2', 'ALGO': 'fast', 'NSW': '0'}
  
  >>> print str2dict('isif=2/nsw=0/dontbe=stupid', pairsep='/')
  {'DONTBE': 'stupid', 'ISIF': '2', 'NSW': '0'}
  
  >>> print str2dict('isif:2/nsw:0/dontbe:stupid', pairsep='/', keyvalsep=':')
  {'DONTBE': 'stupid', 'ISIF': '2', 'NSW': '0'}
  
  Empty string returns an empty dictinary
  >>> print str2dict('')
  {}
  
  You have to use a 'pairsep'arated list of tag value pairs. Not just one word.
  >>> str2dict('justoneword') 
  Traceback (most recent call last):
    ...
  IndexError: list index out of range
  
  Try not to use space (if you do use be sure to use the right amount of space 
  everywhere) as pairsep or keyvalsep. It is hard (yet) to parse space as a 
  separator!!
  """
  d = {}
  if string:
    tagvals = string.strip().split(pairsep)
    for tagval in tagvals:
      d[tagval.strip().split(keyvalsep)[0].strip().upper()] = (
          tagval.strip().split(keyvalsep)[1].strip() )
  return d

def linesbetween(str1, str2, lines):
  """
  >>> lines = 'abc def ghi'.split()

  >>> linesbetween('', '', '')
  []

  >>> linesbetween('', 'ghi', lines)
  ['abc', 'def']

  >>> linesbetween('abc', 'ghi', lines)
  ['abc', 'def']

  >>> linesbetween('abc', 'ghi', [])
  []

  >>> linesbetween('abc', 'ghi', 'string')
  []
  """

  append = False
  chunk = []
  for line in lines:
    if line.find(str1) > -1:
      append = True
    if line.find(str2) > -1:
      append  = False
    if append:
      chunk.append(line)

  return chunk

def translate(t, p):

  """
  t = translation vector [t1 t2 t3]
  p = point/matrix to be translated in the form
  p = [x1 y1 z1
       x2 y2 z2
       x3 y3 z3
       x4 y4 z4]

  # Translating (2,2,2) by (0,0.0)
  >>> translate([0,0,0], [2,2,2])
  array([[ 2.,  2.,  2.]])

  # Translating the points (2,2,2) and (3,3,3) by (0,0,0)
  >>> translate([0,0,0], [[2,2,2],[3,3,3]])
  array([[ 2.,  2.,  2.],
         [ 3.,  3.,  3.]])

  # Translating the points (2,2,2) and (3,3,3) by (1,0,1)
  >>> translate([1,0,1], [[2,2,2],[3,3,3]])
  array([[ 3.,  2.,  3.],
         [ 4.,  3.,  4.]])

  # Returns an empty array if passed an empty list of coordinates
  >>> translate([1,1,1], [])
  []
  
  >>> a = np.array([])
  >>> translate([1,1,1], a)
  []

  """

  if list(p) != []:
    t = np.mat(t).transpose()
    p = np.mat(p).transpose()

    T = np.mat(np.identity(4, float))
    T[:3,3] = t

    prows, pcols = p.shape

    P = np.mat(np.ones((prows+1, pcols), float))
    P[:3] = p

    return np.array(np.dot(T, P)[:-1].transpose())
  else:
    return []

"""Replicates the gvien coordinates to the given super-cell"""

def replicate (supercell, coordinates):

  """
  Builds a 'supercell' of 'coordinates'
  supercell is a tuple of three numbers (n1, n2, n3)
  coordinates is a (nx3) array of coordinates in the form
  (x1, y1, z1) ... (xn, yn, zn)

  >>> coord = [0.0, 0.0, 0.0]
  >>> supercell = (2, 2, 2)
  >>> replicate (supercell, coord)
  array([[ 0.,  0.,  0.],
         [ 0.,  0.,  1.],
         [ 0.,  1.,  0.],
         [ 0.,  1.,  1.],
         [ 1.,  0.,  0.],
         [ 1.,  0.,  1.],
         [ 1.,  1.,  0.],
         [ 1.,  1.,  1.]])

  # Returns an empty list if the input list of coordinates is empty
  >>> coord = []
  >>> supercell = (2, 2, 1)
  >>> replicate (supercell, coord)
  []

  >>> coord = [[1.0, 1.0, 1.0], [0.5, 0.5, 0.5]]
  >>> supercell = (2, 2, 2)
  >>> replicate (supercell, coord)
  array([[ 1. ,  1. ,  1. ],
         [ 0.5,  0.5,  0.5],
         [ 1. ,  1. ,  2. ],
         [ 0.5,  0.5,  1.5],
         [ 1. ,  2. ,  1. ],
         [ 0.5,  1.5,  0.5],
         [ 1. ,  2. ,  2. ],
         [ 0.5,  1.5,  1.5],
         [ 2. ,  1. ,  1. ],
         [ 1.5,  0.5,  0.5],
         [ 2. ,  1. ,  2. ],
         [ 1.5,  0.5,  1.5],
         [ 2. ,  2. ,  1. ],
         [ 1.5,  1.5,  0.5],
         [ 2. ,  2. ,  2. ],
         [ 1.5,  1.5,  1.5]])
  """

  if list(coordinates) != []:
    newcoord = []
    for index in np.ndindex(supercell):
      newcoord.append((translate(index, coordinates)))

    tmp = np.ravel(newcoord)

    return (tmp.reshape(len(tmp)/3, 3))
  else:
    return []

def car2dir(A, c):
  d = np.dot(c, np.linalg.inv(A))
  return d

def dir2car(A, d):
  c = np.dot(d, A)
  return c

def minimg(a, a0):
  dr = a - a0
  dr = dr - np.around(dr)
  return dr

def center_around (configuration=None, point=None, tcoord='direct'):
  from poscar import poscar 

  latvec = configuration.cellvec
  print latvec
  ids    = configuration.ids
  if tcoord == 'cartesian':
    cartcoord = configuration.coord
    point_direct = car2dir(latvec, point)
    directcoord = car2dir(latvec, cartcoord)
  elif tcoord == 'direct':
    point_direct = point
    directcoord = configuration.coord
  else:
    print "%s %s" % ('Unkown type of coordinates', tcoord)
    raise SystemExit
  centered_direct = minimg(directcoord, point_direct)
  centered_cart = dir2car(latvec, centered_direct)
  distances = (centered_cart**2).sum(axis=1)**0.5
  sorted_ions = distances.argsort(kind='mergesort')
  centered_config = poscar()
  centered_config.cellvec = latvec
  sorted_ids = []
  sorted_coord = []
  sorted_distances = []
  for i in sorted_ions:
    sorted_ids.append(ids[i])
    sorted_coord.append(centered_cart[i])
    sorted_distances.append(distances[i])
    #formatted_coord = "%7.3f %7.3f %7.3f" % tuple(centered_cart[i])
    #print "%2.2i %5.3f %6s %s" % (i, distances[i], ids[i], formatted_coord)
  centered_config.coord = np.array(sorted_coord)
  centered_config.ids = sorted_ids
  return zip(centered_config.ids, centered_config.coord, sorted_distances, sorted_ions), centered_config

def center_around_333 (poscar, point=None, tcoord='direct'):
  point = point/3.0
  supercell = poscar.replicate((3,3,3))
  print supercell.cellvec
  return center_around (supercell, point, tcoord)

def strippv(eles):
    stripped_list = []
    for ele in eles:
        stripped_list.append(ele.split('_')[0])
    return stripped_list

def unarchive(tgz, file):
    import tarfile
    if tarfile.is_tarfile(tgz):
        ar = tarfile.open(tgz)
    else:
        print '%s is not an archive' % (tgz)
        exit(1)

    return ar.extractfile(file)

def unarchive_calc(dir):
    import os
    import glob
    import tarfile

    xtgz = [i.split(os.path.sep)[-1] for i in glob.glob(os.path.join(dir, '?.tgz'))]
    xxtgz = [i.split(os.path.sep)[-1] for i in glob.glob(os.path.join(dir, '??.tgz'))]

    print xtgz, xxtgz

    if xxtgz:
        tgz = os.path.join(dir, sorted(xxtgz)[-1])
    else:
        tgz = os.path.join(dir, sorted(xtgz)[-1])

    tar = tarfile.open(tgz)
    tar.extractall(path=dir)
    files_in_archive = tar.getnames()
    for file in files_in_archive:
        src = os.path.join(dir, file)
        dst = os.path.join(dir, file[:-2] + '_FINAL')
        if os.path.exists(dst):
            os.remove(dst)
            os.symlink(src, dst)
        else:
            os.symlink(src, dst)
    
    return None

if __name__ == '__main__':
  import doctest
  doctest.testmod()
