import numpy as np
from translate import translate

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
        print index
    for index in np.ndindex(supercell):
      newcoord.append((translate(index, coordinates)))

    tmp = np.ravel(newcoord)

    return (tmp.reshape(len(tmp)/3, 3))
  else:
    return []

if __name__ == '__main__':
  import doctest
  doctest.testmod()
