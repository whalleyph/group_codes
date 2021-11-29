import numpy as np

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

if __name__ == '__main__':
  import doctest
  doctest.testmod()
