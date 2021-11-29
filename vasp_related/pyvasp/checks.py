import os.path
from StringIO import StringIO
def completion(file):
  """Check if the calclation has completed successfully.
  Does not matter if the calculation has not converged.
  Checks vasp.out for the words writing wavefunctions.
  
  >>> completion('test_files/vasp.out')
  True
  
  >>> completion('test_files/nofile')
  File does not exist ...

  >>> completion('test_files/INCAR')
  False
  
  >>> completion()
  Traceback (most recent call last):
  ...
  TypeError: completion() takes exactly 1 argument (0 given)
  """

  if file:
      try:
        vaspout = open(file)
      except IOError:
        print "File does not exist at %s" % (file)
      else:
        if vaspout.read().find('User time') > -1:
          return True
        else:
          return False

def convergence(file):
    """Check for convergence from vasp.out.
    Returns True if it finds the word required
    
    >>> convergence('test_files/nofile')
    File does not exist at test_files/nofile
    >>> print convergence('test_files/vasp.out')
    True
    >>> convergence('test_files/INCAR')
    False
    
    >>> convergence()
    Traceback (most recent call last):
    ...
    TypeError: convergence() takes exactly 1 argument (0 given)
    """

    if file:
        try:
          vaspout = open(file)
        except IOError:
          print "File does not exist at %s" % (file)
        else:
          if vaspout.read().find('required') > -1:
            return True
          else:
            return False

def tgzconvergence(tgz):
    import tarfile as tar
    index = tgz.split(os.path.sep)[-1].replace('.tgz', '')
    vaspout = 'vasp.out.' + index
    tarfile = tar.open(tgz)
    if tarfile.extractfile(vaspout).read().find('req') > -1:
        return True
    else:
        return False

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags = doctest.ELLIPSIS)
