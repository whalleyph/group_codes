import copy
from StringIO import StringIO

test_kpoints = StringIO(
"""\
Test KPOINTS file
0
Monkhorst
3 9 9
0 0 0
"""
)

class kpoints:

    """
    comment = comment line
    genscheme = k-point generation scheme
    scheme = scheme used for generation (Monkhorst-Pack, etc.)
    kpoints = k-point grid
    shift = shift in the k-point grid

    >>> kpoints().read(test_kpoints).set_kpoints((4,11,11)).set_shift((0.5, 0.5, 0.5)).write('kpoints_testfile')
    """

    def __init__ (self):

        self._file_ = 'KPOINTS'
        self.comment = ''
        self.genscheme = '0'
        self.scheme = 'Monkhorst-Pack'
        self.kpoints = (1, 1, 1)
        self.shift = (0., 0., 0.)

    def read(self, file=None):
        """
        >>> kp = kpoints().read(test_kpoints)

        >>> print kp.comment
        Test KPOINTS file

        >>> print kp.genscheme
        0

        >>> print kp.scheme
        Monkhorst

        >>> print kp.kpoints
        (3, 9, 9)

        >>> print kp.shift
        (0.0, 0.0, 0.0)

        """
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

        lines = self._file_.readlines()
        lines = [line.strip() for line in lines]
        self.comment = lines[0]
        self.genscheme = lines[1]
        self.scheme = lines[2]
        self.kpoints = tuple (int(i) for i in lines[3].split()[:3])
        self.shift = tuple (float(i) for i in lines[4].split()[:3])
        return self

    def __str__ (self):
        return '\n'.join([self.comment, self.genscheme, self.scheme,
                      str(self.kpoints), str(self.shift)])

    def set_comment(self, comment):
        self.comment = comment
        return self

    def set_genscheme(self, genscheme):
        self.genscheme = genscheme
        return self

    def set_scheme(self, scheme):
        self.scheme = scheme
        return self

    def set_kpoints(self, kpoints):
        self.kpoints = kpoints
        return self

    def set_shift(self, shift):
        self.shift = shift
        return self

    def write (self, file, comment = None, genscheme = None, scheme = None, 
      kpoints = None, shift = None):

        """
        >>> kpoints().write('kpoints_testfile')
        """

        if genscheme:
          genscheme_str = genscheme
        else:
          genscheme_str = self.genscheme

        if kpoints:
          kpoints_str = ' '.join([str(i) for i in list(kpoints)])
        else:
          kpoints_str = ' '.join([str(i) for i in list(self.kpoints)])

        if comment:
          comment_str = comment
        elif kpoints:
          comment_str = 'x'.join([str(i) for i in list(kpoints)])
        else:
          comment_str = 'x'.join([str(i) for i in list(self.kpoints)])

        if scheme:
          if scheme.lower()[0] == 'm':
            scheme_str = 'Monkhorst-Pack'
          elif scheme.lower()[0] == 'g':
            scheme_str = 'Gamma'
        else:
          scheme_str = self.scheme

        if shift:
          shift_str = ' '.join([str(i) for i in list(shift)])
        else:
          shift_str = ' '.join([str(i) for i in list(self.shift)])
          
        kpoints = '\n'.join([comment_str,
                            genscheme_str, 
                            scheme_str,
                            kpoints_str,
                            shift_str]) + '\n'

        open(file, 'w').write(kpoints)

        return None

def genkpoints(a, b, c, n, maxkpoints=10000):
    """
    Same as genkpoints but all kpoint triplets are converted to 
    all odd or all even and returned.
    """

    def makealleven(list):
        for i in range(len(list)):
          if list[i]%2 == 1:
            list[i] += 1
        return list

    def makeallodd(list):
        for i in range(len(list)):
          if list[i]%2 == 0:
            list[i] += 1
        return list

    kptlst = []
    k1 = (0, 0, 0)
    for i in range(30,201):
        k2 = (int(round(i/a)), int(round(i/b)), int(round(i/c)))

        diffkpt = (k2 != k1)
        lessthanlimit = (n*k2[0]*k2[1]*k2[2] <= maxkpoints)

        if diffkpt and lessthanlimit:
          bitsum = k2[0]%2 + k2[1]%2 + k2[2]%2
          if bitsum == 1:
            k2 = tuple(makealleven(list(k2)))
          if bitsum == 2:
            k2 = tuple(makeallodd(list(k2)))

          kptlst.append(k2)
          k1 = copy.copy(k2)

    uniqkptlst = []
    for kpt in kptlst:
        if kpt not in uniqkptlst:
            uniqkptlst.append(kpt)

    return uniqkptlst

if __name__ == '__main__':

    import doctest, argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true')
    args = parser.parse_args()

    if args.test:
        doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)
