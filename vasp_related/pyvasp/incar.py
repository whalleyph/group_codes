#! /bin/env python

import os
from sys import stderr
from utils import str2dict
from StringIO import StringIO

test_incar = StringIO(
"""\
ALGO = Fast
EDIFFG = -0.05
ENCUT = 400
IBRION = 2
ICHARG = 1
ISIF = 3
ISMEAR = 1
ISPIN = 2
ISTART = 1
LORBIT = 11
LREAL = Auto
MAGMOM = 4*0.0 4*3.0 8*3.0
NSW = 60
POTIM = 0.1
PREC = accurate
SIGMA = 0.1
SYSTEM = C4Fe7CrFe4 Fe11MC4
VOSKOWN = 1
""")


class incar:
    def __init__(self):

        self._file_ = 'INCAR'
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

        lines = self._file_.readlines()

        lines = [line.strip() for line in lines 
                    if not line.strip().startswith('#') and line.strip() != '']
        lines = [line.strip() for line in lines 
                    if '=' in line]
        # print lines

        for line in lines:
          linesplit = line.strip().split('=')
          self.tags[linesplit[0].strip()] = linesplit[1].strip()

        return self
  
    def __str__ (self):
        string = ''
        keys = sorted(self.tags.keys())
        for key in keys:
            string += (' = '.join([key.upper(), self.tags[key]])+'\n')
        return string.strip()

    def write(self, file='INCAR'):
        inc = open(file, 'w')
        keys = sorted(self.tags.keys())
        for key in keys:
            inc.write(' = '.join([key.upper(), self.tags[key]])+'\n')
        inc.close()
        return None

    def settag(self, string=None, dictionary=None, 
                   pairsep=',', keyvalsep='='):

        """
        >>> print incar().read(test_incar).settag('lplane = .TRUE.')
        ALGO = Fast
        ...
        LPLANE = .TRUE.
        ...
        NSW = 60
        ...

        >>> print incar().read(test_incar).settag('lplane = .TRUE., nsw=49')
        ALGO = Fast
        ...
        LPLANE = .TRUE.
        ...
        NSW = 49
        ...
        
        >>> print incar().read(test_incar).settag('')
        ALGO = Fast
        ...

        >>> print incar().read(test_incar).settag(dictionary={'isif': '96', 'nsw': '900'})
        ALGO = Fast
        ...
        ISIF = 96
        ...
        NSW = 900
        ...
        """

        if string is None and dictionary is None:
            print '''
            You have to supply either a string, dictionary!
            Check help
            '''
            return

        if string is not None:
            d = str2dict(string, pairsep, keyvalsep)

        if dictionary:
            d = dictionary

        for key in d:
            self.tags[key.upper()] = d[key]
        return self

    def settagk(self, **kwargs):

        """
        >>> print incar().read(test_incar).settagk(lplane = '.TRUE.')
        ALGO = Fast
        ...
        LPLANE = .TRUE.
        ...
        VOSKOWN = 1

        >>> incar().read(test_incar).settagk('')
        Traceback (most recent call last):
        ...
        TypeError: settagk() takes exactly 1 argument (2 given)
        """

        for key in kwargs:
            self.tags[key.upper()] = kwargs[key]
        return self

    def unsettag(self, string=None):

        """
        >>> print incar().read(test_incar).unsettag('algo voskown')
        EDIFFG = -0.05
        ...
        SYSTEM = C4Fe7CrFe4 Fe11MC4

        >>> print incar().read(test_incar).unsettag('')
        ALGO = Fast
        ...
        """

        if string:
            args = string.split()
            for key in args:
                if self.tags.has_key(key.upper()):
                    del (self.tags[key.upper()])
        return self

if __name__ == '__main__':

    import argparse, sys

    parser = argparse.ArgumentParser(
                description = 
    """
    Manipulate INCAR from command line
    Examples:

    To set tags:
    incar.py -s 'nsw = 40, ispin = 2, encut = 400'

    To unset tags:
    incar.py -u 'lorbit isif'

    To set and unset tags from the file 4c/Al/INCAR:
    incar.py -f 4c/Al/INCAR -s 'encut = 400' -u 'ibrion lorbit'
    """)

    parser.add_argument('-s', '--set', type=str,
        help="Comma separated key = value pairs like 'ispin = 2, nsw = 20'")
    parser.add_argument('-u', '--unset', type=str,
        help="Space separated list of tags like 'ispin lorbit'")
    parser.add_argument('-f', '--file', type=str)
    parser.add_argument('-t', '--test', action='store_true')
    args = parser.parse_args()

    if args.test:
        import doctest
        doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)

    if args.set and args.unset:
        print incar().read(args.file).settag(string=args.set).unsettag(
          string=args.unset)
    elif args.set:
        print incar().read(args.file).settag(string=args.set)
    elif args.unset:
        print incar().read(args.file).unsettag(string=args.unset)
