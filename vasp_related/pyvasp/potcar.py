#! /bin/env python

import os
from sys import stderr
from StringIO import StringIO

test_paw = StringIO(

"""\
   TITEL  = PAW_GGA C 05Jan2001
   TITEL  = PAW_GGA Fe 03Mar1998
   TITEL  = PAW_GGA Fe 03Mar1998
   TITEL  = PAW_GGA Nb_pv 09Jan2002
""")

class potcar:
    def __init__ (self):

        self._file_ = 'POTCAR'
        self.ions = []
        self.nelectrons = []
        self.xc = ''


    def read(self, file=None):

        """
        >>> print potcar().read(test_paw)
        pw91 ['C', 'Fe', 'Fe', 'Nb_pv']

        >>> print potcar().read('test_files/POTCAR')
        pw91 ['C', 'Fe', 'Fe', 'Nb_pv']


        >>> print potcar().read('nonexistent_file')
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
                self._file_.seek(0)
        else:
            self._file_ = open(self._file_)

        lines = self._file_.readlines()

        # self.xc: get the exchange correlation function used
        for line in lines:
            if line.find('TITEL') > -1:
                linesplit = line.split()
                if linesplit[2] == 'PAW_GGA':
                    self.xc = 'pw91'
                elif linesplit[2] == 'PAW_PBE':
                    self.xc = 'pbe'
                elif linesplit[2] == 'PAW':
                    self.xc = 'lda'
                elif linesplit[2] == 'US':
                    for anotherline in lines:
                        if anotherline.find('gradient') > -1:
                            self.xc = 'usgga'
                        else:
                            self.xc = 'uslda'
                elif linesplit[2] == 'oo':
                    self.xc = 'oo'

                # self.ions: collect iontypes
                self.ions.append(linesplit[3])

            # self.nelectrons: collect the valence electrons for each
            #                  ion type
            if line.find('ZVAL') > -1:
                # POMASS = 28.085; ZVAL = 4.000 mass and valenz
                no_valence_electrons = float(line.split()[5])
                self.nelectrons.append(no_valence_electrons)

        return self


    def __str__ (self):
        return ' '.join([self.xc, str(self.ions)])

    def get_iontypes(self):
        return self.ions

    def write (self, ions, xc, filename='POTCAR'):

        """
        >>> potcar().write(ions = ['H', 'O'], xc = 'pw91')

        >>> ions = 'C Fe Fe'.split()
        >>> potcar().write(ions=ions, xc='usgga')

        >>> potcar().write()
        Traceback (most recent call last):
        ...
        TypeError: write() takes at least 3 arguments (1 given)

        # Writes a blank POTCAR file
        >>> potcar().write(ions = [], xc = '')
        """

        if list:
            self.ions = ions
        if xc:
            self.xc = xc
        potcar = ''
        envvar = 'VASP_PP_PATH'
        if envvar in os.environ and os.environ[envvar] != '':
            root = os.environ[envvar]
        else:
            print 'Looks like %s is not set in ~/.bashrc' % envvar

        if self.xc.lower() == 'pw91':
            pproot = 'potpaw_GGA'
        elif self.xc.lower() == 'pbe':
            pproot = 'potpaw_PBE'
        elif self.xc.lower() == 'lda':
            pproot = 'potpaw'
        elif self.xc.lower() == 'uslda':
            pproot = 'pot'
        elif self.xc.lower() == 'usgga':
            pproot = 'pot_GGA'

        for ele in self.ions:
            pp = os.path.join(root, pproot, ele, 'POTCAR.Z')
            if os.path.exists (pp):
                potcar = potcar + os.popen('zcat '+ pp).read()
            else:
                print 'Pseudopotential for %s does not eixist at %s' % (
                        ele, pp)

        open (os.path.join(filename), 'w').write(potcar)

        return None

    def get_nelectrons(self):
        return self.nelectrons


if __name__ == '__main__':

    import doctest, argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--test', action='store_true')

    args = parser.parse_args()

    if args.test:
        doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)
