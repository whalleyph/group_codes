from utils import linesbetween
from cStringIO import StringIO as sio
import numpy as np
import tarfile
import os, sys
from latconv import latvec2latpar

class outcar:
    def __init__(self):
        self._file_ = 'OUTCAR'
        self.lines = ''

    def read(self, f=None):

        if tarfile.is_tarfile(f):
            num = f.split(os.path.sep)[-1].split('.')[0]
            outcar = '.'.join(['OUTCAR', num])
            foutcar = tarfile.open(f).extractfile(outcar)
            self.lines = foutcar.readlines()
        else:
            try:
                self.lines = f.readlines()
            except AttributeError:
                try:
                    self.lines = open(f).readlines()
                except:
                    print '%s found neither file nor string' % (f)
        
        self.lines = [line.strip() for line in self.lines]

        return self

        #if isinstance(f, file):
        #    self._file_ = f
        #elif isinstance(f, str):
        #    self._file_ = open(f)
        #else:
        #    print '%s is neither file nor string' % (f)
        #    sys.exit(-1)

        #if file is not None:
        #    try:
        #        file.pos
        #    except AttributeError:
        #        self._file_ = open(file)
        #    else:
        #        self._file_ = file
        #        self._file_.seek(0)
        #else:
        #    self._file_ = open(self._file_)

        self.lines = [line.strip() for line in self._file_.readlines()]
        print self.lines

        return self

    @property
    def params(self):
        rawincar = linesbetween('Dimension of arrays:',
        '----------------------------', self.lines)
    
        params = {}
        for line in rawincar:
            stringsplit = line.split()
        
            for i in range(len(stringsplit)-1):
                if stringsplit[i].isupper() and stringsplit[i+1] == '=':
                    params[stringsplit[i]] = stringsplit[i+2].replace(
                                                        ';','').strip()
                if (stringsplit[i].isupper() and stringsplit[i][-1] == '='):
                    params[stringsplit[i][:-1]] = stringsplit[i+1].replace(
                                                        ';','').strip()

            if line.split('=')[0].strip() == 'SYSTEM':
                params['SYSTEM'] = line.split('=')[1]

        for line in self.lines:
            if line.find('new method') > -1:
                params['LREAL'] = 'A'
                break

        return params
  
    def composition(self):
        def nions(self):
            for line in self.lines:
                if line.find('ions per type') > -1:
                    nions = [int(i) for i in line.split('=')[1].split()]
            return nions

        def tions(self):
            tions = [line.split()[3] for line in self.lines if 
                    line.find('TITEL') > -1]
            return tions
        
        if len(tions(self)) != len(nions(self)):
            print 'Screwed OUTCAR'
            sys.exit(-1)
        return zip (tions(self), nions(self))

    def nions(self):
        for line in self.lines:
            if line.find('ions per type') > -1:
                nions = [int(i) for i in line.split('=')[1].split()]
                break
        return nions

    def sumnions(self):
        for line in self.lines:
            if line.find('ions per type') > -1:
                nions = [int(i) for i in line.split('=')[1].split()]
                break
        return np.array(nions).sum()

    def tions(self):
        tions = [line.split()[3] for line in self.lines if 
                line.find('TITEL') > -1]
        return tions

    @property
    def e0(self):
        try:
            e0 = [float(line.split()[-1]) for line in self.lines if 
                    line.find('energy  without') > -1]
        except IndexError:
            print self._file_
            e0 = 0.0
        return e0

    def get_final_e0(self):
        try:
            fin_e0 = self.e0()[-1]
        except IndexError:
            fin_e0 = 0.0
        return fin_e0

    def vcell(self):
        volumes = [float(line.strip().split(':')[1]) for line in self.lines
                if line.find('volume of cell') > -1]
        del volumes[0]

        return volumes

    def get_volume(self):
        volumes = [float(line.strip().split(':')[1]) for line in self.lines
                if line.find('volume of cell') > -1]
        del volumes[0]

        return volumes

    def get_final_volume(self):
        return self.get_volume()[-1]

    def get_forces(self):
        forces = []
        nions = int(self.params['NIONS'])
        for i in range(len(self.lines)):
            if self.lines[i].find('POSITION') > -1:
                data = '\n'.join(self.lines[i+2:i+2+nions])
                data = sio(data)
                break

        forces = np.loadtxt(data, usecols=(3,4,5))

        return forces

    def get_positions(self):
        configurations = []
        nions = int(self.params['NIONS'])
        for i in range(len(self.lines)):
            if self.lines[i].find('POSITION') > -1:
                data = '\n'.join(self.lines[i+2:i+2+nions])
                data = sio(data)
                positions = np.loadtxt(data, usecols=(0,1,2))
                configurations.append(positions)

        return np.array(configurations)

    def get_scaled_positions(self):
        scaled_positions = []
        positions = self.get_positions()
        latvecs = self.get_latvec()
        for l, c in zip(latvecs, positions):
            d = np.linalg.solve(l.T, c.T).T
            scaled_positions.append(d)

        return np.array(scaled_positions)

    def get_mean_max_forces(self):
        forces = []
        nions = int(self.params['NIONS'])
        for i in range(len(self.lines)):
            if self.lines[i].find('POSITION') > -1:
                forces.append([line.strip().split()[-3:] for 
                    line in self.lines[i+2:i+2+nions]])

        forces = np.array(forces).astype(float)

        meanmax = []
        for i in range(len(forces)):
            mean = ((forces**2).sum(axis=2)**0.5)[i].mean()
            max  = ((forces**2).sum(axis=2)**0.5)[i].max()
            min  = ((forces**2).sum(axis=2)**0.5)[i].min()
            meanmax.append(np.hstack([mean,max,min]))

        return np.array(meanmax).astype(float)

    def pressure(self):
        for line in self.lines:
            if line.find('in kB') > -1:
                stress = [float(i) for i in line.split()[2:]]
            if line.find('external pressure') > -1:
                pressure = float(line.split()[3])
        return [pressure, stress]

    def magmom(self):
        magmom = []
        nions = int(self.params['NIONS'])
        for i in range(len(self.lines)):
            if self.lines[i].find('magnetization (x)') > -1:
                magmom.append([float(line.split()[-1]) for line in 
                    self.lines[i+4:i+4+nions]])

        return np.array(magmom)

    def get_magnetic_moments(self):
        return self.magmom()

    def efermi(self):
        efermi = [float(line.split()[2]) for line in 
                self.lines if line.find('E-fermi') > -1]
        if len(efermi) > 1:
            print 'Hope you know what you are doing!!\
            More than one E-fermi found'
        return efermi

    def get_latvec(self):
        latvecs = []
        for i, line in enumerate(self.lines):
            if line.find('direct lattice vectors') > -1:
                s = '\n'.join(self.lines[i+1:i+4])
                latvec = np.loadtxt(sio(s), usecols=(0,1,2))
                latvecs.append(latvec)
        
        return np.array(latvecs[1:])

    def get_latpar(self):
        latvec = self.get_latvec()
        return latvec2latpar(latvec)

    def get_elastic_moduli_total(self):
        for i, line in enumerate(self.lines):
            if line.find('TOTAL ELASTIC') > -1:
                s = '\n'.join(self.lines[i+3:i+9])
                break
        # Matrix from VASP
        # uses indices      6 in place of 4
        #                   4 in place of 5
        #                   5 in place of 6
        v = np.loadtxt(sio(s), usecols=(1,2,3,4,5,6))
        d = {0:0, 1:1, 2:2, 3:4, 4:5, 5:3}
        t = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                t[i,j] = v[d[i],d[j]]
        
        elastic_moduli_total= t

        return elastic_moduli_total

    def get_elastic_moduli_symmetric(self):
        for i, line in enumerate(self.lines):
            if line.find('SYMMETRIZED ELASTIC') > -1:
                s = '\n'.join(self.lines[i+3:i+9])
                break
        # Matrix from VASP
        # uses indices      6 in place of 4
        #                   4 in place of 5
        #                   5 in place of 6
        v = np.loadtxt(sio(s), usecols=(1,2,3,4,5,6))
        d = {0:0, 1:1, 2:2, 3:4, 4:5, 5:3}
        t = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                t[i,j] = v[d[i],d[j]]
        
        elastic_moduli_symmetrized= t

        return elastic_moduli_symmetrized

    def get_elastic_moduli_ionic(self):
        for i, line in enumerate(self.lines):
            if line.find('ELASTIC MODULI CONTR') > -1:
                s = '\n'.join(self.lines[i+3:i+9])
                break
        # Matrix from VASP
        # uses indices      6 in place of 4
        #                   4 in place of 5
        #                   5 in place of 6
        v = np.loadtxt(sio(s), usecols=(1,2,3,4,5,6))
        d = {0:0, 1:1, 2:2, 3:4, 4:5, 5:3}
        t = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                t[i,j] = v[d[i],d[j]]
        
        elastic_moduli_ionic = t

        return elastic_moduli_ionic

    def get_latpar_a(self):
        return self.get_latpar()[0]
    def get_latpar_b(self):
        return self.get_latpar()[1]
    def get_latpar_c(self):
        return self.get_latpar()[2]
    def get_latpar_alpha(self):
        return self.get_latpar()[3]
    def get_latpar_beta(self):
        return self.get_latpar()[4]
    def get_latpar_gamma(self):
        return self.get_latpar()[5]

    def get_vibrational_modes(self):
        frequencies = []
        position_modes = []
        for i, line in enumerate(self.lines):
            if line.find('f  =') > 1 or line.find('f/i=') > 1:
                frequencies.append(line)
                s = '\n'.join(self.lines[i+2:i+2+self.sumnions()])
                v = np.loadtxt(sio(s))
                position_modes.append(v)
        return zip(frequencies, position_modes)

