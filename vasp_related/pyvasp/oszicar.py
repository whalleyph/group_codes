import os
import tarfile

class oszicar:
    def __init__(self):
        self._file_ = 'OSZICAR'
        self.lines = ''

    def read(self, f=None):

        if tarfile.is_tarfile(f):
            num = f.split(os.path.sep)[-1].split('.')[0]
            outcar = '.'.join(['OSZICAR', num])
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
    def e0(self):
        return [float(line.split()[4]) for line in self.lines if line.find('F=') > -1]

    def get_e0(self):
        return self.e0()

    def get_energies(self):
        return self.e0()

    def get_final_energy(self):
        return self.e0()[-1]

    def totmagmom (self):
        return [float(line.split()[-1]) for line in self.lines if line.find('mag=') > -1]

    def get_total_magnetic_moment (self):
        return [float(line.split()[-1]) for line in self.lines if line.find('mag=') > -1]

    def get_temperatures(self):
        temp = []
        for line in self.lines:
            if line.find('T=') > -1:
                try:
                    temp.append(float(line.split()[2]))
                except ValueError: 
                    temp.append(999999.0)

        return temp
