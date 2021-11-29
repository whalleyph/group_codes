import os
import shutil
import logging

class VaspCalculation:
    def __init__(self, directory = '.'):
        self.directory = directory
        self.files_required = 'INCAR POSCAR POTCAR KPOINTS'.split()
        self.incar, self.poscar, self.potcar, self.kpoints  = \
                [os.path.join(self.directory, i) for i in self.files_required]
        self.file_map = dict([('INCAR', 'INCAR'), ('POSCAR', 'POSCAR'), ('POTCAR', 'POTCAR'), ('KPOINTS', 'KPOINTS')])

        ## include KPOINTS if KPOINTS file exists in the directory
        #if os.path.exists(os.path.join(self.directory, 'KPOINTS')):
        #        self.files_required += ['KPOINTS']
        #        self.kpoints = os.path.join(self.directory, 'KPOINTS')
        #        self.file_map['KPOINTS'] =  'KPOINTS'

        # ignore KPOINTS if KSPACING exists in INCAR
        #if 'kspacing' in open(self.incar).read().lower():
        #    self.files_required.remove('KPOINTS')
        #    del self.file_map['KPOINTS']
        #    del self.kpoints

        self.continue_from_directory = '.'

    def read(self):
        pass

    def write(self):
        if self.directory != self.continue_from_directory:
            os.mkdir(self.directory)
            for f in self.files_required:
                src = os.path.join(self.continue_from_directory, f)
                file_root = src.split(os.path.sep)[-1]
                dst = os.path.join(self.directory, self.file_map[file_root])
                logging.info('copying %s to %s' % (src, dst))
                try:
                    shutil.copy(src, dst)
                except:
                    logging.info('%s does not exist, copying %s failed' % (src, src))
        elif self.directory == self.continue_from_directory:
            for f in self.files_required:
                src = os.path.join(self.continue_from_directory, f)
                file_root = src.split(os.path.sep)[-1]
                dst = os.path.join(self.directory, self.file_map[file_root])
                #  move file only if non-empty
                if os.path.exists(src) and open(src).read().strip():
                    try:
                        shutil.move(src, dst)
                        logging.info('moving %s to %s' % (src, dst))
                    except:
                        pass
                else:
                    logging.info('File empty: %s. Not moving.' % (src))
            else:
                logging.info('CONTCAR does not exist, doing nothing')
        return None

    def continue_from(self, vaspcalculation):
        self.continue_from_directory = vaspcalculation.directory
        src_contcar = os.path.join(self.continue_from_directory, 'CONTCAR')
        if os.path.exists(src_contcar):
            self.files_required.remove('POSCAR')
            self.files_required += ['CONTCAR']
            self.file_map['CONTCAR'] = 'POSCAR'
        else:
            logging.info('CONTCAR does not exist in %s. Using POSCAR.' % (self.continue_from_directory))
        return self

    def apply_delta_incar(self, delta_incar):
        # Read INCAR and convert it to a dictionary
        # Will fail if there are comments in the INCAR
        incar = {}
        incar_lines = open(self.incar).readlines()
        lines = [i.strip() for i in incar_lines if i.strip().find('=')>0 and i.strip()[0] != '#']
        for line in lines:
            key, value = [i.strip() for i in line.strip().split('=')]
            incar[key] = value
        
        # Read all the changes to be applied
        # Keys to be removed start with '-'
        # Keys to be changed start with '>'
        # Kesy to be added start with '+'
        changes = delta_incar.strip().split('\n')
        remove_keys = [i for i in changes if i.strip()[0] == '-']
        change_keys = [i for i in changes if i.strip()[0] in ['>', '+']]

        # Remove keys to be removed
        for key in remove_keys:
            k = key.replace('-', '').strip()
            if incar.get(k):
                del incar[k]

        # Change keys to be changed and add keys to be added
        # Using a dictionary implies changes and additions can be handled
        # similarly
        for keyval in change_keys:
            k, v = keyval.replace('+', '').replace('>', '').split('=')
            k, v = k.strip(), v.strip()
            incar[k] = v

        # Build the new INCAR string
        incar_string = ''
        for key in sorted(incar.keys()):
            incar_string += ' = '.join([key, incar[key]])+'\n'

        open(self.incar, 'w').write(incar_string)

        return self

    def apply_delta_kpoints(self, delta_kpoints):
        kpts = open(self.kpoints).readlines()
        kpts[0] = 'x'.join(delta_kpoints.split()) + '\n'
        kpts[3] = delta_kpoints + '\n'
        open(self.kpoints, 'w').write(''.join(kpts))

        return self

    def submit():
        pass

class Simple(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)

class Meta_GGA(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['WAVECAR']
        self.file_map['WAVECAR'] = 'WAVECAR'
        self.wavecar = os.path.join(directory, 'WAVECAR')

class Simple_band(Simple):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['CHGCAR']
        self.file_map['CHGCAR'] = 'CHGCAR'
        self.chgcar = os.path.join(directory, 'CHGCAR')

class Hybrid(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['WAVECAR']
        self.file_map['WAVECAR'] = 'WAVECAR'
        self.wavecar = os.path.join(directory, 'WAVECAR')

class Optics(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['WAVECAR']
        self.file_map['WAVECAR'] = 'WAVECAR'
        self.wavecar = os.path.join(directory, 'WAVECAR')

class GW(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['WAVECAR', 'WAVEDER']
        self.file_map['WAVECAR'] = 'WAVECAR'
        self.file_map['WAVEDER'] = 'WAVEDER'
        self.wavecar = os.path.join(directory, 'WAVECAR')
        self.waveder = os.path.join(directory, 'WAVEDER')

class Hybrid_band(Hybrid):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['CHGCAR']
        self.file_map['CHGCAR'] = 'CHGCAR'
        self.chgcar = os.path.join(directory, 'CHGCAR')

class Dimer(VaspCalculation):
    def __init__(self, directory = '.'):
        VaspCalculation.__init__(self, directory)
        self.files_required += ['MODECAR']
        self.file_map['MODECAR'] = 'MODECAR'
        self.modecar = os.path.join(directory, 'MODECAR')

    def continue_from(self, vaspcalculation):
        VaspCalculation.continue_from(self, vaspcalculation)
        src_modecar = os.path.join(self.continue_from_directory, 'NEWMODECAR')
        if os.path.exists(src_modecar):
            self.files_required.remove('MODECAR')
            self.files_required += ['NEWMODECAR']
            self.file_map['NEWMODECAR'] = 'MODECAR'
        else:
            logging.info('NEWMODECAR does not exist in %s Using MODECAR.' % (self.continue_from_directory))
        return self


#class relaxation(VaspCalculation):
#    def __init__(self):
#        VaspCalculation.__init__(self)
#
#class dftd3(VaspCalculation):
#    def __init__(self):
#        VaspCalculation.__init__(self)

if __name__ == '__main__':
    normal = Simple()
    accurate = Simple().continue_from(normal)
    dimer = Dimer('dimer').continue_from(accurate)
    dimer_accurate = Dimer('dimer/accurate').continue_from(dimer)
    print normal.files_required
    print accurate.files_required
    print dimer.files_required
    print dimer_accurate.files_required
    normal.write()
    accurate.write()
    dimer.write()
    dimer_accurate.write()
