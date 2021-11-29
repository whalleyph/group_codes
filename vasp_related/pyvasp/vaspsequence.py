import os
import yaml
import vaspcalculation
import logging
import shutil

logging.basicConfig(filename='calculation.log',
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)


try:
    nprocs = len(open(os.environ['PBS_NODEFILE']).readlines())
except:
    nprocs = 2

def execute(fsteps):
    pwd = os.getcwd()
    steps = yaml.load(open(fsteps))
    for step in steps:
        logging.info('Running %s' % (step['name']))
        # Set up the starting and calculation directories and set up
        # the calculations
        start_from_directory = step['start_from'].strip('_')
        start_in_directory = step['start_in'].strip('_')

        # Change '.' to name of the present directory
        if start_from_directory == '.':
            start_from_directory = pwd
        if start_in_directory == '.':
            start_in_directory = pwd

        logging.info('Starting directory:  %s' % (start_from_directory))
        logging.info('Starting in directory:  %s' % (start_in_directory))

        if 'type' not in step.keys():

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Simple(start_from_directory)

            logging.info('Simple calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Simple(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'metagga':

            logging.info('Meta GGA calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Meta_GGA(start_from_directory)

            logging.info('Meta GGA calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Meta_GGA(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'hybrid':

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Simple(start_from_directory)

            logging.info('Hybrid calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Hybrid(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'optics':

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Simple(start_from_directory)

            logging.info('Optics calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Optics(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'gw':

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Optics(start_from_directory)

            logging.info('Optics calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.GW(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'band':

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Simple(start_from_directory)

            logging.info('Hybrid calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Simple_band(start_in_directory).continue_from(previous_calculation)

        elif 'type' in step.keys() and step['type'] == 'hybrid_band':

            logging.info('Simple calculation in: %s' % (start_from_directory))
            previous_calculation = vaspcalculation.Simple(start_from_directory)

            logging.info('Hybrid calculation in: %s' % (start_in_directory))
            calculation = vaspcalculation.Hybrid_band(start_in_directory).continue_from(previous_calculation)

        # Write out the files first
        calculation.write()

        # Make any changes that are required
        if 'incar' in step.keys():
            open(calculation.incar, 'w').write(step['incar'])
        if 'dincar' in step.keys():
            calculation.apply_delta_incar(step['dincar'])

        if 'kpoints' in step.keys():
            open(calculation.kpoints, 'w').write(step['kpoints'])
        if 'dkpoints' in step.keys():
            calculation.apply_delta_kpoints(step['dkpoints'])

        # clean the directory if kpoints change
        if (('dkpoints' in step.keys() or ('dincar' in step.keys() and 'kspacing' in step['dincar'].lower())) 
                and start_from_directory == start_in_directory):
            logging.info('KPOINTS changed, cleaning directory')
            os.system('vasp_realclean')

        os.chdir(calculation.directory)
        exe = step.get('exe')
        if not exe: exe = 'vasp'
        cmd = 'mpirun -np %i %s > %s' % (nprocs, exe, 'vasp.out')
        logging.info('Running %s' % (cmd))
        os.system(cmd)

        # keep any required files after the calculation
        if 'keep_files' in step.keys():
            for f in step['keep_files']:
                src = f
                dst = f + '.' + step.keys()
                shutil.copy(src, dst)

        os.chdir(pwd)

if __name__ == '__main__':
    pass
