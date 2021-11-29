#! /usr/bin/env python

#PBS -l nodes=1:ppn=1

exe = 'vasp'

import yaml
import os, sys
import subprocess as sp


def write_incar(s):
    f = open('INCAR', 'w')
    f.write(s)
    return None

def write_kpoints(s):
    f = open('KPOINTS', 'w')
    f.write(s)
    return None

def submit_job(cmd):
    sp.call(cmd, shell=True)
    return None

def check_success(test):
    if test.strip() == 'convergence':
        pass
    elif test.strip() == 'completion':
        pass

def make_archive():
    sp.call('vasp_move.sh', shell=True)
    return None

def move_contcar_to_poscar():
    sp.call('mv CONTCAR POSCAR', shell=True)
    return None

if __name__ == '__main__':
    np = len(os.environ['PBS_NODEFILE'].readlines())
    os.chdir(os.environ['PBS_O_WORKDIR'])
    sp.call('/usr/bin/modulecmd python load vasp/5.2.12', shell=True)
    cmd = 'mpirun -np ' + str(np) + ' ' + exe + ' > vasp.out'

    steps = yaml.load(open('steps.yaml'))

    print len(steps)

    for step in steps:
        write_incar(step['incar'])
        write_kpoints(step['kpoints'])
        #submit_job(cmd)
        #check_success(step['test'])
        #if step['keep']:
        #    make_archive()

        #if not step['keep']:
        #    move_contcar_to_poscar()
