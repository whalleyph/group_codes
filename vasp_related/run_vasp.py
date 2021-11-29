#! /usr/bin/env python

import os
cmd = os.popen('/usr/bin/modulecmd python load vasp/5.2.12')
exec(cmd)
##os.system('which vasp')
#pwd = os.environ['PBS_O_WORKDIR']
nprocs = len(open(os.environ['PBS_NODEFILE']).readlines())
cmd = 'mpirun -np %i %s > %s' % (nprocs, 'vasp', 'vasp.out')
exitcode = os.system(cmd)
#exitcode = os.system('vasp')
