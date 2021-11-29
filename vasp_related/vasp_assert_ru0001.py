#! /usr/bin/env python

import outcar as o
import kpoints as k
import poscar as p
import numpy as np

out = o.outcar().read('OUTCAR')
kpt = k.kpoints().read('KPOINTS')
pos = p.poscar().read('POSCAR')
p = out.params()

#print kpt
#for k, v in p.items():
#    print '%s = %s' % (k, v)

# Check INCAR
#assert p['ENAUG'] == '645.0'
#assert p['ENCUT'] == '400.0'
#assert p['EDIFF'] == '0.1E-04'
#assert p['IBRION'] == '2'
#assert p['IDIPOL'] == '0'
#assert p['ISIF'] == '2'
#assert p['ISMEAR'] == '1'
#assert p['ISPIN'] == '1'
#assert p['LDIPOL'] == 'F'
#assert p['LREAL'] == 'A'
#assert p['PREC'] == 'accura'
#assert p['SIGMA'] == '0.05'
#assert p['POTIM'] == '0.1000'


ass = """ENAUG = 645.0
ENCUT = 400.0
EDIFF = 0.1E-04
IBRION = 2
IDIPOL = 0
ISIF = 2
ISMEAR = 1
ISPIN = 1
LDIPOL = F
LREAL = A
PREC = accura
SIGMA = 0.05
POTIM = 0.1000
"""

assertions = {}

for a in ass.strip().split('\n'):
    key, value = a.split('=')
    assertions[key.strip()] = value.strip()

for key in assertions.keys():
    try:
        assert p[key] == assertions[key]
    except AssertionError:
        print key, p[key], assertions[key]


# Check KPOINTS
assert kpt.kpoints == (16,16,1)

# Check POSCAR
cell =  np.array(  [[5.4524953358405277,    0.0000000000000000,    0.0000000000000000],
                    [2.7262476679202639,    4.7219994748540612,    0.0000000000000000],
                    [0.0000000000000000,    0.0000000000000000,   25.7510009256419536]])

tolerance = 1e-6

check = (pos.cellvec - cell) < tolerance
assert check.all() == True
