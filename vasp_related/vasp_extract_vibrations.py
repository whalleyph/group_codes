#! /usr/bin/env python

import outcar as o
import ase.io

frequencies, positions = zip(*o.outcar().read('OUTCAR').get_vibrational_modes())
chemical_symbols = ase.io.read('POSCAR', format='vasp').get_chemical_symbols()

for i, frequency in enumerate(frequencies):
    print len(chemical_symbols)
    print frequency
    for j, position in enumerate(positions[i]):
        s = '%2s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f' % (chemical_symbols[j],
                position[0], position[1], position[2],
                position[3], position[4], position[5])
        print s
        #print chemical_symbols[j], position
