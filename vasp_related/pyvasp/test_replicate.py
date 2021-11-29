import replicate
import sys
sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
import poscar
import numpy as np

pos = poscar.poscar().read('/home/cande/proj/Fe3C_03/Fe/POSCAR')

supercell_coord = replicate.replicate((2,2,2), pos.coord)
np.savetxt('replicate_check_1.dat', supercell_coord/2, fmt="%20.16f")

scpos = pos.replicate((2,2,2))
np.savetxt('replicate_check_2.dat', scpos.coord, fmt="%20.16f")
scpos.write('POSCAR_Fe3C_3x3x3_2')
