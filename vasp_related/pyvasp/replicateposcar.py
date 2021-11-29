import numpy as np
from cStringIO import StringIO as sio

def replicate(dc, supercell):
    nx, ny, nz = supercell
    # Supercell direct coordinates
    sdc = []
    for x, y, z in dc:
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    sdc.append([(x+i)/nx, (y+j)/ny, (z+k)/nz])

    sdc = np.array(sdc)

    return sdc


if __name__ == '__main__':
    import sys

    cmdargs = sys.argv
    # print cmdargs

    supercell = [int(i) for i in cmdargs[2].split()]
    nx, ny, nz = supercell
    fposcar = cmdargs[1]

    # Total number of cells
    ncells = nx*ny*nz
    
    # Total number of ions
    nionsline = open(fposcar).readlines()[5]
    nions = sum([int(i) for i in nionsline.split()])
    # print nions

    # Direct coordinates
    dc = np.loadtxt(fposcar, skiprows=7, usecols=(0,1,2))
    dc = dc[:nions]

    sdc = replicate(dc, supercell)

    # Make POSCAR for supercell
    poscarlines = open(fposcar).readlines()
    comment = poscarlines[0].strip()
    scale = poscarlines[1].strip()
    cellvec = np.loadtxt(sio(''.join(poscarlines[2:5])))
    cellvec = cellvec*supercell
    s = ''
    for i in cellvec:
        s += "%20.16f %20.16f %20.16f\n" % (i[0], i[1], i[2])
    cellvec = s.strip()
    supercellnions = [int(i)*ncells for i in nionsline.split()]
    supercellnions = ' '.join([str(i) for i in supercellnions])
    direct = 'Direct'

    # String for supercell coordinates
    ssdc = ''
    for i in sdc:
        ssdc += "%18.16f %18.16f %18.16f\n" % (i[0], i[1], i[2])
    ssdc = ssdc.strip()
    
    supercellposcar = '\n'.join([
        comment,
        scale,
        cellvec,
        supercellnions,
        direct,
        ssdc])

    print supercellposcar

                            
