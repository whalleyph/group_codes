import os
from glob import glob

def get_final_outcar(dir):
    # Check if dir exists and then try
    if os.path.exists(dir):
        try:
           if os.path.exists(os.path.join(dir, 'OUTCAR')):
                outcar = os.path.join(dir, 'OUTCAR')
           elif os.path.exists(os.path.join(dir, 'outfinal')):
               outcar = os.path.join(dir, 'outfinal')
           elif os.path.exists(os.path.join(dir, 'OUTCAR_FINAL')):
               outcar = os.path.join(dir, 'OUTCAR_FINAL')
           else:
               tgz = sorted(glob(os.path.join(dir, '?.tgz')))[-1]
               outcar = tgz
        except IndexError:
            print dir
            outcar = None

        return outcar
    else:
        print '%s does not exist' % (dir)

def get_final_oszicar(dir):
    # Check if dir exists and then try
    if os.path.exists(dir):
        try:
           if os.path.exists(os.path.join(dir, 'OSZICAR')):
                oszicar = os.path.join(dir, 'OSZICAR')
           elif os.path.exists(os.path.join(dir, 'OSZFINAL')):
               oszicar = os.path.join(dir, 'OSZFINAL')
           elif os.path.exists(os.path.join(dir, 'OSZICAR_FINAL')):
               oszicar = os.path.join(dir, 'OSZICAR_FINAL')
           else:
               tgz = sorted(glob(os.path.join(dir, '?.tgz')))[-1]
               oszicar = tgz
        except IndexError:
            print dir
            oszicar = None

        return oszicar
    else:
        print '%s does not exist' % (dir)
