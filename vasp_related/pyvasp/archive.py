import tarfile, sys
import os, glob
import shutil
import numpy as np

class archive:
  def __init__ (self, name = None):
    self.name  = name
    self.ext   = 'tgz'
    self.files = 'INCAR POSCAR OSZICAR OUTCAR CONTCAR vasp.out'.split()
    sep = '.'
    if not self.name:
      tgzlist = np.array([int(tgz.split('.')[0]) for tgz in glob.glob('*.tgz')])
      if len(tgzlist) == 0:
        self.name = str(0)
        self.tgzname = sep.join([self.name, self.ext])
      else:
        latest = max(tgzlist)
        self.name = str(latest + 1)
        self.tgzname = sep.join([self.name, self.ext])
    else:
      self.tgzname = sep.join([self.name, self.ext])

  def pack (self):
    if not os.path.exists (self.tgzname):
      tgz = tarfile.open(self.tgzname, 'w:gz')
    else:
      print "%s exists. Aborting ..." % (self.tgzname)
      sys.exit()

    for file in self.files:
      if os.path.exists (file):
        archivefile = '.'.join([file, self.name])
        shutil.copy(file, archivefile)
        tgz.add(archivefile)
        os.remove(archivefile)
      else:
        print "%s does not exist. Unable to add to archive" % (file)
    tgz.close()

  def pack_n_move (self):
    self.pack()
    shutil.move('CONTCAR', 'POSCAR')

def get_fileobjs(tgzs, filetoextract):

    """
    Returns a list of file objects of the <file> extracted from the 
    list of <tgzs>
    """

    fileobjs = []
    for tgz in tgzs:
        if tarfile.is_tarfile(tgz):
            number = tgz[:-4]
            filename = filetoextract + '.' + number
            file = tarfile.open(tgz).extractfile(filename)
            fileobjs.append(file)
        else:
            print "%s not a tar file" % (tarfile)

    return fileobjs

if __name__ == '__main__':
  archive('step0').pack_n_move()
