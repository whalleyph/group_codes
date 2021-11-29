"""
Program to 
    - Check for the existence of any packaged tgz files
    - If they exist, extract the OUTCAR.* files
    - Check if there is an OUTCAR files
    - Find which of the OUTCAR* files is the latest
      and report it
    - If OUTCAR is the latest prompt user to ask to pack the
      latest files into a tgz file
    - Link the latest OUTCAR vasp.out to OUTCAR_FINAL and vasp.out_FINAL
"""

import os
import tarfile
from glob import glob
from numpy import argsort

def ntgzs(dir):
    full_path_tgzs = glob(os.path.join(dir, '*.tgz'))
    tgzs = [file.split(os.path.sep)[-1] for file in full_path_tgzs]
    stgzs = sorted([int(tgz.replace('.tgz', '')) for tgz in tgzs])
    stgzs = [str(i) for i in stgzs]
    return stgzs

def extract_files(dir, file):
    full_path_tgzs = glob(os.path.join(dir, '*.tgz'))
    for ntgz, tgz in zip(ntgzs(dir), full_path_tgzs):
        tar = tarfile.open(tgz)
        tar.extract(file + '.' + ntgz, path=dir)
    return None

def latest_file(dir, file):
    foutcars = glob(os.path.join(dir, file + '*'))
    if os.path.join(dir, file + '_FINAL') in foutcars:
        foutcars.remove(os.path.join(dir, file + '_FINAL'))
    mtimes = [os.stat(file)[-2] for file in foutcars]
    sorted_foutcars = argsort(mtimes, kind='mergesort')
    # print sorted_foutcars
    latest_outcar = foutcars[sorted_foutcars[-1]]
    return latest_outcar

def link_file(dir, src, dst):
    pwd = os.getcwd()
    s = src.split(os.path.sep)[-1]
    d = dst.split(os.path.sep)[-1]
    os.chdir(dir)
    try:
        os.symlink(s, d)
    except OSError:
        os.remove(src)
        os.symlink(s, d)
    os.chdir(pwd)
    os.getcwd()
    return None

def remove_files(dir, file, latestfile):
    # Remove all other OUTCARS except the last one
    files = glob(os.path.join(dir, file + '*'))
    files.remove(latestfile)
    [os.remove(file) for file in files]
    return None

def makefinal(dir, files):
    for file in files:
        extract_files(dir, file)
        latestfile = latest_file(dir, file)
        remove_files(dir, file, latestfile)
        link_file(dir, latestfile, os.path.join(dir, file + '_FINAL'))

if __name__ == '__main__':
    makefinal('.', sys.argv[1])
