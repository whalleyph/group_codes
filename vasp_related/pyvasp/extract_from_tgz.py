#! /bin/env python

from glob import glob
import tarfile
import sys

def extract_file(tgz, files):
    seq = tgz.split('.')[0]
    tgz = tarfile.open(tgz)
    for file in files:
        tgz.extract(file+'.'+seq)
    return None

f = sys.argv[1]
tgzs = sys.argv[2:]
print tgzs
for itgz in tgzs:
    extract_file(itgz, [f])
