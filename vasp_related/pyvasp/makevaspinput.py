#! /usr/bin/python

import sys
sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
import potcar as pot

lines = open('vaspin').readlines()
lines = [line.strip() for line in lines if line.strip() != '' and line.strip()[0] != '#']
d = {}
for line in lines:
    if line[0] == '[' and line[-1] == ']':
        file = line[1:-1]
        d[file] = []
        continue
    else:
        d[file].append(line)

for file in d.keys():
    if file == 'POTCAR':
        xc = d[file][0]
        ions = d[file][1].split()
        pot.potcar().write(ions, xc)
    else:
        f = open(file, 'w')
        f.write('\n'.join(d[file])+'\n')
        f.close()
