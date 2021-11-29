#! /usr/bin/env python

import sys

print 'Enter Nose frequency in Hz'

nf2 = float(sys.argv[1])

foutcar = 'OUTCAR'

lines = open(foutcar).readlines()

for line in lines:
    if 'Nose-frequenzy' in line:
        nf1 = float(line.strip().split()[4])
        break

print 'Existing Nose frequency: %.2e' % (nf1)

for line in lines:
    if 'SMASS' in line:
        sm1 = float(line.strip().split()[2])
        break

sm2 = sm1*(nf1/nf2)**2
print 'Suggested SMASS: %.2f' % (sm2)
