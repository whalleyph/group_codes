#! /usr/bin/env python

import argparse
import re
from cStringIO import StringIO as sio
import numpy as np

def get_nbands(lines):
    for l in lines:
        if 'NBANDS' in l:
            nbands = int(l.split('=')[-1])
            break
    return nbands

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='')
    
    parser.add_argument('-i', '--infile', default='OUTCAR')
    parser.add_argument('-iv', '--ival_band', type=int, required=True)
    parser.add_argument('-ic', '--icon_band', type=int, required=True)
    
    args = parser.parse_args()
    ival_band = args.ival_band - 1
    icon_band = args.icon_band - 1

    lines = open(args.infile).readlines()
    nbands = get_nbands(lines)

    kpoint_exp = re.compile(r'^ k-point.....:')

    kpoint_data = []
    for i, l in enumerate(lines):
        if kpoint_exp.match(l):
            ikpoint_data = ''.join(lines[i+2:i+2+nbands])
            ikpoint_data = np.loadtxt(sio(ikpoint_data))
            kpoint_data.append(ikpoint_data)
            #print ikpoint_data

    kpoint_data = np.array(kpoint_data)
    val_band = kpoint_data[:,ival_band,1]
    con_band = kpoint_data[:,icon_band,1]

    #print val_band
    #print con_band

    print 'Direct band-gap:', np.min(con_band - val_band)
    print 'Indirect band-gap:', np.min(con_band) - np.max(val_band)
