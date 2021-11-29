#!/usr/bin/env python2

#Wirtten by BK

import ase.io
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--inp", type=str, default="./CONTCAR",
                    help="VASP POSCAR/CONTCAR file.")

args = parser.parse_args()


inp=args.inp

print "Reading ",inp
atoms=ase.io.read(inp,format="vasp")#,direct=True)
pos=atoms.get_scaled_positions()
atoms.set_positions(pos)
print "Writing ",inp+".nw"

atoms.write(inp+".nw",format="nwchem")
#ase.io.write(inp+".nw",pos,format='nwchem')


