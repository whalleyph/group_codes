#!/usr/bin/python

import os,argparse
from sys import exit,argv
from re import search
from copy import deepcopy 

if len(argv)==1:
    print "This script is for switching the given pairs of atoms in an input POSCAR file."

parser = argparse.ArgumentParser(description='This script is for switching the given pairs of atoms in an input POSCAR file. \nProper use:  script -l 0:1 12:25 .')
parser.add_argument('-l', '--list', nargs='*', type=str,required=True,
                    help='list of atom pairs which will be switched (Atom index starts from 0). \nProper use:  -l 0:1 12:25 .')

parser.add_argument('-f', '--file', type=str,
		    help='Name of the input file. Default is POSCAR.')

args = parser.parse_args()

if args.file:
    inpf=args.file
else:
    inpf="POSCAR"

lines=open(inpf,"r").readlines()

heading=[]
flag=0
aTypes=lines[0].split()
aCounts={}
tCount=0
cnt=0
coords={}
for i in range(len(lines)):
    ln=lines[i]
    if search('Car',ln) or search('Dir',ln):
        flag=1
        x=lines[i-1].split()
        for j in range(len(aTypes)): at=aTypes[j]; aCounts[at]=x[j];tCount+=int(aCounts[at])
        heading.append(ln)
        continue
    elif flag==1 and cnt<tCount:
        coords[cnt]=ln[0:-1].split()
        cnt+=1
    else: 
        heading.append(ln)

#print coords
if args.list:
    for ls in args.list:
        i=int(ls.split(':')[0])
        j=int(ls.split(':')[1])
	print "Switching atom ID %i with %i."%(i,j)
	print coords[i], coords[j]
        tmp=deepcopy(coords[i])
        coords[i]=deepcopy(coords[j])
        coords[j]=deepcopy(tmp)


print "Output file: %s.new"%inpf
outf=open(inpf+".new",'w')
keys=coords.keys()
keys.sort()

for h in heading: outf.write(h)
for key in keys: outf.write("%20s %20s %20s\n"%(coords[key][0],coords[key][1],coords[key][2]))

outf.close

#        x0=int(x[0])
#        if len(x[1])==0: x1=len(atoms)
#        else: x1=int(x[1])
#        if x0==0 or x1>len(atoms):
#            print "Wrong atom number (atom count starts at 1)";exit()  #atoms_list = range(len(atoms))
#            atoms_list=range(x0-1,x1)
    #print "atoms_list: ", atoms_list

