#!/usr/bin/env python

import ase.io
import argparse,os.path
from ase import Atoms
from re import search


parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                    help='an integer for the accumulator')


parser.add_argument('-i','--inpf',  type=str,required=True, help='Input file')

parser.add_argument('-o','--outf', type=str,required=False, default="PRIM",help='Output file name. Def: PRIM.')

#parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')

args = parser.parse_args()


inpf=open(args.inpf,'r')
outf=open(args.outf,'w')

flag=0;cnt=0;aT=[];aTypes=[];aCounts=[]
aInd=0
for ln in inpf:
    cnt += 1
    #print ln.lower()
    #print aTypes
    if cnt==1: 
        aT=ln.split() #doesn't have to be on the first line.

    elif cnt==6 or cnt==7: #VASP 4 or 5 format.
        words=ln.split()
        try:
            for i in range(len(words)):
                w=words[i]
                aCounts.append(int(w))
                for j in range(aCounts[-1]):
                    aTypes.append(aT[i])
                
        except:
            outf.writelines(ln);#continue
            if search("cart",ln.lower()) or search("dir" ,ln.lower()): 
                flag=1
                continue
    elif search("cart",ln.lower()) or search("dir" ,ln.lower()): flag=1
    #elif ln.lower()[0:2]=="car" or ln.lower()[0:2]=="dir":     flag=1

    elif flag:
        #atom coords
        ln= ln[0:-1]+ "%4s\n"%aTypes[aInd]
        aInd += 1


    outf.writelines(ln)

outf.close();inpf.close()
