#!/usr/bin/env python
import argparse
from os import system
from sys import exit


parser = argparse.ArgumentParser(
            description='For getting the interpolation between two given vectors.')
    
parser.add_argument('-v1', '--v1', type=str,default='[0.0,0.0,0.0]')
parser.add_argument('-v2', '--v2', type=str,default='[0.1,0.1,0.1]')
parser.add_argument('-np', '--nPts',type=int,default=25)
parser.add_argument('-o', '--out',type=str,default="lin_int_out")

args=parser.parse_args()
np=args.nPts-1
v1=[float(x) for x in args.v1.replace("[","").replace("]","").split(",")]
v2=[float(x) for x in args.v2.replace("[","").replace("]","").split(",")]

if len(v1)!=len(v2): print "Sizes of two vectors do not match. Terminating...";exit()


st=[]
for i in range(len(v1)):
    st.append((v2[i]-v1[i])/float(np))

ipV=[v1]
for i in range(np):
    tmp=[]
    for j in range(len(v1)):
        tmp.append(v1[j]+(st[j]*(i+1)))
    ipV.append(tmp)
#ipV.append(v2)
outf=open(args.out,"w")
outf.write("!Line\n%d\nReciprocal\n"%(np+1))
for i in ipV: 
    #st=""
    for s in i: outf.write("%.5f "%s)
    outf.write(" 1\n")
outf.write("\n")
outf.close()
