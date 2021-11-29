#!/usr/bin/env python2
from re import search,match
from sys import exit,argv
from os import system
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inp", type=str, default=["OUTCAR"],nargs="*",
                    help="Input file to read the data from. Def: ./OUTCAR")
parser.add_argument("-o","--out", type=str, #default="OUTCAR.xyz",
                    help="Output file name; a single file to concetanate trajectories from multiple input files. Def: ${inp}.xyz")
parser.add_argument("-l","--lattice", action='store_true', default=False, help="To include the lattice information in the output XYZ file. Def: No lattice info to speed-up the conversion.")

args = parser.parse_args()

initT=time.time()

if args.out: system("rm -f %s"%args.out)

for inp in args.inp:
    print "Working on %s. "%inp,

    if args.out: 
	    out=args.out
	    outf=open(out,"a")
    else: 
	    out=inp.split(".")[0]+".xyz"
	    outf=open(out,"w")

    inpf=open(inp,"r")

    atomTypes=[]
    atoms=[]
    totalCnt=0 #total no of atoms.
    atomCnt=[]
    flag=0; coords=[]; step=0; lattice=[]
    aflag=1; iflag=1;  aCnt=0

    if args.lattice:lflag=0; #Lattice info will be stored.
    else: lflag=-1

	
    #for ln in inpf.xreadlines():
    for ln in inpf:
	if aflag and search("POTCAR:",ln):
	    atomTypes.append(ln.split()[2].split('_')[0])
	elif aflag and search('local pseudopotential read in',ln):
	    aflag=0

	if not lflag and 'direct lattice vectors' in ln:
	    lattice=[];lflag=1
	elif lflag==1:
	    if 'length of vectors' in ln:
		lflag=0
		continue
	    try:
		lattice.extend(ln.split()[0:3])
	    except:
		lflag=0


	if  iflag and search("ions per type ",ln):
	    x=ln.split()
	    for i in range(4,len(x)):
		#cnt=int(x[i])
		atomCnt.append(int(x[i]))
		totalCnt+=int(x[i])

	    atomTypes.pop(-1)#the last one is repeated.
	    for i in range(len(atomTypes)):
		for j in range(atomCnt[i]):
		    atoms.append(atomTypes[i])
	    iflag=0
	    #print totalCnt

	if not flag and search("POSITION",ln): 
	    flag=1;tmp=[];step+=1
	    outf.write("%d\n"%totalCnt)
	    str2=""
	    for lt in lattice: str2+="%9.6f "%float(lt)
	    if args.lattice:        outf.write('Lattice="%s" \n'%str2[0:-1])
	    else:  outf.write('\n')
	    aCnt=0

	elif flag and search("-----------------",ln): continue
	elif flag and aCnt < totalCnt: #collecting the coords
	    outf.write("%-5s %-44s\n"%(atoms[aCnt],ln[0:44]))
	    aCnt+=1

	elif flag and search("total drift:",ln): 
	    flag=0
	    if not aCnt==totalCnt:
		print "Missing data in step %d !!"%step


    print "%d-step trajectory is written to %s."%(step,out)
    
    inpf.close();   outf.close()
    #if not args.out: outf.close()

print "Elapsed time: %.2f sec."%( time.time()-initT)  

#system("molden5.0 OUTCAR.xyz ")



###########
# DELETED #
###########
"""
    if search("POSITION",ln): flag=1;tmp=[];step+=1
    elif flag and search("-----------------",ln): continue
    elif flag: #collecting the coords
        tmp.append(ln[0:44])

    if search("total drift:",ln): 
        #print ln
        flag=0
        tmp.pop(-1)
        #coords.append(tmp)
        #print tmp
        if len(tmp)==totalCnt:
            outf.write("%d\n"%totalCnt)
	    str2=""
	    for lt in lattice: str2+="%9.6f "%float(lt)
	    outf.write('Lattice="%s" \n'%str2[0:-1])
            for i in range(totalCnt):
                outf.write("%-5s %-44s\n"%(atoms[i],tmp[i]))
        else:"Missing data in step %d !!"%step
    """
