#! /usr/bin/env python

import os,time,re,math,argparse
from sys import exit
from copy import deepcopy

initT=time.time()

inpf=open("OUTCAR.xyz",'r')
lines=inpf.readlines()
inpf.close()

try: noAtoms=int(lines[0])
except: print "No number of atoms info in OUTCAR.xyz"; exit()
try: 
    lattice=[]
    x=lines[1].split('"')[1].split()
    #print x
    z=[]
    for j in range(len(x)):
        z.append(float(x[j]))
        if (j+1)%3==0: lattice.append(z);z=[]      

    #print lattice
except: print "No lattice info in OUTCAR.xyz"; exit()

outf=open("OUTCAR_noPBC.xyz",'w')

#coords=[] #This holds the atomic coords along the trajectory given in OUTCAR.xyz.
cstep=[]
pstep=[]
noSteps=0
#print lattice
for line in lines:
    words=line.split()
    if re.search("Lattice",words[0]):
        if noSteps!=0: pstep=deepcopy(cstep);cstep=[] #coords.append(cstep); cstep=[]
        noSteps += 1
        iatom=0
        outf.write("%d\n"%noAtoms+line)
    elif len(words)==4:
        cstep.append([float(words[1]),float(words[2]),float(words[3])])
        #if iatom==noAtoms-1:print cstep
        if noSteps==1: 
            outf.write("%-3s %8.4f %8.4f %8.4f\n"%(words[0],cstep[iatom][0],cstep[iatom][1],cstep[iatom][2]));iatom+=1 ;continue

        for i in range(3):        
            """
            if lattice[0][i]!=0: 
                delta=0.0
                if abs(cstep[iatom][i]-pstep[iatom][i])>=abs(lattice[0][i]): delta=2*lattice[0][i]
                elif abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[0][i]))<5*10**-1:delta=lattice[0][i]
                #print noSteps,iatom,"burda1",cstep[-1],lattice[0][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+delta
                else: cstep[iatom][i]=cstep[iatom][i]-delta

            if lattice[1][i]!=0: 
                delta=0.0
                if abs(cstep[iatom][i]-pstep[iatom][i])>=abs(lattice[1][i]): delta=2*lattice[1][i]
                elif abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[1][i]))<5*10**-1:delta=lattice[1][i]
                #print noSteps,iatom,"burda1",cstep[-1],lattice[0][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+delta
                else: cstep[iatom][i]=cstep[iatom][i]-delta

            if lattice[2][i]!=0: 
                delta=0.0
                if abs(cstep[iatom][i]-pstep[iatom][i])>=abs(lattice[2][i]): delta=2*lattice[2][i]
                elif abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[2][i]))<5*10**-1:delta=lattice[2][i]
                #print noSteps,iatom,"burda1",cstep[-1],lattice[0][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+delta
                else: cstep[iatom][i]=cstep[iatom][i]-delta

            """
            if lattice[0][i]!=0 and abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[0][i]))<5*10**-1: #cstep[iatom][i]-=lattice[2][i]
                #print noSteps,iatom,"burda1",cstep[-1],lattice[0][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+lattice[0][i]
                else: cstep[iatom][i]=cstep[iatom][i]-lattice[0][i]

            if lattice[1][i]!=0 and abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[1][i]))<5*10**-1:
                #print noSteps,iatom,"burda2",cstep[-1],lattice[1][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+lattice[1][i]
                else: cstep[iatom][i]=cstep[iatom][i]-lattice[1][i]

            if lattice[2][i]!=0 and abs(abs(cstep[iatom][i]-pstep[iatom][i])-abs(lattice[2][i]))<5*10**-1: #cstep[iatom][i]-=lattice[2][i]
                #print noSteps,iatom,"burda3",cstep[-1],lattice[2][i]
                if cstep[iatom][i]<=pstep[iatom][i]: cstep[iatom][i]=cstep[iatom][i]+lattice[2][i]
                else: cstep[iatom][i]=cstep[iatom][i]-lattice[2][i]
                
            #else: None #no change in coord.
        #print cstep[-1]
        outf.write("%-3s %8.4f %8.4f %8.4f\n"%(words[0],cstep[iatom][0],cstep[iatom][1],cstep[iatom][2]))
        #if iatom==noAtoms-1: print cstep
        iatom+=1
        
    elif len(words)<=1: continue
    #if line==lines[-1]: pstep=deepcopy(cstep)#coords.append(cstep);

#print len(coords),noSteps
#if noSteps!=len(coords): print 'Problem with the number of steps read  in OUTCAR.xyz.'; exit()
#else: print noSteps

inpf.close();outf.close()
