#!/usr/bin/python

from sys import argv, exit


if len(argv)!=2:
    print "This script is for reordering the atom types in a given xyz file. Atoms are regrouped so that atom types are sorted ascending and the atom order within the same atom type group is not changed. Not compatible with multiple structure in the same file. "
    print 'Proper use: script xyzfile.'
    exit()


try: inpf=open(argv[1],"r")
except: print '%s could not be found in the current directory.'%argv[1];exit()

atoms={}

cnt=0
noAtoms=0
aT=''

for ln in inpf.readlines():
    cnt += 1
    if cnt==1:
        noAtoms=int(ln[0:-1])
    
    elif cnt>=3:
        ws=ln.split()
        if len(ws)==4:
            aT=ws[0]
            if atoms.has_key(aT):
                atoms[aT].append(ln)
            else:
                atoms[aT]=[ln]
	elif len(ws)>4:
            aT=ws[0]
	    spl=ln.split()[0:4]
	    ln = "%-3s %s %s %s\n"%(spl[0],spl[1],spl[2],spl[3])
            if atoms.has_key(aT):
                atoms[aT].append(ln)
            else:
                atoms[aT]=[ln]
 
        else:
            #second structure in xyz file
            break
    else: continue

keys=atoms.keys()
keys.sort()


outf=open('sorted.xyz',"w")
outf.write("%d \n\n"%noAtoms)
for key in keys:
    for ln in atoms[key]:
        outf.write(ln)

inpf.close();outf.close()
    
print "New atom order is stored in 'sorted.xyz'."               
