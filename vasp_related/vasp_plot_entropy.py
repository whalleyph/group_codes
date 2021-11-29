#!/usr/bin/python

from re import search
from os import system
from sys import exit

print "This script plots the entropy (T*S) per ion (in meV) over the trajectory given in OUTCAR. \nSee 'entropies' file for values. \n\nProper absolute values: << 1-2 meV / ion. If not then decrease the SIGMA value. \n"
inpf=open("OUTCAR","r")
outf=open("entropies","w")
outf.write("#Entropy T*S (in meV) per ion...\n")
#entropies=[]
tmp=0.0
nions=0
count=0
for ln in inpf.readlines():
	if search ("NIONS =",ln):
		nions=int(ln.split()[-1])
		print "%d ions found in OUTCAR."%nions
	elif search('entropy T\*S    EENTRO =',ln):
		if nions==0:
			print "NIONS entry could not be found in OUTCAR. Terminatiing..."
			exit()
		tmp=float(ln.split()[-1])
		#print tmp
	elif search("average \(electrostatic\) potential at core",ln):
		count += 1
		outf.write("%6d %10.6e\n"%(count,tmp/float(nions)*1000))
		tmp=0.0
		#entropies.append(tmp)

inpf.close(),outf.close()

system("xmgrace entropies")
#system("rm tmp_bk")
