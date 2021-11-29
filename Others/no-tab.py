#!/usr/bin/python

import re

outf=open('fort.4.new','w')
with open('fort.4','r') as fl:
	for ln in fl:
		#print ln
		ln=ln.replace('\t','')
		#ln=ln.rstrip()
		outf.write(ln)
outf.close()
