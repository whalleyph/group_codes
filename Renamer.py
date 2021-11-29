#!/usr/bin/python

from sys import exit,argv
from os import system,listdir
from re import search,sub



if len (argv)==3:
	skey=argv[1]
	repl=argv[2]
	print skey,repl
else:
	print "This script is for changing file names in a systematic way. No need to use wild cards. \nUsage: script Search_key Replace_with. e.g. script .cell .res will rename *.cell with *.ref"
	exit()

fList=[]
for fl in listdir('.'):
	#if skey in fl: fList.append(fl)  #Does not support wildcards.
	if search(r'%s'%skey,fl):fList.append(fl)
if len(fList)==0:
	print "No files affected. Terminating."; exit()
else:
	print "%d files are affected by the change. An example is: %s to %s"%(len(fList),fList[0],sub(r'%s'%skey,repl,fList[0]))#fList[0].replace(skey,repl))

if raw_input("Would you like to proceeded? (y:yes): ").lower()!="y": exit()

for fl in fList:
	print fl,sub(r'%s'%skey,repl,fl)
	system("mv %s %s"%(fl,sub(r'%s'%skey,repl,fl)))






exit()

"""
from re import search
import os
print "This script is for changing the file name group into a new one. \nFor example old.* will be changed into new.*"

os.system ("ls")
old = raw_input("Input file group name (e.g. old): ")
new = raw_input("Input new name (e.g. new): ")

liste = os.listdir('.')
expr = old
length = len(old)

fileList=[]
for i in liste:
	if search (expr,i):
		fileList.append(i)

print "List of the files being affected: "
for filen in fileList: 
    print filen
    command = "mv " + filen + " " + new + filen[length:]
#    print command
    os.system(command)
    
"""
