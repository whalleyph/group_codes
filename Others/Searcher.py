#!/usr/bin/python

import os,sys

if len(sys.argv)==1:print """
This script is for searching for a file in the current directory and in
the given level of subdirectories of it. You can also use wild-card characters
(i.e. *,? and etc) within "" for advanced search.
Search results are given in 'search-result' file.

Proper use: script "SearchKey" level
"""

os.system ('rm search_result')
os.system ('rm search_result_sizes')


if  os.system("ls */* >gec"):
    print "No subdirectories, listing only the matching files in the current directory."
    for i in range (1,len(sys.argv)):
        try: int(sys.argv[i])
        except:
            arg=sys.argv[i]
            os.system("ls "+arg+" >> search_result")
            os.system("du -ksh "+arg)
    os.system("rm -f gec")
    sys.exit()
os.system("rm -f gec")

if len(sys.argv) == 1:
    key = raw_input("Filename to find : ")
    level = raw_input('How many levels to go through the subdirectories (def 5): ')
elif len(sys.argv) == 2:
    key = sys.argv[1]
    print key
    level = raw_input('How many levels to go through the subdirectories (def 5): ') 
elif len(sys.argv) >= 3:
    key = sys.argv[1]
    print key
    level = int(sys.argv[2])

if level=="": level=5
else: level=int(level)

if level == 1:
    prefix = ''
else:
    prefix = "."

jobs=os.popen("myjobdirs").readlines()
jobdirs=[]
for j in jobs: jobdirs.append(j[0:-1])
print jobdirs
outf=open("search_result",'w')
pwd=os.popen('echo `pwd`','r').read()[0:-1]

if level>=1:
    for i in range (1,level+1):
        if i == 1:
            searchkey = key
        else:
            searchkey = prefix+"/"+key

        print "\nAt "+ str(i) + ". subdirectory level. Searching for: "+ searchkey
        #os.system("find "+searchkey+" >> search_result")
	dirs=os.popen("find "+searchkey).readlines()
	for dir in dirs:
		dir=pwd+dir[1:]
		dirr="/".join(dir.split("/")[0:-1])
		#print dirr
		#for jd in jobdirs: 
			#jd=jd[0:-1]
		if dirr not in  jobdirs: outf.write(dir)
		#else: print dir
        os.system("du -ksh "+searchkey+" >> search_result_sizes")

        prefix += "/*"

outf.close()
		
