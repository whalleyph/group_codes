#!/usr/bin/python 

import sys,re,os

def next_file_name():
    num = 1
    while True:
        file_name = '%06d.cfg' % num
        if not os.path.exists(file_name):
            return file_name
        num += 1
        if num==1000: break

def main():

	#open a directory
	if not os.path.exists("./Cfg"): 
	  os.makedirs("./Cfg")
	

	#open XDATCAR file
	if not os.path.exists("./XDATCAR"): 
	  print "Oh dear, there is no XDATCAR"
	  sys.exit(1)
	f = open("XDATCAR",'r') 

	# name of calculation 
	name = f.readline()
	print name,

	# scale factor 
	scaleFactor = int(f.readline())
	print '%d' %scaleFactor

	# simulation box
	axes = f.readline()
	X = re.findall('.\d+\.\d+', axes)
	print X[0],X[1],X[2]
	axes = f.readline()
	Y = re.findall('.\d+\.\d+', axes)
	print Y[0],Y[1],Y[2]
	axes = f.readline()
	Z = re.findall('.\d+\.\d+', axes)
	print Z[0],Z[1],Z[2]

	# name of elements in the calculation 
	elements = f.readline()
	atoms = re.findall(r'\w+', elements)
	print atoms

	# number of elements in the calculation 
	numbers = f.readline()
	counters = re.findall(r'\d+',numbers)
	print counters
	#convert string to int
	T2 = [int (x) for x in counters]

	#change the directory
	os.chdir("./Cfg")

	# read XDATCAR file and write *.cfg files
	while True:
	  nothing = f.readline()
	  if nothing=='': break

	  fo = open((next_file_name()),'wb')
	  fo.write('Number of particles = %d\n\n' %sum(T2))
	  fo.write('H0(1,1) = %s\n' %X[0])
	  fo.write('H0(1,2) = %s\n' %X[1])
	  fo.write('H0(1,3) = %s\n' %X[2])
	  fo.write('H0(2,1) = %s\n' %Y[0])
	  fo.write('H0(2,2) = %s\n' %Y[1])
	  fo.write('H0(2,3) = %s\n' %Y[2])
	  fo.write('H0(3,1) = %s\n' %Z[0])
	  fo.write('H0(3,2) = %s\n' %Z[1])
	  fo.write('H0(3,3) = %s\n\n' %Z[2])
	  
	  for i in zip(atoms,T2):
	    for j in range(i[1]):
	      coord = f.readline()
	      T1 = re.findall('.\d+\.\d+', coord)
	      fo.write(' 1.00\t')
	      fo.write(i[0])
	      fo.write('\t%s\t' %T1[0])
	      fo.write('%s\t' %T1[1])
	      fo.write('%s\t' %T1[2])
	      fo.write('0.00\t0.00\t0.00\n')
	  fo.close()

	# change the directory
	os.chdir("..")

	# close the XDATCAR
	f.close()

if __name__ == '__main__':
	main()
