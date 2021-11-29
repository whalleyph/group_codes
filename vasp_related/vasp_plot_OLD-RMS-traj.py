#! /usr/bin/env python

###!/usr/bin/python
###!/home/cande/opt/Enthought/Canopy_64/User/bin/python

import os,time
from multiprocessing import Pool,Process, Value, Array,Queue, current_process, freeze_support
#import thread
#from threading import Thread
import argparse
#import random
from sys import exit
from copy import deepcopy

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        #print args
        #x=getDist(args[0],args[1])
        x=func(*args) 
        #ind=int(args[0].split(".")[0].split("_")[1])
        #output.put([ind,x[0:]])
        output.put(x)
        #print current_process().name, [ind,x]
    #return mat
     

def getDist(f1,f2):
    x=os.popen("dist.pl %s %s"%(f1,atomID)).read()[0:-1]
    return x

def getRMS(f1,atoms_list,rcut,bin):
    y="vasp_get_nearest_neighbors_v2.py -i %s  -n %s -s 0 -rcut=%-.2f -b=%-.2f | grep 'RMS=' "%(f1,atoms_list,rcut,bin)
    #print y
    x=os.popen(y).readlines()
    #data=""
    ind=int(f1.split(".")[0].split("_")[1])
    data=[ind]
    for i in x: 
        wo=i.split(); #data += wo[0][:-1]+" "+ wo[-1] + " " #take the RMS data
        data.append(wo[0][:-1])
        data.append(float(wo[-1]))
        #try:data.append(float(i))
        #except:data.append(i)
        
    #print data
    return data

def dist_mat(inp): #Not needed anymore (serial version).
    mat=[]
    for i in range(1,len(inp)): 
        #print i
        res=getDist(inp[0][:-1],inp[i][:-1])    
        mat.append([i+1,res])
    return mat

#############
#Main script#
#############
initT=time.time()


parser = argparse.ArgumentParser(description='Computes the RMS of different bond types for a given list of atoms along a MD trajectory.')
parser.add_argument('-n', '--atoms_list', nargs='*', type=int,#required=True
                    help='list of atoms for whom the nearest neighbors should be computed (count starts at 1)')
parser.add_argument('-nl', '--list',
                    help='atom list, e.g. 1:128')
parser.add_argument('-rcut', '--rcut', default=10.0, type=float,
                    help='cutoff radius within which the nearest neighbors should be calculated')
parser.add_argument('-b', '--bin', default=0.01, type=float,
                    help='width of the bins in the histogram.')

args = parser.parse_args()

if args.atoms_list:
    atoms_list=args.atoms_list
    #atoms_list = np.array(args.atoms_list) - 1
    #for i in range(len(atoms_list)): atoms_list[i] = int(atoms_list[i])-1
elif args.list:
    x=args.list.split(':')
    x0=int(x[0])
    if len(x[1])==0: x1=len(atoms)
    else: x1=int(x[1])
    if x0==0 or x1>len(atoms):
        print "Wrong atom number (atom count starts at 1)";exit()  #atoms_list = range(len(atoms))
    atoms_list=range(x0,x1+1)
    #print "atoms_list: ", atoms_list
else:
    #atoms_list = range(len(atoms))
    print "The list of atoms for which the RMS of bond lengths will be computed must be given. Use -n or -nl for this purpose. See help (-h) for furhter information."
    exit()

aList=""
for i in args.atoms_list:
    aList+=str(i)+" "

bin=args.bin
#print "Bin: ",bin
rcut=args.rcut


#ls=os.listdir(".")
if os.access("./TRAJ",os.R_OK):  print "TRAJ folder located.";os.chdir("./TRAJ")
ls=os.popen("ls config*.vasp").readlines() 
if len(ls)==0: ls=os.popen("ls config*.POSCAR").readlines()
if len(ls)==0: 
    print"config*.POSCAR could not be located. Creating the required files from OUTCAR." 
    #os.system("vasp_outcar2xyz -s -d TRAJ -v"); os.chdir("./TRAJ");os.system("rm config*.xyz");ls=os.popen("ls config*.POSCAR").readlines()
    os.system("vasp_outcar2xyz -s -d TRAJ -v"); os.chdir("./TRAJ");ls=os.popen("ls config*.POSCAR").readlines()

print "%d files (time steps) were found."%(len(ls))

#if ls[0].split(".")[-1][0:4]=="vasp":
#    print "Renaming config*.vasp to config*.POSCAR..."
#    for i in ls: core=i.split(".")[0]; os.system("mv %s.vasp %s.POSCAR"%(core,core))
#    ls=os.popen("ls config*.POSCAR").readlines()

ls.sort()

#nsteps=int(ls[-1].split(".")[0].split("_")[1])
if os.access("./rms-bonds.dat",os.R_OK):  
    if raw_input("rmsd.dat located. Would you like to overwrite it? (y: yes): ")!="y": exit()


#Serial version.
#dist_mat(ls)

#Parallel version.
NUMBER_OF_PROCESSES = 8
#TASKS1 = [(getDist, (ls[0][:-1], ls[i][:-1])) for i in range(len(ls))]
TASKS1 = [(getRMS, (ls[i][:-1],aList,rcut,bin)) for i in range(len(ls))]


freeze_support()
# Create queues
task_queue = Queue()
done_queue = Queue()

# Submit tasks
for task in TASKS1:
    task_queue.put(task)

    # Start worker processes
for i in range(NUMBER_OF_PROCESSES):
    Process(target=worker, args=(task_queue, done_queue)).start()
    #print b.mat


# Get and print results
mat={}
for i in range(len(TASKS1)): #get unordered results from each process.
    #print '\t', done_queue.get()
    res=done_queue.get()
    mat[res[0]]=res[1:]


outf=open("../rms-bonds.dat","w")
keys=mat.keys()
#print keys
keys.sort()
ivals=[]
for key in keys: 
    vals=mat[key]
    str1=""
    if key == 0: 
        ivals=deepcopy(vals);
        print ivals
        str1="#Step "
        for i in range(0,len(vals),2):
            #mat[0][i+1]=0.0
            str1+="%s "%mat[0][i]
        outf.write(str1+"\n")
    else:
        str1="%d "%key
        for i in range(1,len(vals),2): #Determine the relative difference wrt to the first step.
            #if i%2==0:
            #print mat[key][i]
            mat[key][i]=abs(mat[key][i]-ivals[i])/ivals[i]
            str1 += "%.4f "%(mat[key][i])
        outf.write(str1+"\n")
#print mat


# Tell child processes to stop
for i in range(NUMBER_OF_PROCESSES):
    task_queue.put('STOP')


outf.close()

print "%d processes were used. Elapsed time is %.2f sec."%(NUMBER_OF_PROCESSES, time.time()-initT)

#os.system("xmgrace ../rms-bonds.dat&")

#if raw_input("Would you like to keep the trajectory files? (y: yes): ")!="y":
os.system("rm -rf ../TRAJ")



########
#BACKUP#
########

#These do not use the processes as planned (0% usage).

#p = Process(target=dist_mat, args=([ls]))
#p.start()
#p.join()

#p = Pool(4)
#rs=p.apply_async(dist_mat, [ls])
#TASKS=[]
#for i in range(1,100):
#    TASKS.append([getDist,ls[0],ls[i]])
#rs=p.apply_async(TASKS) 
#rs=p.imap(dist_mat,TASKS)
#print rs
#print rs.get()
#print p.map(dist_mat,[ls],1)
