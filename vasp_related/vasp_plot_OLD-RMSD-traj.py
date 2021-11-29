#!/usr/bin/python

import os
from multiprocessing import Pool,Process, Value, Array,Queue, current_process, freeze_support
#import thread
#from threading import Thread
import time
#import random
from sys import exit

#
# Function run by worker processes
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        #print args
        #x=getDist(args[0],args[1])
        x=func(*args) 
        ind=int(args[1].split(".")[0].split("_")[1])
        output.put([ind,x])
        #print current_process().name, [ind,x]
    #return mat
     

def getDist(f1,f2):
    x=os.popen("dist.pl %s %s"%(f1,f2)).read()[0:-1]
    return x

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

#os.access("./TRAJ")

#ls=os.listdir(".")
if os.access("./TRAJ",os.R_OK):  print "TRAJ folder located.";os.chdir("./TRAJ")
ls=os.popen("ls config*.vasp").readlines()
if len(ls)==0: 
    print"config*.vasp could not be located. Creating the required files from OUTCAR." 
    os.system("vasp_outcar2xyz -s -d TRAJ -v"); os.chdir("./TRAJ");os.system("rm config*.xyz");ls=os.popen("ls config*.vasp").readlines()

print "%d files (time steps) were found."%(len(ls))
ls.sort()

#nsteps=int(ls[-1].split(".")[0].split("_")[1])
if os.access("./rmsd.dat",os.R_OK):  
    if raw_input("rmsd.dat located. Would you like to overwrite it? (y: yes): ")!="y": exit()

outf=open("../rmsd.dat","w")

#Serial version.
#dist_mat(ls)

#Parallel version.
NUMBER_OF_PROCESSES = 8
TASKS1 = [(getDist, (ls[0][:-1], ls[i][:-1])) for i in range(len(ls))]


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
    mat[res[0]]=res[1]


keys=mat.keys()
keys.sort()
for key in keys: outf.write("%d %s\n"%(key,mat[key]))
#print mat


# Tell child processes to stop
for i in range(NUMBER_OF_PROCESSES):
    task_queue.put('STOP')


outf.close()

print "%d processes were used. Elapsed time is %.2f sec."%(NUMBER_OF_PROCESSES, time.time()-initT)

os.system("xmgrace ../rmsd.dat&")

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
