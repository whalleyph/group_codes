#!/usr/bin/env python3
import glob,os
import ase
import ase.io
from ase import Atoms
import numpy as np

from ase.io.bader import attach_charges
import matplotlib
matplotlib.use('TkAgg') #needed as OVITO uses the default QT5 backend
import matplotlib.pyplot as plt

from os import system
import argparse,time
import os.path

from ovito.io import import_file
from ovito.io.ase import *

chgs={};nstrs={}
cnt=1
for inpf in glob.iglob("./**/partial_charges_bader.xyz",recursive=1):
    #print(inpf)
    path='/'.join(inpf.split('/')[0:-1])+'/'
    #print(path)
    atoms=ase.io.read(inpf)
    #print(dir(atoms))#list of attrb and fncs
    if 1:
        ichgs=(atoms.todict()['partial_charges_bader'])
        
    else:
        attach_charges(atoms,fileobj=path+'ACF.dat')
        ichgs=atoms.get_initial_charges()
    sel=[at.index for at in atoms if at.symbol =='Li']
    #print('Selection: ',sel)
    key=len(sel)
    if key in chgs:  chgs[key]=np.append(chgs[key],ichgs[sel])
    else:chgs[key]=ichgs[sel]
    if key in nstrs: nstrs[key]+=1
    else:nstrs[key]=1
    cnt+=1

print('%d files were found.'%cnt)

keys=sorted(chgs.keys())
#mins=[];maxs=[];means=[];stds=[]
data=[]
for key in keys:
    vals=chgs[key]
    data.append([key,np.min(vals),np.max(vals),np.mean(vals),np.std(vals),np.count_nonzero(vals<0)/nstrs[key]])
    print (data[-1])

data=np.array(data).T

if 1: err=[data[3]-data[1],data[2]-data[3]] #use min,max
else: err=data[4]#use std dev.

if 1: 
    #fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True) #does not work with single row and column
    fig,ax0=plt.subplots(figsize=(11.69, 8.27));    ax1=ax0.twinx()
else: 
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(11.69, 8.27))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax0 = plt.subplot(gs[0])  #Top subplot
    ax1=ax0.twinx()

ax0.errorbar(data[0],data[3],yerr=err,fmt='-o',capsize=3,ms=3,lw=1)#ls='',marker='.'
#plt.plot(data[0],data[1],marker='-')
#plt.plot(data[0],data[2],marker='-')
ax1.plot(data[0],data[-1],c='r',marker='')
ax0.set_xlabel('#Li atoms')
ax0.set_ylabel('Charge [e]')
ax1.set_ylabel('Avg. No of Li- anions per structure')
fig.tight_layout()
plt.show()
