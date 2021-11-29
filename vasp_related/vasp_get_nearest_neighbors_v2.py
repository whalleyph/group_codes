#! /usr/bin/env python

import ase, ase.io
import numpy as np
import argparse, math
from math import sqrt
from sys import exit
from os import system

def avg(arr):
    sum=0.0
    for i in arr:
        sum += i
    return sum/float(len(arr))

def stddev(arr):
    st=0.0
    av=avg(arr)
    if len(arr)<2: return 0.0
    for i in arr:
        st+=1/float(len(arr)-1.0)*(i-av)**2
    return sqrt(st)

def get_images(position, n):
    images = []
    for i in range(-n, n+1):
        for j in range(-n, n+1):
            for k in range(-n, n+1):
                images.append([i+position[0], j+position[1], k+position[2]])
    return np.array(images)

def get_nearest_neighbors(iatom, atoms, idpairs=[],rcut=5.0 ):
    nearest_neighbors = []
    cell = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    natoms = len(atoms)
    
    #pairs=[]
    for i in range(natoms):
        pair="%d-%d"%(iatom,i)
        if pair in idpairs:continue #print "skipped ", pair; continue # for checking against double count of the given atom pairs
        images = get_images(positions[i], 1) #Should we include images in the analysis?????
        for j in images:
            rij = j - positions[iatom]
            cartesian_coordinates = np.dot(cell.transpose(), rij)
            rij = np.linalg.norm(cartesian_coordinates)
            
            if rij < rcut and i!=iatom:# and i> iatom : #no double-counting
                nearest_neighbors.append([atoms[iatom].symbol, atoms[i].symbol, rij, i])
                idpairs.append(pair)
                idpairs.append("%d-%d"%(i,iatom))

                
    #print idpairs
    #nearest_neighbors=[nearest_neighbors,idpairs]
   
    nearest_neighbors.append(idpairs)
    #print 'bora: ',nearest_neighbors
    return nearest_neighbors

def RMSD(arr1,arr2):
    #This is a function for calculating the root-mean-square deviation between two sets of values. This is a general fnc, two sets can be either Cartesian coords or internal coords or anything of the same size and no of elements.
    #Only accepts 1-D arrays.
    if len(arr1) == 0 or len(arr2) == 0 or len(arr1) != len(arr2):
        print "RMSD: Sizes of the input arrays don't match. Will return 0."
        return 0.0

    summ=0.0
    for i in range(len(arr1)):
        summ += (arr1[i]-arr2[i])**2

    RMSD=math.sqrt(summ/len(arr1))
    return RMSD        

def RMS(arr1):
    #This is a function for calculating the root-mean-square of a  set of values. 
    #Only accepts 1-D arrays.
    if len(arr1) == 0 :
        print "RMS: Length of array is zero. Will return 0."
        return 0.0

    summ=0.0
    for i in range(len(arr1)):
        summ += (arr1[i])**2

    RMS=math.sqrt(summ/len(arr1))
    return RMS

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='script to put a molecule in the center of a cubic box')
    
    parser.add_argument('-i', '--input_file', default='CONTCAR',
            help='input file in POSCAR format,default: CONTCAR')
    parser.add_argument('-n', '--atoms_list', nargs='*', type=int,
            help='list of atoms for whom the nearest neighbors should be computed (count starts at 1)')
    parser.add_argument('-nl', '--list', 
            help='atom list, e.g. 1:128')
    parser.add_argument('-rcut', '--rcut', default=4.0, type=float,
            help='cutoff radius within which the nearest neighbors should be calculated')
    parser.add_argument('-b', '--bin', default=0.01, type=float,
            help='width of the bins in the histogram.')
    parser.add_argument('-t', '--type', default='all',
            help='which type of distances to be considered in the analysis: e.g. distance of each atom in the atoms_list to  every atom inlcuded in atoms_list ("in"), or with every atom not included in atoms_list ("out"), or all atoms ("all"). \nOptions: in/out/all, default: all')
    parser.add_argument('-s', '--save', default=1,type=int,
            help='switch for saving data (*.dat) and plotting (*.png). Default value is 1 (on). To switch it off, use -s 0')


    args = parser.parse_args()
    typ= args.type
    print"Input file:   ",args.input_file
    atoms = ase.io.read(args.input_file)
    if args.atoms_list:
        atoms_list = np.array(args.atoms_list) - 1
    else:
        atoms_list = range(len(atoms))

    if args.list:
        x=args.list.split(':')
        x0=int(x[0])
        if len(x[1])==0: x1=len(atoms)
        else: x1=int(x[1])
        if x0==0 or x1>len(atoms):
             print "Wrong atom number (atom count starts at 1)";exit()  #atoms_list = range(len(atoms))
        atoms_list=range(x0-1,x1)
        print "atoms_list: ", atoms_list

    print "Type: ",typ

    ifSave = args.save
    print ifSave

    symbols=[]
    pairs=[] #pairs of atoms
    for at in atoms:
        if at.symbol not in symbols: symbols.append(at.symbol)
    for i in range(len(symbols)):
        for j in range(i,len(symbols)):
            pairs.append("%s-%s"%(symbols[i],symbols[j]))

    
    bin=args.bin
    print "Bin: ",bin

    #initiate the statistics for each pair type
    stats={}
    stats_in={}
    stats_out={}
    for p in pairs:
        stats[p]=[]
        stats_in[p]=[]
        stats_out[p]=[]

    idpairs=[] #This variable for checking against double count of the given atom pairs
    for i in atoms_list:
        #print '====== %i ======' % (i)
        gec=get_nearest_neighbors(i, atoms,idpairs, rcut=args.rcut)
        idpairs=gec.pop(-1)
        for nn in gec:
            pair="%s-%s"%(nn[0],nn[1])
            if pair not in pairs: pair="%s-%s"%(nn[1],nn[0])
            dist=float(nn[2])
            stats[pair].append(dist)
            at2=nn[-1]
            if at2 in atoms_list: stats_in[pair].append(dist)
            else: stats_out[pair].append(dist)
            

    for pair in pairs:
        if typ=="in": stat=stats_in
        elif typ=="out": stat=stats_out
        else: stat=stats
        if len(stat[pair])==0:
            #print "No entry for this pair." 
            continue
        print pair
        dists=stat[pair] #all distances for this pair type
        dists.sort()
        sts={}
        keys=[dists[0]]
        sts[keys[-1]]=[dists[0]]
        #Create histogram
        for i in range(1,len(dists)):
            key=keys[-1]
            ds=dists[i]
            if abs(ds-key)<=bin:
                sts[key].append(ds)
            else:
                keys.append(ds)
                sts[ds]=[ds]

        #Save data for each pair type to a separate file.
        if ifSave: outf=open("%s.dat"%pair,'w')
        maxc=0
        if ifSave:
            for key in keys:
                dt=sts[key]
                print "%.4f-%.4f A, count:%5d, avg: %.4f, stddev: %.4f"%(min(dt),max(dt),len(dt),avg(dt),stddev(dt))
                outf.write("%.3f-%.3f %5d\n"%(min(dt),max(dt),len(dt)))
                if len(dt)>maxc: maxc=len(dt)
            outf.close()

        print "%s: average= %.4f stddev= %.4f RMS= %.4f"%(pair,avg(dists),stddev(dists),RMS(dists))

        if ifSave:
            norm=1
            if norm:
                outf=open("%s_norm.dat"%pair,"w")
                for i in open("%s.dat"%pair).readlines():
                    x=i.split()
                    normy=float(x[1])/float(maxc)
                    outf.write("%s %.4f\n"%(x[0],normy))

                outf.close()
                system("mv %s_norm.dat %s.dat"%(pair,pair))
                maxc=1
                
            str1="""set terminal png transparent nocrop enhanced 
set output '%s.png'
set boxwidth 0.9 relative
set style data histograms
set style fill solid 1.0 border -1
set auto x
set yrange [0:%d]
set xtics rotate out scale 0,0
set xtics  norangelimit
set xtics   ()
plot '%s.dat' u 2:xtic(1) ti '%s'  """%(pair,maxc,pair,pair)
            system('echo "%s">tmp'%str1)
            system("gnuplot tmp")
            system("rm tmp ")
