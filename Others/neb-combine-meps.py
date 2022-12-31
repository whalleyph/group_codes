#!/usr/bin/env python3

from os import system
import argparse


#############
#Main script#
#############
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Computes the RMS of different bond types for a given list of atoms along a MD trajectory.')

    #General options
    parser.add_argument('-v','--verb',default=False,action='store_true',help='Verbose output.')

    parser.add_argument('-i', '--inpf', type=str,required=1,nargs='*',
                        help='Input neb.dat files. ')

    parser.add_argument('-o', '--outf', type=str,default="neb_new.dat",
            help='Output file (def: neb_new.dat).')

    parser.add_argument('-e0', '--Eref',type=float, default=0.0,help="Energy of the first step (E_ref). Def: 0.0 eV")

    args = parser.parse_args()


    outf=open(args.outf,'w')
    cnt=0 #step cnt
    Sref=0. #reference coord
    Eref=args.Eref
    for i,inpf in enumerate(args.inpf):
        with open(inpf, "r") as f:
            lines=f.readlines()
            lcnt=1
            for ln in lines:
                if i>0 and lcnt==1: lcnt+=1; continue
                w=[float(l) for l in ln.split()] #stepID, coord, energy, force, stepID
                outf.write('%3d%13.6f%13.6f%13.6f%4d\n'%(cnt,Sref+w[1],Eref+w[2],w[3],cnt))
                cnt+=1
                lcnt+=1

            Sref=w[1]
            Eref=args.Eref+w[2]

    outf.close()


    #system('nebspline.pl') # which reads the neb.dat to create mep.jpg
