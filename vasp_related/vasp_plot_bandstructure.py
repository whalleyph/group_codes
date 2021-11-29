#!/usr/bin/env python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import eigenval as er
import argparse
from os import system

if (__name__ == '__main__'):
    #filename = 'EIGENVAL'
    #if len(sys.argv) == 2:
    #    filename = sys.argv[1]
    #else:
    #    print("Defaulting filename to 'EIGENVAL', if you wish to specify another filename on the command prompt do it like this: plot_eigenval.py myfilename")

    parser = argparse.ArgumentParser(description='script to create total DOSplot from a DOSCAR file')

    parser.add_argument('-i', '--input_file', default='EIGENVAL',
                        help='input file in EINGENVAL format')

    parser.add_argument('-f', '--fermi', action='store_true', default=False,
                        help='If Fermi level should be shifted to zero.')

    parser.add_argument('-y1', '--y1', default=-3.0, type=float,
            help='Start point of the y axis')
    parser.add_argument('-y2', '--y2', default=3.0, type=float,
            help='End point of the y axis')

    parser.add_argument('-s', '--shift', default=0.0, type=float,
            help='Amount of shift in the x direction')

    args = parser.parse_args()
    

    fermi=args.fermi


    parser2 = er.EigenvalParser()
    kpoints = parser2.parse(args.input_file)

    


    if fermi:
        f = open("DOSCAR")
        l = f.readlines()

        natoms = int(l[0].strip().split()[1])
        hold = l[5].strip().split()
        ndat, efermi = int(hold[2]), float(hold[3])
        print "Efermi= ", efermi


    #for i, band in enumerate(er.get_bands(kpoints)):
    #    print i, band,len(band)
    #    print
        #plot(range(0,len(band)), array(band)-efermi, 'k')
    kp=er.get_bands(kpoints)

    data=[]
    #outf=open("tmp",'w')
    for j in range(len(kp[0])):
	#str1="%d "%(j+1)
        dt=[]
	for i in range(len(kp)):
            #str1+="%.5f "%kp[i][j]
            if fermi: dt.append(kp[i][j] - efermi)  #No fermi correction needed
	    elif args.shift != 0.0: dt.append(kp[i][j] +  args.shift)
            else: dt.append(kp[i][j])
        #str1+="\n"
        data.append(dt)

        #outf.write(str1)

    y1=args.y1;y2=args.y2
    #outf.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1,len(data)+1), data)

    #xlim(0,119)
    #ylim(-25,2)
    #if fermi: plt.title("Band Structure - Efermi=%.2f eV"%efermi)
    #elif args.shift != 0: plt.title("Band Structure  Efermi=%.2f eV"%efermi)
    #else: plt.title("Band Structure")
    plt.title("Band Structure")

    xlabel="K-Points"
    if fermi: ylabel="Energy - E_f=%.2f [eV]"%efermi
    elif args.shift != 0: ylabel='E + E_shift=%.2f [eV]' % (args.shift)
    else: ylabel="Energy [eV]"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.ylim((y1,y2))
    plt.show()
    plt.savefig('bandstructure.png')

    system("display bandstructure.png &")
