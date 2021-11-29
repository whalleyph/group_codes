#! /usr/bin/env python

import os, sys
from StringIO import StringIO
#import myvasp as v
from copy import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import outcar

#import thesis

#thesis.double()

def nelec(eng, dos):
    l = []
    for i in range(len(eng)):
        l.append(np.trapz(dos[:i+1], eng[:i+1]))
    return np.array(l)

class doscar:

    def __init__ (self):

        self._file_ = 'DOSCAR'
        self.lines = ''


    def read(self, file):

        if file is not None:
            try:
                file.pos
            except AttributeError:
                self._file_ = open(file)
            else:
                self._file_ = file
                self._file_.seek(0)
        else:
            self._file_ = open(self._file_)

        self.lines = self._file_.readlines()

        return self


    def get_total_number_of_ions(self):
        out = outcar.outcar().read('OUTCAR')
        return out.sumnions()


    def get_dos (self, ispin=2, lorbit=11):

        dos = {}

        if ispin==2 and lorbit==11:

            ndata = int(self.lines[5].split()[2])
            dos['efermi'] = float(self.lines[5].split()[3])
            data = StringIO(''.join(self.lines[6:6+ndata]))
            dos['eng'], dos['up'], dos['down'], dos['nelup'], dos['neldown'] = \
                np.loadtxt(data, dtype=float, unpack = True)

            return dos


        elif ispin==1 and lorbit==11:

            ndata = int(self.lines[5].split()[2])

            dos['efermi'] = float(self.lines[5].split()[3])

            data = StringIO(''.join(self.lines[6:6+ndata]))

            dos['eng'], dos['total'], dos['nel'] = \
                np.loadtxt(data, dtype=float, unpack = True)

            return dos



    def idcount(self, ids):
        idcount = {}
        uniqids = []
        for id in ids:
            if id not in uniqids:
                uniqids.append(id)
                idcount[id] = 1
            else:
                idcount[id] += 1

        return idcount



    def add_dos_spin(self, data):
        """\
        Take: Data frame from DOSCAR.

        Return:
        lm decomposed, l decomposed and spin decomposed partial
        density of states

        Data frame format:
        eng sup sdown p(-1)up p(-1)down p(0)up p(0)down p(1)up p(1)down
        0   1   2     3       4         5      6        7      8
        d(-2)up d(-2)down d(-1)up d(-1)down d(0)up d(0)down
        9       10        11      12        13     14
        d(1)up d(1)down d(2)up d(2)down
        15     16       17     18

        >>> ids = 4*['C'] + 8*['Fe8d'] + 4*['Fe4c']
        >>> pdos = doscar().read('test_files/DOSCAR_spin').get_pdos(ids=ids)
        >>> iion = 15

        >>> print pdos[iion]['eng'][18], pdos[iion]['s'][18], 's'
        -6.095 0.03658 s

        >>> print pdos[iion]['eng'][18], pdos[iion]['p'][18], 'p'
        -6.095 0.0032556 p

        >>> print pdos[iion]['eng'][18], pdos[iion]['d'][18], 'd'
        -6.095 0.012130257 d

        >>> print pdos[iion]['eng'][18], pdos[iion]['su'][18], 'su'
        -6.095 0.01839 su

        >>> print pdos[iion]['eng'][18], pdos[iion]['sd'][18], 'sd'
        -6.095 0.01819 sd

        >>> print pdos[iion]['eng'][18], pdos[iion]['pu'][18], 'pu'
        -6.095 0.0016343 pu

        >>> print pdos[iion]['eng'][18], pdos[iion]['pd'][18], 'pd'
        -6.095 0.0016213 pd

        >>> print pdos[iion]['eng'][18], pdos[iion]['du'][18], 'du'
        -6.095 0.006408405 du

        >>> print pdos[iion]['eng'][18], pdos[iion]['dd'][18], 'dd'
        -6.095 0.005721852 dd

        >>> print pdos[iion]['eng'][18], pdos[iion]['up'][18], 'up'
        -6.095 0.026432705 up

        >>> print pdos[iion]['eng'][18], pdos[iion]['down'][18], 'down'
        -6.095 0.025533152 down

        >>> print pdos[iion]['eng'][18], pdos[iion]['total'][18], 'total'
        -6.095 0.051965857 total

        >>> print pdos['efermi'], 'efermi'
        7.2078783 efermi

        >>> print pdos['Fe4c']['eng'][18], pdos['Fe4c']['d'][18], 'Fe4c-d'
        -6.095 0.01113346225 Fe4c-d

        >>> print pdos['Fe8d']['eng'][18], pdos['Fe8d']['d'][18], 'Fe8d-d'
        -6.095 0.0064093125 Fe8d-d

        >>> print pdos['C']['eng'][18], pdos['C']['p'][18], 'C-p'
        -6.095 0.0005635925 C-p
        """

        cols = 'eng \
                s0u s0d \
                pm1u pm1d p0u p0d pp1u pp1d \
                dm2u dm2d dm1u dm1d d0u d0d dp1u dp1d dp2u dp2d'.split()
        
        # Change all the data into a dictionary of pdos

        pdos = {}

        # lm decomposed DoS
        for icolumn, column in enumerate(cols):
            pdos[column] = data[:,icolumn]


        # l decomposed DoS

        pdos['s'] = np.sum(data[:,1:3], axis=1)
        pdos['p'] = np.sum(data[:,3:9], axis=1)
        pdos['d'] = np.sum(data[:,9:19], axis=1)

        # l and spin decomposed DoS

        pdos['su'] = data[:,1]
        pdos['sd'] = data[:,2]


        pdos['pu'] = data[:,3] + data[:,5] + data[:,7]
        pdos['pd'] = data[:,4] + data[:,6] + data[:,8]

        pdos['du'] = data[:,9] + data[:,11] + data[:,13] + \
                     data[:,15] + data[:,17]
        pdos['dd'] = data[:,10] + data[:,12] + data[:,14] + \
                     data[:,16] + data[:,18]

        # Spin decomposed DoS

        pdos['up'] = data[:,1] + data[:,3] + data[:,5] + \
                    data[:,7] + data[:,9] + data[:,11] + \
                    data[:,13] + data[:,15] + data[:,17]
        pdos['down'] = data[:,2] + data[:,4] + data[:,6] + \
                    data[:,8] + data[:,10] + data[:,12] + \
                    data[:,14] + data[:,16] + data[:,18]

        # Total DoS

        pdos['total'] = np.sum(data[:,1:], axis=1)


        return pdos



    def add_dos_nospin(self, data):
        cols = 'eng \
                s0 \
                pm1 p0 p1 \
                dm2 dm1 d0 dp1 dp2'.split()
        
        # Change all the data into a dictionary of pdos

        pdos = {}

        # lm decomposed DoS
        for icolumn, column in enumerate(cols):
            pdos[column] = data[:,icolumn]


        # l decomposed DoS

        pdos['s'] = np.sum(data[:,1:2], axis=1)
        pdos['p'] = np.sum(data[:,2:5], axis=1)
        pdos['d'] = np.sum(data[:,5:10], axis=1)

        # Total DoS

        pdos['total'] = np.sum(data[:,1:], axis=1)

        return pdos



    def get_pdos (self, ids=None, ispin=2, lorbit=11):

        """\
        Valid for ISPIN = 2 and LORBIT = 11
        Assuming the format
        Eng s_up s_down p(-1)_up p(-1)_down p(0)_up ... d(2)_up d(2)_down
        """

        # Skip the total dos part
        if ispin == 2 and lorbit == 11:
            pdos = {}

            pdos['efermi'] = float(self.lines[5].split()[3])
            nions = int(self.lines[0].split()[0])
            ndata = int(self.lines[5].split()[2])
            nskip = 6 + ndata
            lines = self.lines[nskip:]

            if ids and len(ids) != nions:
                print "Something wrong in the ids??", nions, len(ids), ids
                sys.exit(1)

            for i in range(nions):
                pdos[i] = {}
                data_chunk = StringIO(''.join(lines[:ndata+1]))
                #print data.getvalue()
                lines = lines[ndata+1:]
                data = np.loadtxt(data_chunk, dtype=float, skiprows=1)
                data_chunk.close()

                pdos[i] = self.add_dos_spin(data)
                if ids:
                    id = ids[i]
                    if not pdos.has_key(id):
                        pdos[id] = copy(pdos[i])
                    else:
                        for key in pdos[id].keys():
                            pdos[id][key] = pdos[id][key] + pdos[i][key]

            # Average DoS for all unique ids
            if ids:
                for id, count in self.idcount(ids).items():
                    for key in pdos[id].keys():
                        pdos[id][key] = pdos[id][key]/count

            # Sum of all total partial DoS

            pdos['eng'] = pdos[0]['eng']
            pdos['total'] = pdos[0]['total']
            pdos['up'] = pdos[0]['up']
            pdos['down'] = pdos[0]['down']
            for i in range(1,nions):
                pdos['total'] = pdos['total'] + pdos[i]['total']
                pdos['up'] = pdos['up'] + pdos[i]['up']
                pdos['down'] = pdos['down'] + pdos[i]['down']
            return pdos 

        elif ispin == 1 and lorbit == 11:
            pdos = {}

            pdos['efermi'] = float(self.lines[5].split()[3])
            nions = int(self.lines[0].split()[0])
            ndata = int(self.lines[5].split()[2])
            nskip = 6 + ndata
            lines = self.lines[nskip:]

            if ids and len(ids) != nions:
                print "Something wrong in the ids??", nions, len(ids), ids
                sys.exit(1)

            for i in range(nions):
                pdos[i] = {}
                data_chunk = StringIO(''.join(lines[:ndata+1]))
                #print data.getvalue()
                lines = lines[ndata+1:]
                data = np.loadtxt(data_chunk, dtype=float, skiprows=1)
                data_chunk.close()

                pdos[i] = self.add_dos_nospin(data)
                if ids:
                    id = ids[i]
                    if not pdos.has_key(id):
                        pdos[id] = copy(pdos[i])
                    else:
                        for key in pdos[id].keys():
                            pdos[id][key] = pdos[id][key] + pdos[i][key]

            # Average DoS for all unique ids
            if ids:
                for id, count in self.idcount(ids).items():
                    for key in pdos[id].keys():
                        pdos[id][key] = pdos[id][key]/count

            # Sum of all total partial DoS

            pdos['eng'] = pdos[0]['eng']
            pdos['total'] = pdos[0]['total']
            for i in range(1,nions):
                pdos['total'] = pdos['total'] + pdos[i]['total']
            return pdos 



    def plot_dos (self, ispin,lorbit,sfermi=True, show=False, axis_limits=None,shift=0.0):
        import numpy as np
        import pylab as plt
        dos = self.get_dos()
        if sfermi:
            energy = dos['eng'] - dos['efermi']
            xlabel = r'$E - E_{f}\ \ \left[eV\right]$'
        elif shift!=0.0:
            energy = dos['eng'] + shift
            xlabel = r'$E + %.2f\ \ \left[eV\right]$'%shift
        else:
            energy = dos['eng']
            xlabel = r'$E\ \ \left[eV\right]$'

        magnetization = dos['nelup'] - dos['neldown']

        # Instantiate the figure
        fig = plt.figure()
        plot = fig.add_subplot(111)
	plt.minorticks_on()

        # Get the total number of ions
        nions = self.get_total_number_of_ions()

        # Get DoS per atom
        dosup = dos['up']/nions
        dosdown = -dos['down']/nions

        # Plot 
        plot.plot(energy, dosup, 'k', label='_nolegend_')
        plot.plot(energy, dosdown, 'k', label='_nolegend_')
        # plot.plot(energy, magnetization, 'r', label='Magnetization')

        # Set axis limits
        if axis_limits is not None:
          xmin, xmax, ymin, ymax = axis_limits
        else:
          xmax = np.max(energy)
          xmin = np.min(energy)
          ymax = max(np.max(dosup), np.max(-dosdown))
          ymin = -ymax

        # Plot the cross hairs
        plot.plot([xmin, xmax], [0.0, 0.0], 'k:', label='_nolegend_')
        plot.plot([0.0, 0.0], [ymin, ymax], 'k:', label='_nolegend_')

        plot.axis ([xmin, xmax, ymin, ymax])

        # Annotate the upper and the lower halves
        pl.annotate(r'$up$', xycoords='axes fraction', textcoords='axes fraction', xy=(0.8,.8), xytext=(0.8,0.8))
        pl.annotate(r'$down$', xycoords='axes fraction', textcoords='axes fraction', xy=(0.8,.8), xytext=(0.8,0.2))

        pl.xlabel(xlabel)
        pl.ylabel(r'$DoS/atom$')
        pl.legend()

        if show:
          pl.show()
        else:
          pl.savefig('dos.png')
          pl.savefig('dos.eps')

    def plot_pdos(self):
        """Plot the partial DoS for a given DOSCAR
        Plots:
        """
        pass
        
if __name__ == '__main__':

  import argparse,os

  parser = argparse.ArgumentParser(
      description = 'Plotting/saving density of states')
#  parser.add_argument('-sp', '--ifSpin', action='store_true', default=False,
#      help='If a spin-polarized DOS calculation (ISPIN=2)? (def: False)')   
#  parser.add_argument('-lm', '--ifLM', action='store_true', default=False,
#      help='If a LM-decomposed DOS calcualtion (LORBIT=11)? (def: False)') 
  parser.add_argument('-e', '--efermi', action='store_true', default=False,
      help='subtract efermi or not?') 
  parser.add_argument('-f', '--sfermi', action='store_true', default=False,
      help='shift energy by efermi or not?')
  parser.add_argument('-s', '--shift', default=0.0, type=float,
      help='shift energy by given value')
  parser.add_argument('-sh', '--show', action='store_true', default=False,
      help='show the plot or not?')
  parser.add_argument('-i', '--file', default='DOSCAR',
      help='Input file, def: DOSCAR')
  parser.add_argument(
      '-l', 
      '--axis_limits',
      nargs = 4,
      metavar = 'i',
      type = float,
      help='set the limits on the axes in the form xmin xmax ymin ymax')
  parser.add_argument('-t', '--test', action='store_true')
  args = parser.parse_args()

  ispin=2;lorbit=0
#  if args.ifSpin: ispin=2 
#  else: ispin=1
#  if args.ifLM: lorbit=11 
#  else: lorbit=0
#  print ispin,lorbit

  if args.test:
      import doctest
      doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=True)

  doscar().read(args.file).plot_dos(ispin,lorbit,sfermi=args.sfermi, show=args.show,
          axis_limits=args.axis_limits, shift=args.shift)

  os.system("display dos.png &")
