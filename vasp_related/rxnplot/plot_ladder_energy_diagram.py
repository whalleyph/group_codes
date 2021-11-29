#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import argparse
import re

## matplotlib customizations
import matplotlib as mpl
mpl.rc('font', size=12)#, weight='bold')
mpl.rc('axes', linewidth=1)
mpl.rc('figure.subplot', bottom=0.025, top=0.975, left=0.25, right=0.95)

# All measurements are in cm.

def setup_canvas(energies, height=10, emin=None, emax=None):

        if not emin:
                emin = np.min(energies)
        if not emax:
                emax = np.max(energies)

        rng = emax - emin               # range

        scale = height/rng        # scale

        ymin = scale*emin
        ymax = scale*emax

        mypainter = graph.axis.painter.regular(
                outerticklength=graph.axis.painter.ticklength.normal)

        c = graph.axis.pathaxis(path.line(0, ymax, 0, ymin),
                graph.axis.linear(min=emax, max=emin,
                    title=r'E [kJ/mole]', painter=mypainter))

        return c, scale

def plot_ladder(canvas, energies, labels, offsets, xstart=0.5,
                elw=1.5, text_size=-1, label_offset=0.25):

        for eng, label, offset in zip(energies, labels, offsets):

                xend = xstart + elw

                p = path.line(xstart, eng, xend, eng)

                label = '\bf{' + label + '}'

                if label and offset:
                        c.text(xend+label_offset+offset, eng, label, [text.size(text_size)])

                if label and not offset:
                        c.text(xend+label_offset, eng, label, [text.size(text_size)])

                c.stroke(p, [style.linewidth(0.1)])
                #c.stroke(p)

        return c

def state(y, width=1):
    state = Path([(0,y),(width,y)], [Path.MOVETO, Path.LINETO])
    return state

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--files', nargs='*', required=True,
                        help = 'Files containing energies')
        parser.add_argument('-emin', '--emin', type=float,
                        help = 'Energy minimum')
        parser.add_argument('-emax', '--emax', type=float,
                        help = 'Energy maximum')

        args = parser.parse_args()

        files = args.files
        emin, emax = args.emin, args.emax

        outfiles = [f.replace('.dat', '') for f in files]
        outfile = '_'.join(outfiles)
        
        ndata = len(files)
        energies, labels, offsets = [], [], []
        for f in files:
                d = np.genfromtxt(f, delimiter=',',
                                    dtype=[('energy', 'f2'),
                                           ('label', 'a40'),
                                           ('offset', 'f2')],
                                    filling_values=(0.0, '', 0.0))

                e, l, o = d['energy'], d['label'], d['offset']
                # Match all numbers and replace them as subscripts
                nch = re.compile(r'([C,H])([0-9])')
                #l = [nch.sub(r'\1$_{\2}$', label.strip()) for label in l]
                l = [nch.sub(r'\1_{\2}', label.strip()) for label in l]
                l = [r'$\mathbf{' + label + '}$' for label in l]
                energies.append(e)
                labels.append(l)
                offsets.append(o)

        energies = np.concatenate(energies)
        labels = np.concatenate(labels)
        offsets = np.concatenate(offsets)

        print energies
        print labels
        print offsets

        states = []
        for energy in energies:
            states.append(state(energy, width=0.5))

        fig = plt.figure(figsize=(3.5,7))
        ax = fig.add_subplot(111)
        for state in states:
            patch = patches.PathPatch(state, lw=2)
            ax.add_patch(patch)

        # add labels
        for e, l in zip(energies, labels):
            ax.text(.55, e, l)#, size='x-large')
        ax.set_xlim(-0.25, 1.25)
        ax.set_ylim(emin, emax)
        ax.set_ylabel(r'E [kJ/mole]')
        ax.set_xticks([])
        plt.savefig(outfile+'.png')
        plt.savefig(outfile+'.svg')
        plt.show()
        #print states

        #c, scale = setup_canvas(energies, height=12, emin=emin, emax=emax)
        #plot_ladder(c, scale*energies, labels, offsets)
        #c.writeEPSfile(outfile, paperformat=document.paperformat.A4, fittosize=1)
