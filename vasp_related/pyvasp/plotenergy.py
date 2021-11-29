#! /usr/bin/python

import sys
sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
import outcar as o
import archive as ar
import matplotlib.pyplot as plt

def get_energies (outcars):
    energies = []
    for outcar in outcars:
        out = o.outcar().read(outcar)
        vol = out.e0()
        energies.append(vol)

    return energies

def plot_energy(ax, energies):
    i = 0
    for energy in energies:
        x = range(i, i+len(energy))
        y = energy
        i = i + len(energy)
        ax.plot(x, y, 'o-')

    ax.set_xlabel('Step')
    ax.set_ylabel(r'$E0 \,\left[eV\right]$')
    return ax

if __name__ == '__main__':

    tgzs = sys.argv[1:]
    energies = get_energies(tgzs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_energy(ax, energies)
    plt.show()
