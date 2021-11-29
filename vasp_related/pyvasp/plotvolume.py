#! /usr/bin/python

import sys
sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
import outcar as o
import archive as ar
import matplotlib.pyplot as plt
import tarfile

def get_volumes (outcars):
    volumes = []
    for outcar in outcars:
        out = o.outcar().read(outcar)
        vol = out.get_volume()
        volumes.append(vol)

    return volumes

def plot_volume(ax, volumes):
    i = 0
    for volume in volumes:
        x = range(i, i+len(volume))
        y = volume
        i = i + len(volume)
        ax.plot(x, y, 'o-')

    ax.set_xlabel('Step')
    ax.set_ylabel(r'$Volume\,\left[\AA^{3}\right]$')
    return ax

if __name__ == '__main__':

    tgzs = sys.argv[1:]

    volumes = get_volumes(tgzs)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_volume(ax, volumes)
    plt.show()
