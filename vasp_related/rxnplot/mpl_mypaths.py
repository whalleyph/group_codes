#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

# line
def line(ax, start, end, **kwargs):
    line = Path([start, end], [Path.MOVETO, Path.LINETO])
    patch = patches.PathPatch(line, **kwargs)
    ax.add_patch(patch)
    print 'added line from (%.2f, %.2f) to (%.2f, %.2f)' % (
            start[0], start[1], end[0], end[1])
    return ax

# Testing
if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = line(ax, (0.0, 0.0), (1.0, 0.0), ls='dotted')
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 1)
    plt.show()
