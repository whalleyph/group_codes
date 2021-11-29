#!/bin/env python3
# Import OVITO modules.
#from ovito.io import *
#from ovito.modifiers import *
from sys import exit
# Import NumPy module.
#import numpy

from ovito.io import import_file, export_file
#from ovito.data import SurfaceMesh, SimulationCell
from ovito.modifiers import ConstructSurfaceModifier,VoronoiAnalysisModifier

#from os import system
from sys import argv

# Load a simulation snapshot of a Cu-Zr metallic glass.
#pipeline = import_file("input/simulation.0.dump")
pipeline = import_file(argv[1])

# Set atomic radii (required for polydisperse Voronoi tessellation).
atom_types = pipeline.source.data.particles['Particle Type'].types
#atom_types[0].radius = 1.35   # Cu atomic radius (atom type 1 in input file)
#atom_types[1].radius = 1.55   # Zr atomic radius (atom type 2 in input file)
atom_types[0].radius = 0.70   # C atomic radius (atom type 1 in input file)
#get it from neighborlist.natural_cutoffs

# Set up the Voronoi analysis modifier.
voro = VoronoiAnalysisModifier(
    compute_indices = True,
    use_radii = True,
    edge_threshold = 0.1
)
pipeline.modifiers.append(voro)
pipeline.modifiers.append(ConstructSurfaceModifier(radius = 2.9))

# Let OVITO compute the results.
data = pipeline.compute()
#mesh = data.surfaces['surface']

# Query computed surface properties:
print("Surface area: %f" % data.attributes['ConstructSurfaceMesh.surface_area'])
print("Solid volume: %f" % data.attributes['ConstructSurfaceMesh.solid_volume'])
fraction = data.attributes['ConstructSurfaceMesh.solid_volume'] / data.cell.volume
print("Solid volume fraction: %f" % fraction)

# Export the surface triangle mesh to a VTK file.
#export_file(mesh, "output/surface.vtk", "vtk/trimesh")

exit()

# Access computed Voronoi indices.
# This is an (N) x (M) array, where M is the maximum face order.
voro_indices = data.particles['Voronoi Index']

# This helper function takes a two-dimensional array and computes a frequency
# histogram of the data rows using some NumPy magic.
# It returns two arrays (of equal length):
#    1. The list of unique data rows from the input array
#    2. The number of occurences of each unique row
# Both arrays are sorted in descending order such that the most frequent rows
# are listed first.
def row_histogram(a):
    ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
    counts = numpy.bincount(inverse)
    sort_indices = numpy.argsort(counts)[::-1]
    return (a[indices[sort_indices]], counts[sort_indices])

# Compute frequency histogram.
unique_indices, counts = row_histogram(voro_indices)

# Print the ten most frequent histogram entries.
for i in range(10):
    print("%s\t%i\t(%.1f %%)" % (tuple(unique_indices[i]),
                                 counts[i],
                                 100.0*float(counts[i])/len(voro_indices)))
