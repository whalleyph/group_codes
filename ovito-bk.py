#!/bin/env python3
from sys import exit

from ovito.io import import_file, export_file
#from ovito.data import SurfaceMesh, SimulationCell

from os import system
from sys import exit,stdout,argv,version_info

import argparse,os.path,time#,re
#from ase.neighborlist import NeighborList
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #needed as OVITO uses the default QT5 backend
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec


initT=time.time()

parser = argparse.ArgumentParser(description='Script for converting between different file formats using the ASE interface.')

parser.add_argument('-t','--typ', type=str,required=True, help='Analysis type, options: coord (coordination analysis including g(r)), voro (Solid volume based on Voronoi), XX')
parser.add_argument('-i','--inpf', nargs='*', type=str,required=True, help='Input file(s)')
#parser.add_argument('-ot','--otype', type=str,required=True, help='Output file type')
parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: input name is used as root.')
parser.add_argument('-skip','--skip', type=int,default=1,required=False, help='To skip steps while reading a trajectory. Def: All steps')
parser.add_argument('-all','--ifAll', default=False,action='store_true', help='Consider all the geometries in the input file (e.g. trajectories). Def: only last geometry is considered.')
parser.add_argument('-co','--cutoff', type=float,default=1.85,help='Cutoff distance for coordiantion analysis and Voronoi volume caclulations')
parser.add_argument('-ns','--noSubs', default=False,action='store_true', help='Do not include the substrate (atoms in the inital frame) in the coord analysis.')
parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='Overwrite if output file exists. Def: No')
args = parser.parse_args()


# Load a simulation snapshot of a Cu-Zr metallic glass.
#pipeline = import_file("input/simulation.0.dump")
pipeline = import_file(args.inpf) #can import multiple files at once
print('%d steps were read from %s, skipping %d steps.'%(pipeline.source.num_frames,args.inpf[0],args.skip));stdout.flush()

"""
ext=args.otype.split("/")[-1]
if args.outf: outf=args.outf+"."+ext
        else: outf=inpf.split(".")[0]+"."+ext

if args.typ in ["file",'f']:
    for frame in range(0,pipeline.source.num_frames,args.skip):
        """

if args.typ in ["coord",'c']:
    from ovito.modifiers import CoordinationAnalysisModifier
    # Insert the modifier into the pipeline:
    modifier = CoordinationAnalysisModifier(cutoff = args.cutoff, number_of_bins = 200,partial=0)
    pipeline.modifiers.append(modifier)
    
    font = {'family': 'serif', 'size': 18}
    plt.rc('font', **font)
    fig = plt.figure(figsize=(11.69, 8.27))
    gs = GridSpec(1, 1)#, width_ratios=[1, 2],height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0])  #Top subplot
    #ax2 = plt.subplot(gs[1]  , sharey=ax1)
    ax2=ax1.twinx()
    
    #if args.ifAll:
    btype=[]
    #atoms=ase.io.read(args.inpf[0])
    data = pipeline.compute(0)#initial step
    if args.noSubs:subs=range(data.particles.count) #subs=[at.index for at in atoms[0]]
    else:subs=[]
    #print(subs)
    
    # Initialize array for accumulated RDF histogram to zero:
    total_rdf = np.zeros((modifier.number_of_bins, 2))
    btype=[]
    conv=1.6605391 #amu/A^3 to g/cm^3
    for frame in range(0,pipeline.source.num_frames,args.skip):
        data = pipeline.compute(frame)
        # Accumulate RDF histograms:
        total_rdf += data.tables['coordination-rdf'].xy()

        #z_coords=data.particles.position.z
        vol = abs(np.linalg.det(data.cell[0:3,0:3])) #data.cell.volume is equivalent.
        dens= data.particles.count*12.01/vol*conv #based on box volume #needs to take mass from the input, but input is always H for Lammps trajectory.
        coord=data.particles.coordination #for each atom
        sp=0;sp2=0;sp3=0;oth=0
        #"""
        if args.noSubs:
            for at in range(len(coord)):
                if at in subs:  continue
                if coord[at] == 2: sp+=1
                elif coord[at] == 3: sp2+=1
                elif coord[at] == 4: sp3+=1
                elif coord[at] <2 or coord[at]>4: oth+=1
        else:
            sp=np.count_nonzero(coord==2)
            sp2=np.count_nonzero(coord==3)
            sp3=np.count_nonzero(coord==4)
            oth=np.count_nonzero((coord<=1) | (coord>4 ))
            
        tot=float(sp+sp2+sp3+oth)

        if tot==0: btype.append([frame,0.,0.,0.,0.,dens])#frame*args.skip
        else:btype.append([frame,sp3/tot,sp2/tot,sp/tot,oth/tot,dens])
        
    # Averaging:
    total_rdf /= pipeline.source.num_frames
    # Export the average RDF to a text file:
    np.savetxt("rdf.txt", total_rdf)
    
    inpf=args.inpf[0]
    with open('%s.data'%(inpf.split('.')[0]), 'w') as f:
        #with open('hybrid.data', 'w') as f:
        f.write('#timeStep, sp3, sp2, sp, others\n')
        for b in btype:
            f.write('%8d %.3f %.3f %.3f %.3f \n'%(b[0],b[1],b[2],b[3],b[4]))
            
    print("Elapsed time: %.2f sec."%( time.time()-initT))

    btype=np.array(btype)
    ax1.plot(btype.T[0],btype.T[1],label='sp3')
    ax1.plot(btype.T[0],btype.T[2],label='sp2')
    ax1.plot(btype.T[0],btype.T[3],label='sp')
    ax1.plot(btype.T[0],btype.T[4],label='others')
    ax1.set_ylabel('Fraction')
    ax1.set_xlabel('Timestep')
    ax2.set_ylabel('Density [g/cc]')
    ax2.plot(btype.T[0],btype.T[-1],label='density',color='k')
    
    ax1.legend()
    #plt.legend()
    #plt.ylim(0,1)
    plt.show()

    
elif args.typ in ["voro",'v']:
    from ovito.modifiers import ConstructSurfaceModifier,VoronoiAnalysisModifier
    
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
    pipeline.modifiers.append(ConstructSurfaceModifier(radius = args.cutoff)) #Should 1.82 for matching N2-BET

    # Let OVITO compute the results.
    for frame in range(0,pipeline.source.num_frames,args.skip):
        data = pipeline.compute(frame)
        #data = pipeline.compute()
        #mesh = data.surfaces['surface']

        # Query computed surface properties:
        print("Surface area: %f" % data.attributes['ConstructSurfaceMesh.surface_area'])
        print("Solid volume: %f" % data.attributes['ConstructSurfaceMesh.solid_volume'])
        fraction = data.attributes['ConstructSurfaceMesh.solid_volume'] / data.cell.volume
        print("Solid volume fraction: %f" % fraction)

        # Export the surface triangle mesh to a VTK file.
        #export_file(mesh, "output/surface.vtk", "vtk/trimesh")
    
print("Elapsed time: %.2f sec."%( time.time()-initT))

exit()




###############
#Deleted Parts#
###############
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
