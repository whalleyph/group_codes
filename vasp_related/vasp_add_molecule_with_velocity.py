#! /usr/bin/env python

import ase.io
import ase
import argparse
import numpy as np
import sys,os

from ase.geometry import find_mic

def minDist_ASE(arr1,arr2,latt): #needs np.array inputs
    if len(arr1)==4: arr1=arr1[1:4]
    if len(arr2)==4: arr2=arr2[1:4]
    D,x=find_mic(np.array([arr1-arr2]),latt,pbc=True)
    return x[0]

def get_speed_of_molecule(molecule, temperature,ntimes=1):
    # calculate the velocity of the molecule
    temperature = args.temperature
    gas_constant = 8.3144621 # kg m^2/s^2-K-mol
    avogadro_number = 6.022e23
    kB = gas_constant/avogadro_number
    #velocity = args.velocity
    mol_mass = ntimes* np.sum(molecule.get_masses()) * 1e-3 / avogadro_number # kg
    #mol_mass = np.sum(molecule.get_masses()) * 1e-3 / avogadro_number # kg
    if temperature:
        #print np.sum(molecule.get_masses())
        print(ntimes*np.sum(molecule.get_masses()))
        avg_velocity = np.sqrt(8*kB*temperature/(np.pi*mol_mass))     # m/s
        avg_velocity = avg_velocity*1e-5    # angstrom/fs
    elif args.velocity:
        avg_velocity = args.velocity
        temp=(avg_velocity*1e5)**2/(8*kB/(np.pi*mol_mass))
        print("Corresponding temp (K): ",temp)
    else:
        print('Either the temperature or the velocity should be supplied')
        sys.exit(-1)
    print('%.3e Angstrom / femtosecond' % (avg_velocity))
    ke = 0.5*mol_mass*((avg_velocity*1e5)**2) # Joule
    ke = ke*6.2415093e+18
    print('Kinetic energy: %.2f meV' % (ke*1e3))

    return avg_velocity

def get_random_ab(cell, ntimes, tol):
    ab = []
    ab.append(np.random.rand(2))
    for i in range(1,ntimes):
        not_accepted = True
        cnt=0
        while not_accepted:
            if cnt>100:print('not possible to add the adatoms with the given parameters....');break
            cnt+=1
            # generate_random_coordinate
            atmp, btmp = np.random.rand(2)
            #xtmp, ytmp, dummy = np.dot(cell.transpose(), np.array([atmp, btmp, 0.0]))
            tmp = np.dot(cell.transpose(), np.array([atmp, btmp, 0.0]))
            xtmp, ytmp, dummy = tmp
            #print('xtmp ytmp dummy', xtmp, ytmp, dummy)
            # compute distance with all previous coordinates
            r_prev = []
            for ir in ab:
                r1x, r1y, r1z = np.dot(cell.transpose(), np.array([ir[0], ir[1], 0.0]))
                r1 = np.array([r1x, r1y,0.])
                r2 = np.array([xtmp, ytmp,0.])
                #print('r2', r2)
                #r12 = np.linalg.norm(r2-r1)
                r12=minDist_ASE(r2,r1,cell)
                r_prev.append(r12)

            within_tol = True
            for ir in r_prev:
                if ir < tol:
                    within_tol = False
                    break

            if within_tol:
                not_accepted = False
                ab.append(np.array([atmp, btmp]))

    return np.array(ab)


def add_molecule_to_surface(surface, molecule, ntimes=1, nlayers=1,random=False, ab=None, tol=2.0):
    # add molecule to the surface
    new_surface = surface.copy()

    for l in range(1,nlayers+1):
        cell = new_surface.cell
        zmax = np.max(new_surface.get_positions()[:,2])
        mol_z = zmax + args.das
        mol_c = mol_z/cell[2,2]

        if random:
            ab = get_random_ab(cell, ntimes=ntimes, tol=tol)
        
        for i in range(ntimes):
            mol_a, mol_b = ab[i]
            mol_pos = np.dot(cell.transpose(), np.array([mol_a, mol_b, mol_c]))
            # translate molecule to new position
            molecule.translate(mol_pos - molecule.get_positions()[0])
            #Rotate randomly
            axes=['x','y','z']
            axis=axes[np.random.randint(3, size=1)[0]]
            angle=360*np.random.random()-180 #sample btw [-180,180)
            print(axis,angle)
            molecule.rotate(angle,axis,center='COP')
            
            new_molecule = ase.Atoms(molecule, cell=surface.cell)
            new_molecule.wrap()
            new_surface += new_molecule

    return new_surface

def add_velocities(surface, molecule, temperature, input_file, tmp_output_file,ntimes=1):
    ns = len(surface)
    nm = len(molecule)
    velocities=[]
    #This part is for getting the existing velocities from the input file.
    input_lines = open(input_file).readlines()
    for i, line in enumerate(input_lines):
        if 'cart' in line.lower() or 'direct' in line.lower():
            positions = input_lines[i+1:i+ns]
            velocities = input_lines[i+ns+2:i+ns+2+ns]
            if len(velocities)!=0: print("Velocities read from input file:");#print(velocities)
    
    tmp_output_lines = open(tmp_output_file).readlines()
    if len(velocities)==0: tmp_output_lines.append('\n')

    if  args.velocity: speed=args.velocity
    else: speed = get_speed_of_molecule(molecule, temperature,ntimes)

    if len(velocities)==0: #If no velocities found in the input.
        vel1 = ' %15.8e %15.8e %15.8e\n' % (0.0, 0.0, 0.0) #keep surface atoms fixed
        velocities.append(ns*vel1)
    vel = ' %15.8e %15.8e %15.8e\n' % (0.0, 0.0, -speed)
    velocities.append(ntimes*nm*vel)
    velocities.insert(0,"\n")
    #print velocities
    tmp_output_lines.extend(velocities)
    #print tmp_output_lines
    new_output_file_string = ''.join(tmp_output_lines)
    return new_output_file_string



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='')

    parser.add_argument('-i', '--input-file', default='POSCAR')
    parser.add_argument('-o', '--output-file', default='new.vasp')
    parser.add_argument('-m', '--molecule', required=1) #default='molecule.xyz')
    parser.add_argument('-n', '--nmolecules', type=int, default=1,
            help='Number of molecules to add per each layer')
    parser.add_argument('-nl', '--nlayers', type=int, default=1,
            help='Number of layers of molecules to add')
    parser.add_argument('-v', '--velocity', type=float, default=None,
            help='in Ang/fs')
    parser.add_argument('-t', '--temperature', type=float, default=None,
            help='temperature in K')
    parser.add_argument('-d', '--das', type=float, default=4.5,
            help='distance above surface')
    parser.add_argument('-ab', '--ab', nargs=2, type=float, default=None,
            help='position where the molecule should be positioned in the ab plane \
                  will be placed randomly if not specified')
    parser.add_argument('-tol', '--tolerance', type=float, default=2.0,
            help='minimum distance between each of the added molecules')

    args = parser.parse_args()

    surface = ase.io.read(args.input_file) #, format='vasp')
    molecule = ase.io.read(args.molecule) #, format='xyz')
    nmolecules = args.nmolecules
    ab = args.ab
    tol = args.tolerance
    print(surface,molecule)
    print(get_speed_of_molecule(molecule, args.temperature,nmolecules))

    if args.ab:
        new_surface = add_molecule_to_surface(surface, molecule, ntimes=nmolecules,nlayers=args.nlayers, ab=ab, tol=tol)
    else:
        new_surface = add_molecule_to_surface(surface, molecule, ntimes=nmolecules, nlayers=args.nlayers,random=True, tol=tol)

    ase.io.write('tmp_output', new_surface, format='vasp', vasp5=True, direct=False)

    # now append velocities to this newly written file
    #molecule_speed = get_speed_of_molecule(molecule, args.temperature)
    new_output_file_string = add_velocities(surface, molecule, args.temperature, args.input_file, 'tmp_output',nmolecules)

    open(args.output_file, 'w').write(new_output_file_string)
    os.system('rm tmp_output')
