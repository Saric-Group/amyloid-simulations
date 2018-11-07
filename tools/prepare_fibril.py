# encoding: utf-8
'''
A module that can be called as a program to generate a file with locations and rotations
of rods for a preassembled fibril. This file can be read to generate a fibril by the
"create_rods" method of a lammps_multistate_rods.Simulation class instance.

Created on 10 Aug 2018

@author: Eugen Rožić
'''
import numpy as np
import pyquaternion
import lammps_multistate_rods as rods
import os

def prepare_fibril(model, N, phi, theta, r0, output):
    '''
    This method outputs a file with data to create a preformed fibril to "output".
    
    model : A lammps_multistate_rods.Model clas instance
    phi : azimuth angle (with respect to y-axis)
    theta : elevation angle (with respect to x-y plane)
    r0 : centre of the fibril (a triplet)
    output : path/name of the output file
    
    returns : a triplet of (min,max) coordinate pairs, one for each dimension, of
    extremal locations among all centers of the rods
    '''
    
    rod_radius = model.rod_radius
    r0 = np.array(r0)

    R_z = pyquaternion.Quaternion(axis=[0,0,1], degrees=phi)
    R_x = pyquaternion.Quaternion(axis=[1,0,0], degrees=theta)
    if theta > 0:
        R_x_inv = pyquaternion.Quaternion(axis=[1,0,0], degrees=theta-180)
    else:
        R_x_inv = pyquaternion.Quaternion(axis=[1,0,0], degrees=theta+180)

    # correct composite rotation is to first rotate around z then around x', which is equivalent to
    # rotations first around x then around z for the same angles
    R_tot = R_z * R_x
    R_tot_inv = R_z * R_x_inv

    locations = [None]*N # list of position 3-vectors
    rotations = [None]*N # list of tuples of a position 3-vector and angle
    mins = [np.infty]*3
    maxs = [-np.infty]*3

    for i in range(N):
        if i % 2 == 0:
            locations[i] = np.array([0, (i-N/2)*rod_radius, -rod_radius])#+0.27692])
            locations[i] = R_tot.rotate(locations[i]) + r0
            if R_tot.angle == 0.0:
                rotations[i] = (R_tot.angle, [1.0, 0.0, 0.0])
            else:
                rotations[i] = (R_tot.angle, R_tot.axis)
        else:
            locations[i] = np.array([0, (i-N/2)*rod_radius, +rod_radius])#-0.27692])
            locations[i] = R_tot.rotate(locations[i]) + r0
            if R_tot_inv.angle == 0.0:
                rotations[i] = (R_tot_inv.angle, [1.0, 0.0, 0.0])
            else:
                rotations[i] = (R_tot_inv.angle, R_tot_inv.axis)
        
        for j in range(3):
            if locations[i][j] > maxs[j]:
                maxs[j] = locations[i][j]
            if locations[i][j] < mins[j]:
                mins[j] = locations[i][j]

    with open(output, 'w') as output_file:
        output_file.write('monomers: {:d}\n\n'.format(N))
        for loc, rot in zip(locations, rotations):
            output_file.write('{:.2f} {:.2f} {:.2f} {:.2f} {:.3f} {:.3f} {:.3f}\n'.format(
                loc[0], loc[1], loc[2], rot[0], *rot[1]))
        output_file.write('\n')
    
    return zip(mins, maxs)
            
if __name__ == '__main__': 

    import argparse

    parser = argparse.ArgumentParser(
        description='Application for preparing a fibril formation for rods/monomers',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config_file',
                        help='path to the "lammps_multistate_rods" model config file')
    parser.add_argument('phi', type=float,
                        help='azimuth angle (from y-axis); in degrees [0-360>')
    parser.add_argument('theta', type=float,
                        help='elevation angle (from x-y plane); in degrees [-90,90]')
    parser.add_argument('output', type=str,
                        help='path/name for the output file')
    parser.add_argument('--r0', nargs=3, type=float, default=[0.,0.,0.],
                        help='the location of the center of the fibril')
    parser.add_argument('-N', type=int, default=20,
                        help='number of monomers in the fibril')
    parser.add_argument('--test', type=int, default=None,
                        help='the length of a test run of the fibril (nve+langevin)')
    args = parser.parse_args()

    model = rods.Model(args.config_file)

    fibril_edges = prepare_fibril(model, args.N, args.phi, args.theta, args.r0, args.output)
    
    if args.test == None:
        quit()
    
    from lammps import PyLammps
    
    seed = 12345
    T = 1.0
    damp = 0.1
    output_dir = os.path.dirname(args.output)
    log_path = os.path.join(output_dir, "prepare_fibril.log")
    output_freq = 10
    dump_path = os.path.join(output_dir, 'prepare_fibril.dump')
    
    py_lmp = PyLammps(cmdargs=['-screen','none'])
    simulation = rods.Simulation(py_lmp, model, seed, output_dir, log_path=log_path, clusters=1.0)

    py_lmp.units("lj")
    py_lmp.dimension(3)
    py_lmp.boundary("p p p")
    
    xmin = fibril_edges[0][0] - model.rod_length
    xmax = fibril_edges[0][1] + model.rod_length
    ymin = fibril_edges[1][0] - model.rod_length
    ymax = fibril_edges[1][1] + model.rod_length
    zmin = fibril_edges[2][0] - model.rod_length
    zmax = fibril_edges[2][1] + model.rod_length
    py_lmp.region("box", "block", xmin, xmax, ymin, ymax, zmin, zmax)
    simulation.setup("box")

    simulation.create_rods(state_ID=model.num_states-1, file=args.output)

    py_lmp.fix("thermostat", "all", "langevin", T, T, damp, seed)
    simulation.set_rod_dynamics("nve")

    py_lmp.neigh_modify("every 1 delay 1")

    # OUTPUT
    dump_elems = "id x y z type mol c_"+simulation.cluster_compute
    py_lmp.dump("dump_cmd", "all", "custom", output_freq, '"'+dump_path+'"', dump_elems)
    py_lmp.dump_modify("dump_cmd", "sort id")
    py_lmp.thermo_style("custom", "step atoms", "pe temp")
    py_lmp.thermo(1000)

    py_lmp.command('run {:d}'.format(args.test))


