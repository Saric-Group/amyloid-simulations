# encoding: utf-8
'''
Created on 10 Aug 2018

@author: Eugen Rožić
'''
import numpy as np
import pyquaternion
import lammps_multistate_rods as rods

import argparse

parser = argparse.ArgumentParser(description='''
        Application for preparing a fibril formation for rods/monomers''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('config_file', help='path to the "lammps_multistate_rods" model config file')
parser.add_argument('phi', type=float, help='azimuth angle (from y-axis); in degrees [0-360>')
parser.add_argument('theta', type=float, help='elevation angle (from x-y plane); in degrees [-90,90]')
parser.add_argument('output', type=str, help='path/name for the output file')
parser.add_argument('--center', nargs=3, type=float, default=[0.,0.,0.],
                    help='the location of the center of the fibril')
parser.add_argument('-N', type=int, default=50, help='number of monomers in the fibril')
args = parser.parse_args()

model = rods.Model(args.config_file)
rod_radius = model.rod_radius
phi = args.phi
theta = args.theta
output = args.output
x0, y0, z0 = r0 = np.array(args.center)
N = args.N

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

for i in range(N):
    if i % 2 == 0:
        locations[i] = np.array([0, (i-N/2)*rod_radius, -rod_radius+0.27692])
        locations[i] = R_tot.rotate(locations[i]) + r0
        if R_tot.angle == 0.0:
            rotations[i] = (R_tot.angle, [1.0, 0.0, 0.0])
        else:
            rotations[i] = (R_tot.angle, R_tot.axis)
    else:
        locations[i] = np.array([0, (i-N/2)*rod_radius, +rod_radius-0.27692])
        locations[i] = R_tot.rotate(locations[i]) + r0
        if R_tot_inv.angle == 0.0:
            rotations[i] = (R_tot_inv.angle, [1.0, 0.0, 0.0])
        else:
            rotations[i] = (R_tot_inv.angle, R_tot_inv.axis)

with open(output, 'w') as output_file:
    output_file.write('monomers: {:d}\n\n'.format(N))
    for loc, rot in zip(locations, rotations):
        output_file.write('{:.2f} {:.2f} {:.2f} {:.2f} {:.3f} {:.3f} {:.3f}\n'.format(
            loc[0], loc[1], loc[2], rot[0], *rot[1]))

import time
seed = int((time.time() % 1)*1000000)
temp = 1.0
damp = 0.1
   
from lammps import PyLammps
    
py_lmp = PyLammps()
    
simulation = rods.Simulation(py_lmp, model, seed, temp, 'data/test')
    
py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p p")
py_lmp.lattice("sc", 1./(10.0**3))
box_size = 10.0
py_lmp.region("box", "block", -box_size / 2, box_size / 2, -box_size / 2, box_size / 2, -box_size / 2, box_size / 2)
     
simulation.setup("box") # a lot of customisation options available here
    
simulation.create_rods(state_ID=1, file=output)
  
py_lmp.fix("thermostat", "all", "langevin", temp, temp, damp, seed)#, "zero yes")
simulation.set_rod_dynamics("nve")
  
# OUTPUT
dump_elems = "id x y z type mol c_"+simulation.cluster_compute
py_lmp.dump("dump_cmd", "all", "custom", 100, 'data/test/fibril.dump', dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
  
py_lmp.thermo_style("custom", "step atoms", "pe temp")
py_lmp.thermo(1000)
  
py_lmp.command('run 10000')
    