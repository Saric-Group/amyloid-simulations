# encoding: utf-8
'''
A module for testing...

Created on 10 Aug 2018

@author: Eugen Rožić
'''
import os
from lammps import PyLammps
import lammps_multistate_rods as rods
import lammps_multistate_rods.tools as rods_tools
            
model = rods.Model('../data/test/5p_cross-grad_v3.cfg')
    
output_dir = os.path.dirname('../data/test/test_out/')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
data_in = os.path.join(output_dir, 'data.in')
    
r0 = [0.0, 0.0, 0.0]
z_axis = [0.0, 0.0, 1.0]
phi = 45.0
theta = 35.0
N = 20

#box_size = 2*model.rod_length
#rods_tools.prepare_single(r0, 0, 0, data_in)
fibril_edges = rods_tools.prepare_fibril(model, N, phi, theta, r0, data_in)
rods_state = 1
    
seed = 12345
T = 1.0
damp = 0.1
log_path = os.path.join(output_dir, "test.log")
output_freq = 10
dump_path = os.path.join(output_dir, 'test.dump')
    
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

simulation.create_rods(state_ID=rods_state, file=data_in)

py_lmp.fix("thermostat", "all", "langevin", T, T, damp, seed)
simulation.set_rod_dynamics("nve")

py_lmp.neigh_modify("every 1 delay 1")

# OUTPUT
dump_elems = "id x y z type mol c_"+simulation.cluster_compute
py_lmp.dump("dump_cmd", "all", "custom", output_freq, '"'+dump_path+'"', dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
py_lmp.thermo_style("custom", "step atoms", "pe temp")
py_lmp.thermo(1000)

py_lmp.timestep(0.01)
run_length = 1000
py_lmp.command('run {:d}'.format(run_length))
