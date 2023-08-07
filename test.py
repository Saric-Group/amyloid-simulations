# encoding: utf-8
'''
A module for testing...

Created on 10 Aug 2018

@author: Eugen Rožić
'''
import os, sys

from lammps import PyLammps

import lammps_multistate_rods as rods
import lammps_multistate_rods.tools as rods_tools

cfg_file_loc = sys.argv[1]
script_loc = os.path.abspath(os.path.dirname())

params = rods.Rod_params().from_file(cfg_file_loc)
model = rods.Rod_model(params)
    
output_dir = os.path.join(os.path.dirname(cfg_file_loc), 'test_out')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
data_in = os.path.join(output_dir, 'data.in')
    
r0 = [0.0, 0.0, 0.0]
z_axis = [0.0, 0.0, 1.0]
phi = 45.0
theta = 35.0
N = 20

#box_size = 2*model.rod_length
fibril_rods = []
#fibril_rods.append(rods_tools.preparation.prepare_single(r0, phi, theta, out_path = data_in))
fibril_edges = rods_tools.preparation.prepare_fibril(model, N, phi, theta, r0,
                                                     data = fibril_rods, out_path = data_in)
    
seed = 12345
T = 1.0
damp = 0.1
log_path = os.path.join(output_dir, "test.log")
output_freq = 10
dump_path = os.path.join(output_dir, 'test.dump')
    
py_lmp = PyLammps(cmdargs = ['-echo','both'])
py_lmp.log('"' + log_path + '"')

simulation = rods.Simulation(py_lmp, model, T, seed, output_dir)

py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p p")
#py_lmp.neighbor(2.5, 'bin')
    
xmin = fibril_edges[0][0] - model.rod_length
xmax = fibril_edges[0][1] + model.rod_length
ymin = fibril_edges[1][0] - model.rod_length
ymax = fibril_edges[1][1] + model.rod_length
zmin = fibril_edges[2][0] - model.rod_length
zmax = fibril_edges[2][1] + model.rod_length
py_lmp.region("box", "block", xmin, xmax, ymin, ymax, zmin, zmax)
simulation.setup("box")

simulation.create_rods(state_ID = 1, exact = [fibril_rods])

py_lmp.fix("thermostat", "all", "langevin", T, T, damp, seed)
simulation.set_rod_dynamics("nve", opt = ["mol", model.rod_states[0]])

overlap = 2.1 * model.rod_radius
conc = 1. / model.rod_length**3 # one per cube of side = rod length
N_sol_approx = int((xmax - xmin) * (ymax - ymin) * (zmax - zmin) * conc + 1)
mc_every = 100

simulation.set_state_concentration(0, conc, mc_every, int(0.01 * N_sol_approx + 1),
                                   opt = ["overlap_cutoff", overlap])
simulation.set_state_transitions(mc_every, N + N_sol_approx)#, opt = ['auto_skin'])#, opt = ['full_energy'])

# OUTPUT
dump_elems = "id x y z type mol"# c_"+simulation.cluster_compute
py_lmp.dump("dump_cmd", "all", "custom", output_freq, '"' + dump_path + '"', dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
py_lmp.thermo_style("custom", "step atoms", "pe temp",
                    " ".join(["v_{}".format(group_var)
                              for group_var in simulation.state_group_vars]),
                    "f_{}[2]".format(simulation.state_trans_fix), # state change successes
                    "f_{}[1]".format(simulation.state_trans_fix)) # state change attempts
py_lmp.thermo(mc_every)

py_lmp.neigh_modify("every 1 delay 1")
py_lmp.timestep(0.01)
py_lmp.run(1000)
