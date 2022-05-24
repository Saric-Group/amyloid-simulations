# encoding: utf-8
'''
This application does an NVE+Langevin LAMMPS simulation of spherocylinder-like rods
(defined in a .cfg file) using the "lammps_multistate_rods" library, with some rods
preassembled in a fibril (using "tools/prepare_fibril.py).
The initial locations of the rods are at SC lattice points defined by the input params,
excluding the fibril region, and their orientations are randomly determined at each
insertion point.

Created on 16 Mar 2018

@author: Eugen Rožić
'''
from mpi4py import MPI

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()

import os
import argparse

parser = argparse.ArgumentParser(description='Program for NVE+Langevin hybrid LAMMPS'\
                                 ' simulation of spherocylinder-like rods, using the'\
                                 ' "lammps_multistate_rods" library, with a preassembled fibril.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('cfg_file',
                    help='path to the "lammps_multistate_rods" model configuration file')
parser.add_argument('run_file',
                    help='path to the run configuration file')
parser.add_argument('simlen', type=int,
                    help='the length of the simulation')

parser.add_argument('--seed', type=int,
                    help='the seed for random number generators')
parser.add_argument('--out', type=str, default=None,
                    help='name/path for the output folder (defaults to cfg_file path w/o ext)')

parser.add_argument('-o', '--output_freq', type=int,
                    help='configuration output frequency (in MD steps);'\
                    ' default behavior is after every batch of MC moves')
parser.add_argument('-v', '--verbose', action='store_true',
                    help="print all LAMMPS output to stdout")

args = parser.parse_args()

if not args.cfg_file.endswith('.cfg'):
    raise Exception('Model configuration file (first arg) has to end with ".cfg"!')

if not args.run_file.endswith('.run'):
    raise Exception('Run configuration file (second arg) has to end with ".run"!')

if args.seed is None:
    seed = 0
    if mpi_rank == 0:
        import time
        seed = int((time.time() % 1)*1000000)
        print "WARNING: no seed given explicitly; using:", seed
    seed = mpi_comm.bcast(seed, root = 0)
else:
    seed = args.seed

if args.out is None:
    output_folder = os.path.splitext(args.cfg_file)[0]
else:
    output_folder = args.out
#========================================================================================
from lammps import PyLammps
import lammps_multistate_rods as rods
import lammps_multistate_rods.tools as rods_tools

if not os.path.exists(output_folder) and mpi_rank == 0:
    os.makedirs(output_folder)
    
run_filename = os.path.splitext(os.path.basename(args.run_file))[0]
sim_ID = '{:s}_{:d}'.format(run_filename, seed)
    
dump_filename = sim_ID+'.dump'
dump_path = os.path.join(output_folder, dump_filename)

log_filename = '{:d}.lammps'.format(seed)
log_path = os.path.join(output_folder, log_filename)

run_args = rods.Rod_params() # any object with __dict__ would do
if mpi_rank == 0:
    execfile(args.run_file,
             {'__builtins__' : None, 'True' : True, 'False' : False, 'None' : None},
             vars(run_args))
run_args = mpi_comm.bcast(run_args, root = 0)

out_freq = args.output_freq if args.output_freq != None else run_args.mc_every

if mpi_rank == 0:
    py_lmp = PyLammps(cmdargs = ['-echo', 'both'], comm = mpi_comm, verbose = args.verbose)
else:
    py_lmp = PyLammps(cmdargs = ['-echo', 'both'], comm = mpi_comm)
py_lmp.log('"' + log_path + '"')

rod_params = rods.Rod_params()
if mpi_rank == 0:
    rod_params.from_file(args.cfg_file);
rod_params = mpi_comm.bcast(rod_params, root = 0)

# CREATE BASE OBJECTS
model = rods.Rod_model(rod_params)
simulation = rods.Simulation(py_lmp, model, run_args.temp, seed, output_folder)

py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p p")
py_lmp.lattice("sc", 1/(run_args.cell_size**3))
py_lmp.region("box", "block", -run_args.num_cells / 2, run_args.num_cells / 2,
                              -run_args.num_cells / 2, run_args.num_cells / 2,
                              -run_args.num_cells / 2, run_args.num_cells / 2)
simulation.setup("box")

# create fibril
fibril_edges = []
fibril_rods = []
if mpi_rank == 0:
    fibril_edges = rods_tools.prepare_fibril(model, run_args.seed_size, run_args.seed_phi,
                                             run_args.seed_theta, run_args.seed_r0,
                                             data = fibril_rods)
fibril_edges = mpi_comm.bcast(fibril_edges, root = 0)
fibril_rods = mpi_comm.bcast(fibril_rods, root = 0)

simulation.create_rods(state_ID = model.num_states - 1, exact = [fibril_rods])

# create other rods
box_size = run_args.num_cells * run_args.cell_size
xmin = fibril_edges[0][0] - model.rod_length / 2
xmax = fibril_edges[0][1] + model.rod_length / 2
ymin = fibril_edges[1][0] - model.rod_length / 2
ymax = fibril_edges[1][1] + model.rod_length / 2
zmin = fibril_edges[2][0] - model.rod_length / 2
zmax = fibril_edges[2][1] + model.rod_length / 2
py_lmp.region("fibril", "block", xmin, xmax, ymin, ymax, zmin, zmax, "units box")
#TODO py_lmp.region("box_minus_fibril", "subtract", 2, "box", "fibril")
py_lmp.region("left", "block",
              'EDGE', xmin,
              'EDGE', 'EDGE',
              'EDGE', 'EDGE',
              "units box")
py_lmp.region("right", "block",
               xmax, 'EDGE',
              'EDGE', 'EDGE',
              'EDGE', 'EDGE',
              "units box")
py_lmp.region("front", "block",
              'EDGE', 'EDGE',
              'EDGE', ymin,
              'EDGE', 'EDGE',
              "units box")
py_lmp.region("back", "block",
              'EDGE', 'EDGE',
               ymax, 'EDGE',
              'EDGE', 'EDGE',
              "units box")
py_lmp.region("down", "block",
              'EDGE', 'EDGE',
              'EDGE', 'EDGE',
              'EDGE', zmin,
              "units box")
py_lmp.region("up", "block",
              'EDGE', 'EDGE',
              'EDGE', 'EDGE',
               zmax, 'EDGE',
              "units box")
py_lmp.region("box_minus_fibril", "union", 6, "up", "down", "front", "back", "left", "right")
simulation.create_rods(region = ["box_minus_fibril"])

# DYNAMICS
py_lmp.fix("thermostat", "all", "langevin",
           run_args.temp, run_args.temp, run_args.damp, seed)#, "zero yes")

simulation.set_rod_dynamics("nve", opt = ["mol", model.rod_states[0]])

py_lmp.neigh_modify("every 1 delay 1")
py_lmp.timestep(run_args.dt)

# THERMALIZE INITIAL CONFIGURATION
simulation.deactivate_state(0, vx_eps=5.0)
py_lmp.run(10000)
simulation.activate_state(0)
py_lmp.reset_timestep(0)

# STATE-CHANGING & CONCENTRATION FIXES
if model.num_states <= 1 or run_args.mc_tries <= 0:
    raise Exception("Multiple states need to exist and MC moves need to be made for fibrils to grow!")
mc_tries = int(run_args.mc_tries * simulation.rods_count())
mc_exchange_tries = int(0.01 * mc_tries + 1)
concentration = run_args.conc / run_args.cell_size**3
overlap = (2.1 * model.rod_radius) / run_args.cell_size
simulation.set_state_transitions(run_args.mc_every, mc_tries)#, opt = ["full_energy"])    
simulation.set_state_concentration(0, concentration, run_args.mc_every, mc_exchange_tries,
                                   opt = ["overlap_cutoff", overlap])

# GROUPS & COMPUTES
if hasattr(run_args, 'label_fibrils'):
    fibril_group = 'beta_patches'
    beta_active_patch_types = sorted(filter(lambda t: (t in model.active_bead_types) and\
                                                      (t not in model.body_bead_types),
                                            model.state_bead_types[1]))
    py_lmp.variable(fibril_group, 'atom', '"' + 
                                          '||'.join(['(type == {:d})'.format(t)
                                                     for t in beta_active_patch_types]) + 
                                          '"')
    py_lmp.group(fibril_group, 'dynamic', simulation.rods_group, 'var', fibril_group,
                 'every', out_freq)
    fibril_compute = "fibril_ID"
    if hasattr(run_args, 'fibril_cutoff'):
        fibril_cutoff = run_args.fibril_cutoff
    else:
        fibril_cutoff = 0
        i = -1
        for t1 in beta_active_patch_types:
            i += 1
            for t2 in beta_active_patch_types[i:]:
                try:
                    int_key = model.eps[(t1,t2)][1]
                except:
                    continue
                int_range = model.int_types[int_key][1]
                cutoff = model.bead_radii[t1] + model.bead_radii[t2] + int_range*2/3
                if cutoff > fibril_cutoff:
                    fibril_cutoff = cutoff
    py_lmp.compute(fibril_compute, fibril_group, 'aggregate/atom', fibril_cutoff)

# OUTPUT
dump_elems = "id x y z type mol"
try:
    dump_elems += " c_"+fibril_compute
except:
    pass
py_lmp.dump("dump_cmd", "all", "custom", out_freq, '"' + dump_path + '"', dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
py_lmp.thermo_style("custom", "step atoms", "pe temp",
                    " ".join(["v_{}".format(group_var)
                              for group_var in simulation.state_group_vars]),
                    "f_{}[2]".format(simulation.state_trans_fix), # state change successes
                    "f_{}[1]".format(simulation.state_trans_fix)) # state change attempts)
py_lmp.thermo(out_freq)

# RUN ...
py_lmp.run(args.simlen)

MPI.Finalize()
