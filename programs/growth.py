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
parser.add_argument('-s', '--silent', action='store_true',
                    help="doesn't print anything to stdout")

args = parser.parse_args()

if not args.cfg_file.endswith('.cfg'):
    raise Exception('Model configuration file (first arg) has to end with ".cfg"!')

if not args.run_file.endswith('.run'):
    raise Exception('Run configuration file (second arg) has to end with ".run"!')

if args.seed is None:
    import time
    seed = int((time.time() % 1)*1000000)
    print "WARNING: no seed given explicitly; using:", seed
else:
    seed = args.seed

if args.out is None:
    output_folder = os.path.splitext(args.cfg_file)[0]
else:
    output_folder = args.out
#========================================================================================

#from mpi4py import MPI #TODO make MPI work...
from lammps import PyLammps
import lammps_multistate_rods as rods
import lammps_multistate_rods.tools as rods_tools

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
run_filename = os.path.splitext(os.path.basename(args.run_file))[0]
sim_ID = '{:s}_{:d}'.format(run_filename, seed)
    
dump_filename = sim_ID+'.dump'
dump_path = os.path.join(output_folder, dump_filename)

log_filename = '{:d}.lammps'.format(seed)
log_path = os.path.join(output_folder, log_filename)

run_args = rods.rod_model.Params()
execfile(args.run_file,
         {'__builtins__' : None, 'True' : True, 'False' : False, 'None' : None},
         vars(run_args))

out_freq = args.output_freq if args.output_freq != None else run_args.run_length

py_lmp = PyLammps(cmdargs=['-screen','none'])
py_lmp.log('"'+log_path+'"')
model = rods.Rod_model(args.cfg_file)
simulation = rods.Simulation(py_lmp, model, seed, output_folder)

py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p p")
py_lmp.lattice("sc", 1/(run_args.cell_size**3))
py_lmp.region("box", "block", -run_args.num_cells / 2, run_args.num_cells / 2,
                              -run_args.num_cells / 2, run_args.num_cells / 2,
                              -run_args.num_cells / 2, run_args.num_cells / 2)
simulation.setup("box")

# create fibril
fibril_temp_file = os.path.join(output_folder, '{:d}_fibril.dat'.format(seed))
fibril_edges = rods_tools.prepare_fibril(model, run_args.seed_size, run_args.seed_phi,
                                         run_args.seed_theta, run_args.seed_r0, fibril_temp_file)
simulation.create_rods(state_ID=model.num_states-1, file=[fibril_temp_file])
os.remove(fibril_temp_file)
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
simulation.set_rod_dynamics("nve")

py_lmp.neigh_modify("every 1 delay 1")
py_lmp.timestep(run_args.dt)

# RANDOMISE INITIAL CONFIGURATION
simulation.deactivate_state(0, vx_eps=5.0)
py_lmp.command('run 10000')
simulation.activate_state(0)
py_lmp.reset_timestep(0)

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
py_lmp.thermo_style("custom", "step atoms", "pe temp")
dump_elems = "id x y z type mol"
try:
    dump_elems += " c_"+fibril_compute
except:
    pass
py_lmp.dump("dump_cmd", "all", "custom", out_freq, dump_path, dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
py_lmp.thermo(out_freq)

# RUN...
if model.num_states == 1 or run_args.mc_moves == 0:
    raise Exception("Multiple states need to exist and MC moves need to be made for fibrils to grow!")
mc_moves_per_run = int(run_args.mc_moves * simulation.rods_count())

py_lmp.command('run {:d} post no'.format(run_args.run_length-1)) #so output happens after state changes
remaining = args.simlen - run_args.run_length + 1
while True:
    success = simulation.state_change_MC(mc_moves_per_run)#, replenish=("box", 2*model.rod_radius, 10)) TODO
    if not args.silent:
        base_count = simulation.state_count(0)
        beta_count = simulation.state_count(1)
        print 'step {:d} / {:d} :  beta-to-soluble ratio = {:d}/{:d} = {:.5f} (accept rate = {:.5f})'.format(
              (i+1)*run_args.run_length, args.simlen, beta_count, base_count, float(beta_count)/base_count,
              float(success)/mc_moves_per_run)
    
    if remaining / run_args.run_length > 0:
        py_lmp.command('run {:d} post no'.format(run_args.run_length))
        remaining -= run_args.run_length
    else:
        py_lmp.command('run {:d} post no'.format(remaining))
        break
