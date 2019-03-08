# encoding: utf-8
'''
This application does an NVE+Langevin LAMMPS simulation of spherocylinder-like rods
(defined in a .cfg file) using the "lammps_multistate_rods" library, along with an
NPH+Langevin simulation of a 3-bead bilayer membrane (Cooke et al. 2008).
The initial locations of the rods are at SC lattice points, defined by the input
params, above the membrane, and their orientations are randomly determined at each
insertion point.

Created on 16 Mar 2018

@author: Eugen Rožić
'''

import argparse

parser = argparse.ArgumentParser(description=
                                 'Program for NVE+Langevin hybrid LAMMPS simulation of\
spherocylinder-like rods above a 3-bead bilayer membrane.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('config_file',
                    help='path to the "lammps_multistate_rods" model config file')
parser.add_argument('output_folder',
                    help='name for the folder that will be created for output files')
parser.add_argument('cell_size', type=float,
                    help='size of an SC cell (i.e. room for one rod)')
parser.add_argument('num_cells', type=float,
                    help='the number of cells in the z dimension')
parser.add_argument('sim_length', type=int,
                    help='the total number of MD steps to simulate')

parser.add_argument('--seed', type=int,
                    help='the seed for random number generators')

parser.add_argument('--Nx', type=int, default=30,
                    help='length of membrane in lipids')
parser.add_argument('--Ny', type=int, default=30,
                    help='width of membrane in lipids')
parser.add_argument('--mem_sigma', type=float, default=1.0,
                    help='the lipid characteristic length (bead diameter)')
parser.add_argument('--mem_wc', type=float, default=1.4,
                    help='the width of the membrane potential')
parser.add_argument('--mem_eps', type=float, default=1.0,
                    help='the strength of the membrane potential')
parser.add_argument('--int_eps', type=float, default=5.0,
                    help='the strength of the membrane-lipid interaction')

parser.add_argument('-T', '--temp', default=1.0, type=float,
                    help='the temperature of the system (e.g. for Langevin)')
parser.add_argument('-D', '--damp', default=0.1, type=float,
                    help='viscous damping (for Langevin)')

parser.add_argument('-t', '--timestep', default=0.01, type=float,
                    help='timestep length (in lj units)')
parser.add_argument('-R', '--runlen', default=200, type=int,
                    help='number of MD steps between MC moves')
parser.add_argument('--MC_moves', default=1.0, type=float,
                    help='number of MC moves per rod between MD runs')

parser.add_argument('--clusters', default=1.0, type=float,
                    help='the max distance (in rod radii) for two rods to be\
in the same cluster (put to 0.0 to turn cluster tracking off)')

parser.add_argument('-o', '--output_freq', type=int,
                    help='configuration output frequency (in MD steps);\
default behavior is after every batch of MC moves')
parser.add_argument('-s', '--silent', action='store_true',
                    help="doesn't print anything to stdout")

args = parser.parse_args()

#========================================================================================

import numpy as np

class Membrane(object):
    '''
    TODO
    '''
    
    def __init__(self, Nx, Ny, zmax, below, sigma, wc, eps):
        self.Nx = Nx
        self.Ny = Ny
        self.zmax = zmax
        self.below = below
        self.sigma = sigma
        self.wc = wc
        self.eps = eps
        
        self.side_length = 1.1*sigma
        self.xmin = -Nx*self.side_length/2
        self.xmax = Nx*self.side_length/2
        self.ymin = -Ny*self.side_length/2
        self.ymax = Ny*self.side_length/2
        self.zmin = zmax - 6*self.side_length - below*sigma
    
        self.bead_mass = 1.0
        self.bead_types = [1,2,3]
        self.bond_type = 1
        self.angle_type = 1
    
        self.atoms = []
        self.atom_id = 0
        self.mol_id = 0
        self.bonds = []
        self.bond_id = 0
        self.angles = []
        self.angle_id = 0
    
    #TODO LAMMPS setup etc.
    
    def _add_lipid(self, x, y, reverse=False):
        self.mol_id += 1
        for i in range(3):
            self.atom_id += 1
            if reverse:
                z = self.zmax - (5.5 - i)*self.side_length
            else:
                z = self.zmax - (i + 0.5)*self.side_length
            self.atoms.append("\t{:d} {:d} {:d} {:f} {:f} {:f}".format(
                self.atom_id, self.mol_id, self.bead_types[i], x, y, z))
        self.bond_id += 1
        self.bonds.append("\t{:d} {:d} {:d} {:d}".format(
            self.bond_id, self.bond_type, self.atom_id-2, self.atom_id-1))
        self.bond_id += 1
        self.bonds.append("\t{:d} {:d} {:d} {:d}".format(
            self.bond_id, self.bond_type, self.atom_id-1, self.atom_id))
        self.angle_id += 1
        self.angles.append("\t{:d} {:d} {:d} {:d} {:d}".format(
            self.angle_id, self.angle_type, self.atom_id-2, self.atom_id-1, self.atom_id))
        
    def generate(self):
        '''
        Fills the internal variables of the object and makes it ready to output to file.
        '''
        for x in np.arange(self.xmin+self.side_length/2, self.xmax, self.side_length):
            for y in np.arange(self.ymin+self.side_length/2, self.ymax, self.side_length):
                self._add_lipid(x, y) # upper layer lipid (1-2-3)
                self._add_lipid(x, y, reverse=True) # lower layer lipid (3-2-1)
    
    def clear(self):
        '''
        Clears the internal variables filled with "generate_membrane"
        '''
        self.atoms = []; self.atom_id = 0
        self.mol_id = 0
        self.bonds = []; self.bond_id = 0
        self.angles = []; self.angle_id = 0
    
    def output(self, filepath):
        '''
        Outputs the generated membrane to a file at the given path.
        '''
        header = ["LAMMPS Description\n",
                  "\t{:d} atoms".format(len(self.atoms)),
                  "\t{:d} bonds".format(len(self.bonds)),
                  "\t{:d} angles\n".format(len(self.angles)),
                  "\t{:d} atom types".format(len(self.bead_types)),
                  "\t1 bond types",
                  "\t1 angle types\n",
                  "\t{:.5f} {:.5f} xlo xhi".format(self.xmin, self.xmax),
                  "\t{:.5f} {:.5f} ylo yhi".format(self.ymin, self.ymax),
                  "\t{:.5f} {:.5f} zlo zhi\n".format(self.zmin, self.zmax)]
        
        with open(filepath, 'w') as dat_file:
            dat_file.write('\n'.join(header))
            dat_file.write("\nMasses\n\n")
            for bead_type in self.bead_types:
                dat_file.write("\t{:d} {:f}\n".format(bead_type, self.bead_mass))
            dat_file.write("\nAtoms\n\n")
            dat_file.write('\n'.join(self.atoms))
            dat_file.write("\n\nBonds\n\n")
            dat_file.write('\n'.join(self.bonds))
            dat_file.write("\n\nAngles\n\n")
            dat_file.write('\n'.join(self.angles))

#========================================================================================

import os
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

#from mpi4py import MPI #TODO make MPI work...
from lammps import PyLammps
import lammps_multistate_rods as rods

if args.seed is None:
    import time
    args.seed = int((time.time() % 1)*1000000)
    print "WARNING: no seed given explicitly; using:", args.seed
    
dump_path = str(args.mem_eps)+'-'+str(args.int_eps)+'_'+str(args.seed)+'.dump'
dump_path = os.path.join(args.output_folder, dump_path)
log_path = os.path.join(args.output_folder, str(args.seed)+'_lammps.log')
mem_dat_path = os.path.join(args.output_folder, str(args.seed)+'_membrane.dat')

py_lmp = PyLammps(cmdargs=['-screen','none'])
py_lmp.log('"'+log_path+'"')
model = rods.Model(args.config_file)
simulation = rods.Simulation(py_lmp, model, args.seed, args.output_folder,
                             clusters=args.clusters)

# various vars...
lj_factor = pow(2,1./6)

# LAMMPS setup
py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p f")
py_lmp.lattice("sc", 1/(args.cell_size**3))
x_cells = args.Nx*args.mem_sigma/args.cell_size
y_cells = args.Nx*args.mem_sigma/args.cell_size
z_cells = args.num_cells
py_lmp.region("rod_box", "block",
                  -x_cells/2, x_cells/2,
                  -y_cells/2, y_cells/2,
                   0, z_cells)

# setup for rods (+ everything the simulations need for the membrane)
simulation.setup("rod_box", atom_style="molecular", type_offset=3,
                 extra_pair_styles=[], overlay=False,
                 bond_offset=1, extra_bond_styles=['fene'], everything_else = [
                 "extra/bond/per/atom", 2,
                 "angle/types", 1, "extra/angle/per/atom", 1,
                 "extra/special/per/atom", 3])
# setup for membrane particles
py_lmp.angle_style("harmonic")
py_lmp.angle_coeff(1, 5.0, 180)
py_lmp.bond_coeff(1, "fene", 30.0, 1.5*args.mem_sigma, args.mem_eps, args.mem_sigma)
py_lmp.special_bonds("lj", 0.0, 1.0, 1.0)
py_lmp.pair_coeff(1, '1*3', "cosine/squared", args.mem_eps, 0.95*args.mem_sigma*lj_factor,
                  0.95*args.mem_sigma*lj_factor, "wca") #head-head & head-tail
py_lmp.pair_coeff('2*3', '2*3', "cosine/squared", args.mem_eps, 1.0*args.mem_sigma*lj_factor,
                  (1.0*lj_factor + args.mem_wc)*args.mem_sigma, "wca") #tail-tail
# membrane-rod interaction...
sol_lipid_eps = args.int_eps
sol_lipid_contact = 0.5*args.mem_sigma + model.rod_radius 
sol_lipid_cutoff = sol_lipid_contact + 0.5*args.mem_sigma
sol_body_type = model.state_structures[0][0][0] + simulation.type_offset
sol_tip_type = model.state_structures[0][0][-1] + simulation.type_offset
py_lmp.pair_coeff('1*3', sol_body_type, "lj/cut", sol_lipid_eps,
                  sol_lipid_contact, sol_lipid_contact) #this is just volume-exclusion
py_lmp.pair_coeff(1, sol_tip_type, "lj/cut", sol_lipid_eps,
                  sol_lipid_contact, sol_lipid_cutoff)
py_lmp.pair_coeff(2, sol_tip_type, "lj/cut", 0.5*sol_lipid_eps,
                  sol_lipid_contact, sol_lipid_cutoff)
py_lmp.pair_coeff(3, sol_tip_type, "lj/cut", 0.25*sol_lipid_eps,
                  sol_lipid_contact, sol_lipid_cutoff)

# create membrane particles
mem_zmax = -args.cell_size/2
mem_below = 3.0
membrane = Membrane(args.Nx, args.Ny, mem_zmax, mem_below,
                    args.mem_sigma, args.mem_wc, args.mem_eps)
membrane.generate()
membrane.output(mem_dat_path)
py_lmp.read_data('"'+mem_dat_path+'"', "add append")
py_lmp.group("membrane", "type", *range(1, simulation.type_offset+1))
# create rods
py_lmp.region("rod_init", "block",
              membrane.xmin, membrane.xmin + (x_cells-0.01)*args.cell_size,
              membrane.ymin, membrane.ymin + (y_cells-0.01)*args.cell_size,
              0.0, z_cells*args.cell_size, "units box")
simulation.create_rods(region = "rod_init")

# DYNAMICS
# -> all particles on the same temparature (langevin thermostat)
py_lmp.fix("thermostat", "all", "langevin", args.temp, args.temp, args.damp, args.seed)#, "zero yes")
# -> nve for the rods (has to come before membrane nph - WHY?!?!)
simulation.set_rod_dynamics("nve")
py_lmp.fix("wall", "all", "wall/reflect", "zhi EDGE")
# -> nph for the membrane particles
# --> with intra-lipid pair-interaction exclusion (bonds take care of everything)
py_lmp.fix("mem_dyn", "membrane", "nph", "x 0.0 0.0 10", "y 0.0 0.0 10", "couple xy", "dilate membrane")
py_lmp.neigh_modify("exclude", "molecule/intra", "membrane")

py_lmp.neigh_modify("every 1 delay 1")

# OUTPUT
py_lmp.variable("area", "equal", "lx*ly")
py_lmp.thermo_style("custom", "step atoms", "pe temp", "press lx ly v_area")
dump_elems = "id x y z type mol"
if args.clusters > 0.0:
    dump_elems += " c_"+simulation.cluster_compute
if (args.output_freq != None):
    py_lmp.dump("dump_cmd", "all", "custom", args.output_freq, dump_path, dump_elems)
    py_lmp.dump_modify("dump_cmd", "sort id")
    py_lmp.thermo(args.output_freq)
else:
    py_lmp.variable("out_timesteps", "equal", "stride(1,{:d},{:d})".format(args.sim_length+1, args.runlen))
    py_lmp.dump("dump_cmd", "all", "custom", 1, dump_path, dump_elems)
    py_lmp.dump_modify("dump_cmd", "every v_out_timesteps", "sort id")
    py_lmp.thermo(args.runlen)

# RUN...
mc_moves_per_run = 0
if model.num_states > 1:
    mc_moves_per_run = int(args.MC_moves * simulation.rods_count())
    
py_lmp.timestep(args.timestep)

if mc_moves_per_run == 0:
    py_lmp.command('run {:d}'.format(args.sim_length))
else:
    for i in range(int(args.sim_length/args.runlen)-1):   
        py_lmp.command('run {:d} post no'.format(args.runlen))
        success = simulation.state_change_MC(mc_moves_per_run)
        if not args.silent:
            base_count = simulation.state_count(0)
            beta_count = simulation.state_count(1)
            print 'step {:d} / {:d} :  beta-to-soluble ratio = {:d}/{:d} = {:.5f} (accept rate = {:.5f})'.format(
                    (i+1)*args.runlen, args.sim_length, beta_count, base_count,
                        float(beta_count)/base_count, float(success)/mc_moves_per_run)
            
    py_lmp.command('run {:d} post no'.format(args.runlen))
