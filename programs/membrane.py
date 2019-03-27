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

parser.add_argument('--clusters', default=3.0, type=float,
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
    This class holds information about a bilayer membrane made of 3-bead lipid models
    (Cooke at al. 2005).
    It has methods to generate a data file that can be read in by LAMMPS via "read_data"
    to setup a bilayer membrane simulation, and methods to setup LAMMPS with proper
    commands for this kind of simulation.
    All interaction between membrane/lipid particles and other particles in the simulation
    have to specified elsewhere (or manually), but that can be done with the help of
    information held by an instance of this class. 
    '''
    
    top_layer_group = 'mem_top_layer'
    bottom_layer_group = 'mem_bottom_layer'
    membrane_group = 'membrane'
    
    def __init__(self, Nx, Ny, zmax, above, below, sigma, wc, eps,
                 bead_mass=1.0, type_offset=0, bond_offset=0, angle_offset=0):
        '''
        Nx: length of the membrane in lipids
        Ny: width of the membrane in lipids
        zmax: the upper bound on the membrane area (including "above" space)
        above: the empty space buffer below zmax and the top of the upper layer
        below: the empty space buffer below the end of the bottom layer
        sigma: the characteristic length of the membrane particles (diameter of beads)
        wc: the tail-tail interaction (cosine/squared) range
        eps: the depth (strength) of the interaction between lipid particles.
        '''
        self.Nx = Nx
        self.Ny = Ny
        self.zmax = zmax
        self.above = above
        self.below = below
        self.sigma = sigma
        self.wc = wc
        self.eps = eps
        self.bead_mass = bead_mass
        self.type_offset = type_offset
        self.bond_offset = bond_offset
        self.angle_offset = angle_offset
        
        self.side_length = 1.1*sigma
        self.xmin = -Nx*self.side_length/2
        self.xmax = Nx*self.side_length/2
        self.ymin = -Ny*self.side_length/2
        self.ymax = Ny*self.side_length/2
        self.zmin = zmax - above*sigma - 6*self.side_length - below*sigma
        self.zmid = zmax - above*sigma - 3*self.side_length
    
        self.bead_types = [i + type_offset for i in (1,2,3)]
        self.head_type = self.bead_types[0]
        self.tail_types = self.bead_types[1:]
        self.bond_type = 1 + bond_offset
        self.angle_type = 1 + angle_offset
    
        self.atoms = []
        self.atom_id = 0
        self.mol_id = 0
        self.bonds = []
        self.bond_id = 0
        self.angles = []
        self.angle_id = 0
    
    def lammps_setup(self, py_lmp, angle_coeff=5.0, fene_coeff=30.0, full=False):
        '''
        Sets up the LAMMPS parameters relevant for membrane simulation, i.e.
            - harmonic angle bond in the lipid (strength 1/2 * 10 = 5)
            - fene bond between beads of the lipid (30.0, r = 1.5*R)
            - "cosine/squared" pair style coefficients between membrane particles
        
        full: if True also sets up atom, bond and pair styles (minimum necessary for
        simulating the membrane)
        '''
        self.py_lmp = py_lmp
        if full:
            py_lmp.atom_style("molecular")
            py_lmp.pair_style('cosine/squared', 3*self.sigma)
            py_lmp.bond_style('fene')
        py_lmp.angle_style("harmonic")
        py_lmp.angle_coeff(self.angle_type, angle_coeff, 180)
        py_lmp.bond_coeff(self.bond_type, "fene", fene_coeff, 1.5*self.sigma, self.eps, self.sigma)
        #py_lmp.special_bonds("lj", 0.0, 1.0, 1.0)
        
        lj_factor = pow(2,1./6)
        # head-head & head-tail interaction
        for bead_type in self.bead_types:
            py_lmp.pair_coeff(self.head_type, bead_type, "cosine/squared", self.eps,
                              0.95*self.sigma*lj_factor, 0.95*self.sigma*lj_factor, "wca")
        # tail-tail interaction
        for i in range(len(self.tail_types)):
            for j in range(i, len(self.tail_types)):
                py_lmp.pair_coeff(self.tail_types[i], self.tail_types[j], "cosine/squared",
                                  self.eps, 1.0*self.sigma*lj_factor,
                                  (1.0*lj_factor + self.wc)*self.sigma, "wca")
    
    def _add_lipid(self, x, y, reverse=False):
        self.mol_id += 1
        for i in range(3):
            self.atom_id += 1
            if reverse:
                z = self.zmid - (2.5 - i)*self.side_length
            else:
                z = self.zmid + (2.5 - i)*self.side_length
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
        
    def generate_data(self):
        '''
        Fills the internal variables of the object and makes it ready to output to file.
        '''
        for x in np.arange(self.xmin+self.side_length/2, self.xmax, self.side_length):
            for y in np.arange(self.ymin+self.side_length/2, self.ymax, self.side_length):
                self._add_lipid(x, y) # upper layer lipid (1-2-3)
                self._add_lipid(x, y, reverse=True) # lower layer lipid (3-2-1)
    
    def clear_data(self):
        '''
        Clears the internal variables filled with "generate_data"
        '''
        self.atoms = []; self.atom_id = 0
        self.mol_id = 0
        self.bonds = []; self.bond_id = 0
        self.angles = []; self.angle_id = 0
    
    def output_data(self, filepath):
        '''
        Outputs the generated membrane data to a file at the given path.
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
            
    def create_membrane(self, append=True):
        '''
        Creates the membrane using the "read_data" command. If not already done, the
        membrane data for this object will be generated.
        
        append: if False this will create the simulation box through the "read_data",
        otherwise it will just add the membrane particles in the already existent box.
        '''
        try:
            isinstance(py_lmp, object)
        except NameError:
            raise Exception('"Membrane.create_membrane" can only be called after "Membrane.lammps_setup"!')
        
        if self.atom_id == 0:
            self.generate_data()
        temp_dat_file = 'temp-membrane.dat'
        self.output_data(temp_dat_file)
        
        all_atom_mol_ids = self.py_lmp.lmp.gather_atoms_concat("molecule", 0, 1)
        if len(all_atom_mol_ids) > 0:
            mol_offset = max(all_atom_mol_ids)
        else:
            mol_offset = 0
        
        read_data_args = []
        if append:
            read_data_args.append("add append")
            read_data_args.append("offset {:d} {:d} {:d} 0 0".format(
                self.type_offset, self.bond_offset, self.angle_offset))
        self.py_lmp.read_data('"'+temp_dat_file+'"', ' '.join(read_data_args))
        
        self.py_lmp.group(Membrane.top_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+1, mol_offset+membrane.mol_id))
        self.py_lmp.group(Membrane.bottom_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+2, mol_offset+membrane.mol_id))
        self.py_lmp.group(Membrane.membrane_group, "type", *membrane.bead_types)
        
        os.remove(temp_dat_file)

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

out_freq = args.output_freq if args.output_freq != None else args.runlen

py_lmp = PyLammps(cmdargs=['-screen','none'])
py_lmp.log('"'+log_path+'"')
model = rods.Model(args.config_file)
simulation = rods.Simulation(py_lmp, model, args.seed, args.output_folder,
                             clusters=args.clusters)
membrane = Membrane(args.Nx, args.Ny, 0.0, model.rod_length/2, 5.0,
                    args.mem_sigma, args.mem_wc, args.mem_eps)
#membrane.generate()
#mem_dat_path = os.path.join(args.output_folder, str(args.seed)+'_membrane.dat')
#membrane.output(mem_dat_path)

# LAMMPS setup
py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.boundary("p p f")
py_lmp.lattice("sc", 1/(args.cell_size**3))

# setup for rods (+ everything the simulations need for the membrane)
x_cells = args.Nx*membrane.sigma/args.cell_size
y_cells = args.Nx*membrane.sigma/args.cell_size
z_cells = args.num_cells
py_lmp.region("rod_box", "block", -x_cells/2, x_cells/2,
                                  -y_cells/2, y_cells/2,
                                           0, z_cells)
simulation.setup("rod_box", atom_style="molecular", type_offset=max(membrane.bead_types),
                 extra_pair_styles=[], overlay=False,
                 bond_offset=membrane.bond_type, extra_bond_styles=['fene'],
                 everything_else = ["extra/bond/per/atom", 2,
                 "angle/types", membrane.angle_type, "extra/angle/per/atom", 1,
                 "extra/special/per/atom", 3])
membrane.lammps_setup(py_lmp)

# membrane-rod interaction...
sol_lipid_eps = args.int_eps
sol_lipid_contact = 0.5*membrane.sigma + model.rod_radius 
sol_lipid_cutoff = sol_lipid_contact + 0.5*membrane.sigma
sol_body_type = model.state_structures[0][0][0] + simulation.type_offset
sol_tip_type = model.state_structures[0][0][-1] + simulation.type_offset
# lipid-protein volume-exclusion
for mem_bead_type in membrane.bead_types:
    py_lmp.pair_coeff(mem_bead_type, sol_body_type, "lj/cut", sol_lipid_eps,
                      sol_lipid_contact, sol_lipid_contact)
# lipid-protein tip interaction
int_factors = (1.0, 0.5, 0.25)
for bead_type, k in zip(membrane.bead_types, int_factors):
    py_lmp.pair_coeff(bead_type, sol_tip_type, "lj/cut", k*sol_lipid_eps,
                      sol_lipid_contact, sol_lipid_cutoff)

# create membrane
membrane.create_membrane(append=True)
# create rods
py_lmp.region("rod_init", "block",
              membrane.xmin, membrane.xmin + (x_cells-0.01)*args.cell_size,
              membrane.ymin, membrane.ymin + (y_cells-0.01)*args.cell_size,
              0.0, z_cells*args.cell_size, "units box")
simulation.create_rods(region = "rod_init")

# a group and a compute that flag adsorbed proteins (their active beads)
py_lmp.group("mem+active", "union", simulation.active_beads_group, membrane.membrane_group)
py_lmp.compute("mem_cluster", "mem+active", "aggregate/atom", 0.9*sol_lipid_cutoff)
# this counts adsorbed proteins natively (some non-watertight assumptions, like cluster ID)
# from itertools import groupby
# # this is A (type, occurrences) pair for AN active body bead in the base (0-index) state
# an_active_sol_body_bead = [(body_bead_type, len(list(group)))
#                            for body_bead_type, group in groupby(model.state_structures[0][0])
#                            if body_bead_type in model.active_bead_types][0]
# py_lmp.variable("adsorbed", "atom", '"((c_mem_ads == 1) && (type == {:d}))"'.format(
#                     an_active_sol_body_bead[0] + simulation.type_offset))
# py_lmp.group("adsorbed", "dynamic", simulation.rods_group, "var", "adsorbed", "every", out_freq)
# py_lmp.variable("nads", "equal", "count(adsorbed)/{:d}".format(an_active_sol_body_bead[1]))

# MSD computes
py_lmp.compute("mem_top_msd", membrane.top_layer_group, "msd")
py_lmp.compute("mem_bottom_msd", membrane.bottom_layer_group, "msd")
py_lmp.compute("mem_full_msd", membrane.membrane_group, "msd")

# DYNAMICS
# -> all particles on the same temparature (langevin thermostat)
py_lmp.fix("thermostat", "all", "langevin", args.temp, args.temp, args.damp, args.seed)#, "zero yes")
# -> nve for the rods (has to come before membrane nph - WHY?!?!)
simulation.set_rod_dynamics("nve")
# -> nph for the membrane particles
# --> with intra-lipid pair-interaction exclusion (bonds take care of everything)
py_lmp.fix("mem_dyn", "membrane", "nph", "x 0.0 0.0 10", "y 0.0 0.0 10", "couple xy", "dilate membrane")
py_lmp.neigh_modify("exclude", "molecule/intra", "membrane")
# -> reflecting soft walls in z dimension (so no particles are lost)
py_lmp.fix("zwalls", "all", "wall/lj126",
           "zlo EDGE", membrane.eps, membrane.sigma, membrane.sigma*pow(2,1./6),
           "zhi EDGE", membrane.eps, membrane.sigma, membrane.sigma*pow(2,1./6))

py_lmp.neigh_modify("every 1 delay 1")

# OUTPUT
py_lmp.variable("area", "equal", "lx*ly")
py_lmp.thermo_style("custom", "step atoms", "pe temp", "press lx ly v_area")
dump_elems = "id x y z type mol c_mem_cluster"
if args.clusters > 0.0:
    dump_elems += " c_"+simulation.cluster_compute
if (args.output_freq != None):
    py_lmp.dump("dump_cmd", "all", "custom", args.output_freq, dump_path, dump_elems)
    py_lmp.dump_modify("dump_cmd", "sort id")
else:
    py_lmp.variable("out_timesteps", "equal", "stride(1,{:d},{:d})".format(args.sim_length+1, args.runlen))
    py_lmp.dump("dump_cmd", "all", "custom", 1, dump_path, dump_elems)
    py_lmp.dump_modify("dump_cmd", "every v_out_timesteps", "sort id")
py_lmp.fix("mem_top_msd", "all", "ave/time", 1, 1, out_freq, "c_mem_top_msd[4]",
           "file", os.path.join(args.output_folder, str(args.seed)+'mem_top.msd'))
py_lmp.fix("mem_bottom_msd", "all", "ave/time", 1, 1, out_freq, "c_mem_bottom_msd[4]",
           "file", os.path.join(args.output_folder, str(args.seed)+'mem_bottom.msd'))
py_lmp.fix("mem_full_msd", "all", "ave/time", 1, 1, out_freq, "c_mem_full_msd[4]",
           "file", os.path.join(args.output_folder, str(args.seed)+'mem_full.msd'))
# py_lmp.fix("nads_out", "all", "ave/time", 1, 1, out_freq, "v_nads",
#            "file", os.path.join(args.output_folder, 'adsorbed.dat'))
py_lmp.thermo(args.output_freq)

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
