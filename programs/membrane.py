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

import os
import argparse

parser = argparse.ArgumentParser(description=
                                 'Program for NVE+Langevin hybrid LAMMPS simulation of spherocylinder-like\
rods, using the "lammps_multistate_rods" library, above a 3-bead bilayer membrane.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('cfg_file',
                    help='path to the "lammps_multistate_rods" model configuration file')
parser.add_argument('run_file',
                    help='path to the run configuration file')

parser.add_argument('--seed', type=int,
                    help='the seed for random number generators')
parser.add_argument('--out', type=str, default=None,
                    help='name/path for the output folder (defaults to cfg_file path w/o ext)')

parser.add_argument('-o', '--output_freq', type=int,
                    help='configuration output frequency (in MD steps);\
default behavior is after every batch of MC moves')
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
    main_group = 'membrane'
    
    def __init__(self, sigma, wc, eps, Nx, Ny, zmax, above, below,
                 bead_mass=1.0, type_offset=0, bond_offset=0, angle_offset=0):
        '''
        sigma: the characteristic length of the membrane particles (diameter of beads)
        wc: the tail-tail interaction (cosine/squared) range
        eps: the depth (strength) of the interaction between lipid particles.
        Nx: length of the membrane in lipids
        Ny: width of the membrane in lipids
        zmax: the upper bound on the membrane area (including "above" space)
        above: the empty space buffer below zmax and the top of the upper layer
        below: the empty space buffer below the end of the bottom layer
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
        
        self.side_length = 1.1*sigma
        self.xmin = -Nx*self.side_length/2
        self.xmax = Nx*self.side_length/2
        self.ymin = -Ny*self.side_length/2
        self.ymax = Ny*self.side_length/2
        self.zmin = zmax - above*sigma - 6*self.side_length - below*sigma
        self.zmid = zmax - above*sigma - 3*self.side_length
        
        self.type_offset = type_offset
        self.bond_offset = bond_offset
        self.angle_offset = angle_offset
    
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
            
    def create_membrane(self, py_lmp, seed, append=True):
        '''
        Creates the membrane using the "read_data" command. If not already done, the
        membrane data for this object will be generated.
        
        seed: a random seed (for the temp membrane.dat file)
        
        append: if False this will create the simulation box through the "read_data",
        otherwise it will just add the membrane particles in the already existent box.
        '''
        if self.atom_id == 0:
            self.generate_data()
        temp_dat_file = str(seed)+'-temp-membrane.dat'
        self.output_data(temp_dat_file)
        
        all_atom_mol_ids = py_lmp.lmp.gather_atoms_concat("molecule", 0, 1)
        if len(all_atom_mol_ids) > 0:
            mol_offset = max(all_atom_mol_ids)
        else:
            mol_offset = 0
        
        read_data_args = []
        if append:
            read_data_args.append("add append")
#             read_data_args.append("offset {:d} {:d} {:d} 0 0".format(
#                 self.type_offset, self.bond_offset, self.angle_offset))
        py_lmp.read_data('"'+temp_dat_file+'"', ' '.join(read_data_args))
        
        py_lmp.group(Membrane.top_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+1, mol_offset+membrane.mol_id))
        py_lmp.group(Membrane.bottom_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+2, mol_offset+membrane.mol_id))
        py_lmp.group(Membrane.main_group, "type", *membrane.bead_types)
        
        os.remove(temp_dat_file)

# =======================================================================================

#from mpi4py import MPI #TODO make MPI work...
from lammps import PyLammps
import lammps_multistate_rods as rods

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
run_filename = os.path.splitext(os.path.basename(args.run_file))[0]
sim_ID = '{:s}_{:d}'.format(run_filename, seed)
    
dump_filename = sim_ID+'.dump'
dump_path = os.path.join(output_folder, dump_filename)

log_filename = '{:d}.lammps'.format(seed)
log_path = os.path.join(output_folder, log_filename)

run_args = rods.rod_model.Params()
execfile(args.run_file, {'__builtins__': None}, vars(run_args))

out_freq = args.output_freq if args.output_freq != None else run_args.run_length

model = rods.Rod_model(args.cfg_file)
model.generate_mol_files(output_folder)
type_offset = max(model.all_bead_types)
bond_offset = 1
membrane = Membrane(run_args.mem_sigma, run_args.mem_wc, run_args.mem_eps, 
                    run_args.mem_Nx, run_args.mem_Ny, 0.0, model.rod_length/2, 5.0,
                    type_offset=type_offset, bond_offset=bond_offset)
max_type = max(membrane.bead_types)
zmax = run_args.num_cells*run_args.cell_size

# ===== LAMMPS setup ====================================================================
py_lmp = PyLammps(cmdargs=['-screen','none'])
py_lmp.log('"'+log_path+'"')
py_lmp.units('lj')
py_lmp.dimension(3)
py_lmp.boundary('p p f')
py_lmp.atom_style('molecular')
py_lmp.angle_style('harmonic')
py_lmp.bond_style('hybrid fene zero')
py_lmp.pair_style('cosine/squared', model.global_cutoff)

py_lmp.region('box', 'block',
              membrane.xmin, membrane.xmax,
              membrane.ymin, membrane.ymax,
              0.0, zmax)
py_lmp.create_box(max_type, 'box',
                  'bond/types', 2, 'extra/bond/per/atom', 2,
                  'angle/types', 1, 'extra/angle/per/atom', 1,
                  'extra/special/per/atom', 6)

zwalls_fix = 'zwalls'
py_lmp.fix(zwalls_fix, 'all', 'wall/lj126',
           'zlo EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6),
           'zhi EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6))

for state_name in model.rod_states:
    py_lmp.molecule(state_name, '"'+os.path.join(output_folder, state_name+'.mol')+'"')

all_type_range = '{:d}*{:d}'.format(1, max_type)
py_lmp.mass(all_type_range, model.rod_mass*10**-10)
for bead_type in model.body_bead_types:
    py_lmp.mass(bead_type, model.rod_mass/model.num_beads[0])
for bead_type in membrane.bead_types:
    py_lmp.mass(bead_type, membrane.bead_mass)

py_lmp.angle_coeff(membrane.angle_type, 5.0, 180)
py_lmp.bond_coeff(1, 'zero')
py_lmp.bond_coeff(membrane.bond_type, 'fene',
                  30.0, 1.5*membrane.sigma, membrane.eps, membrane.sigma)
py_lmp.special_bonds('lj', 0.0, 1.0, 1.0) #not really necessary because of "neigh_modify exclude molecule/intra"

# set interactions (initially to 0.0, of whatever interaction, between all pairs of types)
py_lmp.pair_coeff(all_type_range, all_type_range, 0.0, model.rod_radius, model.rod_radius, "wca")
lj_factor = pow(2,1./6)
max_range = 0.0
# protein-protein interaction
for bead_types, (eps_val, int_type_key) in model.eps.iteritems():
    sigma = 0
    for bead_type in bead_types:
        if bead_type in model.body_bead_types:
            sigma += model.rod_radius
        else:
            for k in range(model.num_patches):
                if bead_type in model.patch_bead_types[k]:
                    sigma += model.patch_bead_radii[k]
                    break
    type_1 = bead_types[0]
    type_2 = bead_types[1]
    int_type = model.int_types[int_type_key]
    max_range = int_type[1] if int_type[1] > max_range else max_range
    py_lmp.pair_coeff(type_1, type_2, eps_val, sigma, sigma+int_type[1], "wca")
    
# head-head & head-tail inter-lipid interaction
for bead_type in membrane.bead_types:
    py_lmp.pair_coeff(membrane.head_type, bead_type, membrane.eps,
                      0.95*membrane.sigma*lj_factor,
                      0.95*membrane.sigma*lj_factor, 'wca')
# tail-tail inter-lipid interaction
for i in range(len(membrane.tail_types)):
    for j in range(i, len(membrane.tail_types)):
        py_lmp.pair_coeff(membrane.tail_types[i], membrane.tail_types[j], membrane.eps,
                          1.0*membrane.sigma*lj_factor,
                          (1.0*lj_factor + membrane.wc)*membrane.sigma, 'wca')

# lipid-protein interaction
sol_lipid_eps = run_args.mem_int_eps
sol_lipid_contact = 0.5*membrane.sigma*pow(2,1./6) + model.rod_radius 
sol_lipid_cutoff = sol_lipid_contact + run_args.mem_int_range
sol_body_types = filter(lambda t: t not in model.active_bead_types,
                       model.state_bead_types[0])
sol_tip_types = filter(lambda t: t in model.active_bead_types,
                      model.state_bead_types[0])
# lipid-protein volume-exclusion
for sol_body_type in sol_body_types:
    for mem_bead_type in membrane.bead_types:
        py_lmp.pair_coeff(sol_body_type, mem_bead_type, sol_lipid_eps,
                          sol_lipid_contact, sol_lipid_contact, 'wca')
# lipid-protein tip interaction
int_factors = (1.0, 0.5, 0.25)
for sol_tip_type in sol_tip_types:
    for mem_bead_type, k in zip(membrane.bead_types, int_factors):
        py_lmp.pair_coeff(sol_tip_type, mem_bead_type, k*sol_lipid_eps,
                          sol_lipid_contact, sol_lipid_cutoff, 'wca')
    
# ===== RODS ============================================================================
# GROUPS & COMPUTES
rods_group = 'rods'
cluster_group = 'active_rod_beads'
py_lmp.group(rods_group, 'empty')
py_lmp.variable('active', 'atom', '"' + 
                ' || '.join(['(type == {:d})'.format(t)
                             for t in model.active_bead_types])
                + '"')
py_lmp.group(cluster_group, 'dynamic', rods_group,
             'var', 'active', 'every', out_freq)
cluster_compute = rods.Simulation.cluster_compute
py_lmp.compute(cluster_compute, cluster_group, 'aggregate/atom', 2*model.rod_radius + max_range)

# FIXES & DYNAMICS
thermo_fix = 'thermostat'
py_lmp.fix(thermo_fix, 'all', 'langevin', run_args.temp, run_args.temp, run_args.damp, seed)#, "zero yes")
rod_dyn_fix = 'rod_dynamics'
py_lmp.fix(rod_dyn_fix, rods_group, 'rigid/nve/small', 'molecule', 'mol', model.rod_states[0])
py_lmp.neigh_modify("exclude", "molecule/intra", rods_group)
#py_lmp.fix_modify(rod_dyn_fix, 'dynamic/dof yes') 
py_lmp.compute_modify("thermo_temp", "dynamic/dof yes")

py_lmp.region('gcmc_init', 'block',
              membrane.xmin + model.rod_length/2, membrane.xmax - model.rod_length/2,
              membrane.ymin + model.rod_length/2, membrane.ymax - model.rod_length/2,
              0.0 + model.rod_length/2, zmax - model.rod_length/2)
rod_gcmc_fix = 'rod_gcmc'
py_lmp.fix(rod_gcmc_fix, rods_group, 'gcmc', 100, 200, 0,
           0, seed, run_args.temp, run_args.mu, 0.0,
           'region', 'gcmc_init', 'mol', model.rod_states[0],
           'rigid', rod_dyn_fix, 'tfac_insert', 100)

# TEST DUMP...
py_lmp.thermo_style('custom', 'step atoms', 'pe temp')
py_lmp.thermo(out_freq)
py_lmp.dump('test_dump', 'all', 'custom', out_freq, dump_path+'_test',
            'id x y z type mol c_'+cluster_compute)
py_lmp.dump_modify('test_dump', 'sort id')

# GENERATING INITIAL CONFIGURATION
py_lmp.neigh_modify('every', 1, 'delay', 1)
py_lmp.timestep(run_args.dt)
py_lmp.command('run 2000')
py_lmp.unfix(rod_gcmc_fix)
py_lmp.reset_timestep(0)

# ===== MEMBRANE ========================================================================

# create membrane (box update, create membrane & groups, ...)
membrane.create_membrane(py_lmp, seed, append=True)

py_lmp.fix(zwalls_fix, 'all', 'wall/lj126',
           'zlo EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6),
           'zhi EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6))
    
# GROUPS & COMPUTES
adsorbed_group = 'mem_active'
py_lmp.variable('mem_active', 'atom', '"' + 
                ' || '.join(['(type == {:d})'.format(t)
                             for t in list(model.active_bead_types) + membrane.bead_types])
                + '"')
py_lmp.group(adsorbed_group, 'dynamic', 'all',
             'var', 'mem_active', 'every', out_freq)
adsorbed_compute = 'mem_cluster'
py_lmp.compute(adsorbed_compute, adsorbed_group, 'aggregate/atom', sol_lipid_cutoff)
top_msd_compute = 'mem_top_msd'
py_lmp.compute(top_msd_compute, membrane.top_layer_group, 'msd')
bottom_msd_compute = 'mem_bottom_msd'
py_lmp.compute(bottom_msd_compute, membrane.bottom_layer_group, 'msd')
full_msd_compute = 'mem_full_msd'
py_lmp.compute(full_msd_compute, membrane.main_group, 'msd')

# FIXES & DYNAMICS
mem_dyn_fix = 'mem_dynamics'
py_lmp.fix(mem_dyn_fix, membrane.main_group, 'nph',
           'x 0.0 0.0 10', 'y 0.0 0.0 10', 'couple xy',
           'dilate', membrane.main_group)
py_lmp.neigh_modify('exclude', 'molecule/intra', membrane.main_group)

# MEMBRANE EQUILIBRATION
py_lmp.command('run 1000')
py_lmp.reset_timestep(0)

# ===== FINAL ===========================================================================
py_lmp.variable('temp_x', 'equal', '"xhi - {:f}"'.format(model.rod_length/2))
py_lmp.variable('temp_y', 'equal', '"yhi - {:f}"'.format(model.rod_length/2))
py_lmp.region('gcmc_box', 'block',
              '-${temp_x}', '${temp_x}',
              '-${temp_y}', '${temp_y}',
              0.0 + 2*model.rod_length, zmax - model.rod_length/2)
py_lmp.fix(rod_gcmc_fix, rods_group, 'gcmc', 1000, 10, 0,
           0, seed, run_args.temp, run_args.mu, 0.0,
           'region', 'gcmc_box', 'mol', model.rod_states[0],
           'rigid', rod_dyn_fix, 'tfac_insert', 100)

# OUTPUT
py_lmp.variable("area", "equal", "lx*ly")
py_lmp.thermo_style("custom", "step atoms", "pe temp", "press lx ly v_area")
py_lmp.thermo(out_freq)
dump_elems = 'id x y z type mol c_{:s} c_{:s}'.format(cluster_compute, adsorbed_compute)
py_lmp.dump("dump_cmd", "all", "custom", out_freq, dump_path, dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
# MSD data (only in x & y directions)
py_lmp.variable("top_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(top_msd_compute))
py_lmp.fix("mem_top_msd", "all", "ave/time", 1, 1, out_freq, "v_top_msd",
           "file", os.path.join(output_folder, sim_ID+'_mem_top.msd'))
py_lmp.variable("bottom_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(bottom_msd_compute))
py_lmp.fix("mem_bottom_msd", "all", "ave/time", 1, 1, out_freq, "v_bottom_msd",
           "file", os.path.join(output_folder, sim_ID+'_mem_bottom.msd'))
py_lmp.variable("full_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(full_msd_compute))
py_lmp.fix("mem_full_msd", "all", "ave/time", 1, 1, out_freq, "v_full_msd",
           "file", os.path.join(output_folder, sim_ID+'_mem_full.msd'))

# RUN...
py_lmp.command('run {:d}'.format(run_args.sim_length))
