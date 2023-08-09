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
from mpi4py import MPI

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()

import os
import argparse

parser = argparse.ArgumentParser(description='Program for NVE+Langevin hybrid LAMMPS'\
                                 ' simulation of spherocylinder-like rods, using the'\
                                 ' "lammps_multistate_rods" library, above a 3-bead bilayer membrane.',
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
        print("WARNING: no seed given explicitly; using:", seed)
    seed = mpi_comm.bcast(seed, root = 0)
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
                 bead_mass = 1.0, type_offset = 0, bond_offset = 0, angle_offset = 0):
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
    
    def _add_lipid(self, x, y, reverse = False):
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
            self.bond_id, self.bond_type, self.atom_id - 2, self.atom_id-1))
        self.bond_id += 1
        self.bonds.append("\t{:d} {:d} {:d} {:d}".format(
            self.bond_id, self.bond_type, self.atom_id - 1, self.atom_id))
        self.angle_id += 1
        self.angles.append("\t{:d} {:d} {:d} {:d} {:d}".format(
            self.angle_id, self.angle_type, self.atom_id - 2, self.atom_id-1, self.atom_id))
        
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
            
    def create_membrane(self, py_lmp, seed, append = True):
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
        if mpi_rank == 0:
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
        py_lmp.read_data('"' + temp_dat_file + '"', ' '.join(read_data_args))
        
        py_lmp.group(Membrane.top_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+1, mol_offset+membrane.mol_id))
        py_lmp.group(Membrane.bottom_layer_group, "molecule",
                     "{:d}:{:d}:2".format(mol_offset+2, mol_offset+membrane.mol_id))
        py_lmp.group(Membrane.main_group, "type", *membrane.bead_types)
        
        if mpi_rank == 0:
            os.remove(temp_dat_file)

# =======================================================================================
from lammps import PyLammps
import lammps_multistate_rods as rods

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
    with open(args.run_file) as f:
        compiled_file = compile(f.read(), args.run_file, 'exec')
    exec(compiled_file,
         {'__builtins__' : None, 'True' : True, 'False' : False, 'None' : None},
         vars(run_args))
run_args = mpi_comm.bcast(run_args, root = 0)

out_freq = args.output_freq if args.output_freq != None else run_args.mc_every

rod_params = rods.Rod_params()
if mpi_rank == 0:
    rod_params.from_file(args.cfg_file);
rod_params = mpi_comm.bcast(rod_params, root = 0)

model = rods.Rod_model(rod_params)
membrane = Membrane(run_args.mem_sigma, run_args.mem_wc, run_args.mem_eps, 
                    run_args.mem_Nx, run_args.mem_Ny, 0.0, model.rod_length/2, 5.0)

# ===== LAMMPS setup ====================================================================
if mpi_rank == 0:
    py_lmp = PyLammps(cmdargs = ['-echo', 'both'], comm = mpi_comm, verbose = args.verbose)
else:
    py_lmp = PyLammps(cmdargs = ['-echo', 'both'], comm = mpi_comm)
py_lmp.log('"' + log_path + '"')

simulation = rods.Simulation(py_lmp, model, run_args.temp, seed, output_folder)

py_lmp.units('lj')
py_lmp.dimension(3)
py_lmp.boundary('p p f')
py_lmp.region('box', 'block',
              membrane.xmin, membrane.xmax,
              membrane.ymin, membrane.ymax,
              0.0, run_args.Lz)

simulation.setup('box', type_offset = max(membrane.bead_types),
                 extra_pair_styles = [('cosine/squared', model.global_cutoff)],
                 bond_offset = membrane.bond_type, extra_bond_styles=['fene'],
                 opt = ['angle/types', 1, 'extra/angle/per/atom', 1])
py_lmp.angle_style('harmonic')

zwalls_fix = 'zwalls'
py_lmp.fix(zwalls_fix, 'all', 'wall/lj126',
           'zlo EDGE', 1.0, model.rod_radius, model.rod_radius * pow(2, 1./6),
           'zhi EDGE', 1.0, model.rod_radius, model.rod_radius * pow(2, 1./6))

for bead_type in membrane.bead_types:
    py_lmp.mass(bead_type, membrane.bead_mass)
py_lmp.angle_coeff(membrane.angle_type, 5.0 * membrane.eps, 180)
py_lmp.bond_coeff(membrane.bond_type, 'fene', 30.0 * membrane.eps / membrane.sigma**2,
                  1.5 * membrane.sigma, membrane.eps, membrane.sigma)

lj_factor = pow(2, 1./6)
# head-head & head-tail inter-lipid interaction
for bead_type in membrane.bead_types:
    py_lmp.pair_coeff(membrane.head_type, bead_type, membrane.eps,
                      0.95 * membrane.sigma * lj_factor,
                      0.95 * membrane.sigma * lj_factor, 'wca')
# tail-tail inter-lipid interaction
for i in range(len(membrane.tail_types)):
    for j in range(i, len(membrane.tail_types)):
        py_lmp.pair_coeff(membrane.tail_types[i], membrane.tail_types[j], membrane.eps,
                          1.0 * membrane.sigma * lj_factor,
                          (1.0 * lj_factor + membrane.wc) * membrane.sigma, 'wca')

# lipid-protein interaction
sol_lipid_eps = run_args.mem_int_eps
int_factors = (1.0, 0.5, 0.25)
# lipid-protein initial 0 interaction between all
rod_type_range = "{:d}*{:d}".format(simulation.get_min_rod_type(),
                                    simulation.get_max_rod_type())
mem_type_range = "{:d}*{:d}".format(min(membrane.bead_types),
                                    max(membrane.bead_types))
py_lmp.pair_coeff(mem_type_range, rod_type_range, 0.0, model.global_cutoff - 0.01, model.global_cutoff)
# lipid-protein volume-exclusion (all states)
vx_types = [t for t in model.all_bead_types if t not in model.active_bead_types]
for vx_type in vx_types:
    bead_lipid_contact = 0.5 * membrane.sigma * pow(2, 1./6) + model.bead_radii[vx_type]
    for mem_bead_type in membrane.bead_types:
        py_lmp.pair_coeff(mem_bead_type, vx_type + simulation.type_offset,
                          sol_lipid_eps, bead_lipid_contact, bead_lipid_contact, 'wca')
# lipid-protein tip interaction (all states)
tip_types = [state_struct[0][-1] + simulation.type_offset
             for state_struct in model.state_structures]
tip_lipid_contact = 0.5 * membrane.sigma * pow(2, 1./6) + model.rod_radius
tip_lipid_cutoff = tip_lipid_contact + run_args.mem_int_range
for tip_type in tip_types:
    for mem_bead_type, k in zip(membrane.bead_types, int_factors):
        py_lmp.pair_coeff(mem_bead_type, tip_type,
                          k*sol_lipid_eps, tip_lipid_contact, tip_lipid_cutoff, 'wca')
    
# ===== RODS ============================================================================
# GROUPS & COMPUTES
if hasattr(run_args, 'label_micelles'): 
    micelle_group = 'sol_tips'
    sol_tip_bead_type = model.state_structures[0][0][-1]
    py_lmp.variable(micelle_group, 'atom', '"type == {:d}"'.format(sol_tip_bead_type))
    py_lmp.group(micelle_group, 'dynamic', simulation.rods_group, 'var', micelle_group,
                 'every', out_freq)
    micelle_compute = "micelle_ID"
    if hasattr(run_args, 'micelle_cutoff'):
        micelle_cutoff = run_args.micelle_cutoff
    else:
        SS_tip_int_key = model.eps[(sol_tip_bead_type, sol_tip_bead_type)][1]
        SS_tip_int_range = model.int_types[SS_tip_int_key][1]
        micelle_cutoff = 2*model.rod_radius + SS_tip_int_range
    py_lmp.compute(micelle_compute, micelle_group, 'cluster/atom', micelle_cutoff)

#TODO label_fibrils ??

# FIXES & DYNAMICS
thermo_fix = 'thermostat'
py_lmp.fix(thermo_fix, 'all', 'langevin', run_args.temp, run_args.temp, run_args.damp, seed)#, "zero yes")

simulation.set_rod_dynamics('nve', opt = ['mol', model.rod_states[0]])

py_lmp.region('gcmc_init', 'block',
              membrane.xmin + model.rod_length/2, membrane.xmax - model.rod_length/2,
              membrane.ymin + model.rod_length/2, membrane.ymax - model.rod_length/2,
              0.0 + model.rod_length/2, run_args.Lz - model.rod_length/2)
concentration = run_args.conc / model.rod_length**3
simulation.set_state_concentration(0, concentration, 10, 10,
                                   opt = ['region', 'gcmc_init',
                                          'tfac_insert', 1.65])

# TEST DUMP...
# py_lmp.thermo_style('custom', 'step atoms', 'pe temp')
# py_lmp.variable('thermo_var', 'equal', '"stagger({:d}, 1)"'.format(out_freq))
# py_lmp.thermo('v_thermo_var')
# py_lmp.dump('test_dump', 'all', 'custom', out_freq, dump_path+'_init',
#             'id x y z type mol c_'+cluster_compute)
# py_lmp.dump_modify('test_dump', 'sort id')

# GENERATING INITIAL CONFIGURATION
py_lmp.neigh_modify('every', 1, 'delay', 1)
py_lmp.timestep(run_args.dt)
py_lmp.run(1000)
simulation.unset_state_concentration(0)
py_lmp.unfix(zwalls_fix)
py_lmp.reset_timestep(0)

# ===== MEMBRANE ========================================================================

# create membrane (box update, create membrane & groups, ...)
membrane.create_membrane(py_lmp, seed, append = True)

py_lmp.fix(zwalls_fix, 'all', 'wall/lj126',
           'zlo EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6),
           'zhi EDGE', 1.0, model.rod_radius, model.rod_radius*pow(2,1./6))
    
# GROUPS & COMPUTES
adsorbed_group = 'mem_and_tips'
py_lmp.variable('mem_and_tips', 'atom', '"' + 
                ' || '.join(['(type == {:d})'.format(t)
                             for t in tip_types + membrane.bead_types])
                + '"')
py_lmp.group(adsorbed_group, 'dynamic', 'all', 'var', 'mem_and_tips',
             'every', out_freq)
adsorbed_compute = 'mem_cluster'
py_lmp.compute(adsorbed_compute, adsorbed_group, 'aggregate/atom', tip_lipid_cutoff)
top_msd_compute = 'mem_top_msd'
py_lmp.compute(top_msd_compute, membrane.top_layer_group, 'msd')
bottom_msd_compute = 'mem_bottom_msd'
py_lmp.compute(bottom_msd_compute, membrane.bottom_layer_group, 'msd')
full_msd_compute = 'mem_full_msd'
py_lmp.compute(full_msd_compute, membrane.main_group, 'msd')

# FIXES & DYNAMICS
mem_dyn_fix = 'mem_dynamics'
press_eq_t = 1000 * run_args.dt
py_lmp.fix(mem_dyn_fix, membrane.main_group, 'nph',
           'x 0.0 0.0', press_eq_t, 'y 0.0 0.0', press_eq_t, 'couple xy',
           'dilate', 'all')
py_lmp.compute_modify(mem_dyn_fix + '_temp', 'dynamic/dof yes')
py_lmp.neigh_modify('exclude', 'molecule/intra', membrane.main_group)
# -> saves neighbour calculating without consequences (because beads 1 & 3 will
#    never be in range anyway & this doesn't affect bonds and angles)
py_lmp.special_bonds('fene') # not really necessary because of "neigh_modify exclude"

# MEMBRANE EQUILIBRATION
py_lmp.run(2000) #2 x press_eq_t, for good measure
py_lmp.reset_timestep(0)

# ===== FINAL ===========================================================================
py_lmp.variable('temp_x', 'equal', '"xhi - {:f}"'.format(model.rod_length / 2))
py_lmp.variable('temp_y', 'equal', '"yhi - {:f}"'.format(model.rod_length / 2))
py_lmp.region('gcmc_box', 'block',
              '-${temp_x}', '${temp_x}',
              '-${temp_y}', '${temp_y}',
              0.0 + 2 * model.rod_length, run_args.Lz - model.rod_length / 2)
N_rod_approx = int((membrane.xmax - membrane.xmin) *
                   (membrane.ymax - membrane.ymin) * 
                   run_args.Lz * concentration + 1)
mc_tries = int(0.01 * N_rod_approx + 1)
simulation.set_state_concentration(0, concentration, run_args.mc_every, mc_tries,
                                   opt = ['region', 'gcmc_box',
                                          'tfac_insert', 1.65])
# OUTPUT
dump_elems = "id x y z type mol"
try:
    dump_elems += " c_" + micelle_compute
except:
    pass
try:
    dump_elems += " c_" + adsorbed_compute
except:
    pass
py_lmp.dump("dump_cmd", "all", "custom", out_freq, dump_path, dump_elems)
py_lmp.dump_modify("dump_cmd", "sort id")
py_lmp.variable("area", "equal", "lx*ly")
py_lmp.thermo_style("custom", "step atoms", "pe temp", "press lx ly v_area")
py_lmp.thermo(out_freq)
# MSD data (only in x & y directions)
py_lmp.variable("top_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(top_msd_compute))
py_lmp.fix("mem_top_msd", "all", "ave/time", 1, 1, out_freq, "v_top_msd",
           "file", os.path.join(output_folder, sim_ID + '_mem_top.msd'))
py_lmp.variable("bottom_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(bottom_msd_compute))
py_lmp.fix("mem_bottom_msd", "all", "ave/time", 1, 1, out_freq, "v_bottom_msd",
           "file", os.path.join(output_folder, sim_ID + '_mem_bottom.msd'))
py_lmp.variable("full_msd", "equal", '"c_{0}[1] + c_{0}[2]"'.format(full_msd_compute))
py_lmp.fix("mem_full_msd", "all", "ave/time", 1, 1, out_freq, "v_full_msd",
           "file", os.path.join(output_folder, sim_ID + '_mem_full.msd'))

# RUN...
py_lmp.run(args.simlen)

MPI.Finalize()
