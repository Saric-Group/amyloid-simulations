# encoding: utf-8
'''
Created on 16 May 2018

@author: Eugen Rožić

TODO
'''

import re
import os
import argparse

import numpy as np

import lammps_multistate_rods

def wrap_periodic((elem, L)):
    '''
    Wraps position (distance) along a periodic boundary.
    '''
    if elem > L/2.:
        elem -= L
    elif elem < -L/2.:
        elem += L
    return elem

def parse_dump_file(dump_file_path, every, model, particle_offset=0, type_offset=0):
    '''
    Reads the dump file given by the path and extracts data from it.
    
    model : the lammps_multistate_rods.Model class instance which was used to generate the data
    
    returns : a list of box dimensions, a list of timesteps and a corresponding list of snapshot_data,
    where "snapshot_data" is a dictionary by cluster ID's whose values are lists of
    (rod/mol ID, rod state ID) pairs
    
    NOTE: cluster ID's are set to the lowest rod/mol ID of each cluster, not the ones that are in
    the dump file
    '''
    states_types = [ [ elem + type_offset 
                        for patch in state_struct for elem in patch]
                            for state_struct in model.state_structures]
    
    def state_types_to_id(state_types):
        '''
        Returns rod state id from its structure (list of atom types)
        '''
        for i in range(model.num_states):
            if state_types == states_types[i]:
                return i
        return None
    
    with open(dump_file_path, 'r') as dump_file:
        
        box_size = []
        timesteps = []
        raw_data = []
        #TODO get info on which column is what, don't force exact format...
        pattern = re.compile(r'(\d+) ([-\+\d\.eE]+) ([-\+\d\.eE]+) ([-\+\d\.eE]+) (\d+) (\d+) (\d+)')
        snapshot_data = {}
        current_mol_id = None
        current_cluster_id = None
        current_rod = []
        count = 0
        skip = False
        min_row = max_row = 0
        i = 0
        for line in dump_file:
            i += 1
            
            if skip and i != 1:
                if i == max_row:
                    i = 0
                continue
            
            if i == 1:
                count += 1
                if (count-1) % every != 0:
                    skip = True
                    continue
                else:
                    skip = False
                    
            elif i == 2:
                timesteps.append(int(line))
            
            elif i == 4:
                atoms_to_read = int(line)
                max_row = 9 + atoms_to_read
                min_row = 10 + particle_offset
            
            elif i in (6,7,8) and len(box_size) < 3:
                l_bound, r_bound = [float(match) for match in re.compile(r'([-\+\d\.eE]+) ([-\+\d\.eE]+)').match(line).groups()]
                box_size.append(r_bound - l_bound)
               
            elif i >= min_row and i <= max_row:
                bead_id, x, y, z, bead_type, mol_id, cluster_id = pattern.match(line).groups()
                mol_id = int(mol_id)
                cluster_id = int(cluster_id)
                if current_mol_id is None:
                    current_mol_id = mol_id
                elif mol_id != current_mol_id:
                    current_rod_state = state_types_to_id(current_rod)
                    try:
                        snapshot_data[current_cluster_id].append((current_mol_id, current_rod_state))
                    except KeyError:
                        snapshot_data[current_cluster_id] = [(current_mol_id, current_rod_state)]
                    current_mol_id = mol_id
                    current_cluster_id = 0
                    current_rod = []
                current_rod.append(int(bead_type))
                if cluster_id > current_cluster_id:
                    current_cluster_id = cluster_id
            
            if i == max_row:
                current_rod_state = state_types_to_id(current_rod, model)
                try:
                    snapshot_data[current_cluster_id].append((current_mol_id, current_rod_state))
                except KeyError:
                    snapshot_data[current_cluster_id] = [(current_mol_id, current_rod_state)]
                
                #switch keys to correspond to lowest mol_id in each cluster
                new_snapshot_data = {}
                for value in snapshot_data.values():
                    new_snapshot_data[value[0][0]] = value
                
                raw_data.append(new_snapshot_data)
                snapshot_data = {}
                current_mol_id = None
                current_cluster_id = None
                current_rod = []
                i = 0
    
    return box_size, timesteps, raw_data

def output_raw_data(box_size, timesteps, raw_data, output_path):
    '''
    Outputs the data returned by "parse_dump_file"
    '''
    with open(output_path, 'w') as out_file:
        out_file.write('{:f} {:f} {:f}\n'.format(*box_size))
        for timestep, snapshot_data in zip(timesteps, raw_data):
            out_file.write('{:^10d} | {:s}\n'.format(timestep, str(snapshot_data)))

# if called as a program just parses and outputs raw rod cluster data from the given dump files
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=
        '''Application for the analysis of clusters of lammps_multistate_rods as molecules
        from LAMMPS dump files with the output structure: id x y z type mol c_''' +
        lammps_multistate_rods.Simulation.cluster_compute,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config_file',
                        help='path to the "lammps_multistate_rods" model config file')
    parser.add_argument('in_files', nargs='+',
                        help='path(s) of the dump file(s) to analyse')
    parser.add_argument('-n', '--every', type=int, default=1,
                        help='every which snapshot to analyse')
    args = parser.parse_args()
    
    model = lammps_multistate_rods.Model(args.config_file)
    
    for in_file in args.in_files:
        output_path = os.path.splitext(in_file)[0]+"_cluster_data"
        box_size, timesteps, raw_data = parse_dump_file(in_file, args.every, model)
        output_raw_data(box_size, timesteps, raw_data, output_path)

def read_raw_data(input_path):
    '''
    Reads what was output with "output_raw_data"
    
    returns : a box size triplet, a list of timesteps and a corresponding list of snapshot_data,
    where "snapshot_data" is a dictionary by cluster ID's whose values are lists of 
    (rod/mol ID, rod state ID) pairs.
    '''
    timesteps = []
    raw_data = []
    
    with open(input_path, 'r') as in_file:
        box_size = map(float, in_file.readline().split())
        for line in in_file:
            timestep, snapshot_data = line.split(' | ')
            timesteps.append(int(timestep))
            raw_data.append(eval(snapshot_data))
            
    return box_size, timesteps, raw_data

def cluster_sizes_by_type(raw_data):
    '''
    raw_data : a list of "snapshot_data", dictionaries by cluster ID whose values are lists of
    (rod/mol ID, rod state ID) pairs
    
    return : a list of "cluster_sizes", dictionaries by cluster type (same as rod state ID if
    homogeneous, otherwise -1) whose values are pairs of (cluster_sizes, occurrences) lists, and
    a number equal to the maximum cluster size across all timesteps and types.
    '''
    ret = [None]*len(raw_data)
    max_size = 0
    i = 0
    for snapshot_data in raw_data:
        cluster_sizes = {}
        for cluster in snapshot_data.values():
            cluster_type = cluster[0][1] # cluster_type is the same as state id of rods if they are all in the same state
            for rod_id, rod_state in cluster:
                if rod_state != cluster_type:
                    cluster_type = -1
            
            cluster_size = len(cluster)
            try:
                cluster_sizes[cluster_type].append(cluster_size)
            except KeyError:
                cluster_sizes[cluster_type] = [cluster_size]
            
            if cluster_size > max_size:
                max_size = cluster_size
    
        for cluster_type, value in cluster_sizes.iteritems():
            cluster_sizes[cluster_type] = np.unique(value, return_counts=True)
        
        ret[i] = cluster_sizes
        i += 1
        
    return ret, max_size

def free_monomers(raw_data, total=True):
    '''
    raw_data : a list of "snapshot_data", dictionaries by cluster ID whose values are lists of
    (rod/mol ID, rod state ID) pairs
    
    return : a list of numbers of free monomers (of any type) in each of snapshot_data. If
    total=True a list of (free monomers, total monomers) pairs is returned instead
    '''
    ret = [None]*len(raw_data)
    i = 0
    for snapshot_data in raw_data:
        free_monomers = 0
        total_monomers = 0
        for cluster in snapshot_data.values():
            cluster_size = len(cluster)
            total_monomers += cluster_size
            if cluster_size == 1:
                free_monomers += 1
        if total:
            ret[i] = (free_monomers, total_monomers)
        else:
            ret[i] = free_monomers
        i += 1
    
    return ret

#TODO provide more advanced analysis (like changes in clusters between snapshots)

