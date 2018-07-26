# encoding: utf-8
'''
Created on 16 May 2018

@author: Eugen Rožić
'''

import re
import os
import argparse

import numpy as np

import lammps_multistate_rods

def state_struct_to_id(state_struct):
    '''
    Returns rod state id from its structure (list of atom types)
    '''
    for i in range(model.num_states):
        if state_struct == model.state_structures[i].replace('|',''):
            return i
    return None

def wrap_periodic((elem, L)):
    '''
    Wraps position (distance) along a periodic boundary.
    '''
    if elem > L/2.:
        elem -= L
    elif elem < -L/2.:
        elem += L
    return elem

def parse_dump_file(dump_file_path):
    '''
    Reads the dump file given by the path and extracts data from it.
    
    returns : a list of [timestep, snapshot_data] pairs, where "snapshot_data" is a dictionary by
    cluster ID's whose values are lists of (molecule ID, rod state ID) pairs
    '''
    with open(dump_file_path, 'r') as dump_file:
        
        data = []
        #box_size = []
        #TODO get info on which column is what, don't force exact format...
        pattern = re.compile(r'(\d+) ([-\+\d\.eE]+) ([-\+\d\.eE]+) ([-\+\d\.eE]+) (\d+) (\d+) (\d+)')
        snapshot_data = {}
        current_mol_id = None
        current_cluster_id = None
        current_rod_structure = []
        i = 0
        max_row = 0
        for line in dump_file:
            i += 1
            if i == 2:
                timestep = int(line)
            
            elif i == 4:
                atoms_to_read = int(line)
                max_row = 9 + atoms_to_read
            
#             elif i in (6,7,8) and len(box_size) < 3:
#                 l_bound, r_bound = [float(match) for match in re.compile(r'([-\+\d\.eE]+) ([-\+\d\.eE]+)').match(line).groups()]
#                 box_size.append(r_bound - l_bound)
            
            elif i >= 10 and i <= max_row:
                bead_id, x, y, z, bead_type, mol_id, cluster_id = pattern.match(line).groups()
                mol_id = int(mol_id)
                cluster_id = int(cluster_id)
                if current_mol_id is None:
                    current_mol_id = mol_id
                elif mol_id != current_mol_id:
                    current_rod_state = state_struct_to_id(''.join(bead_type for bead_type in current_rod_structure))
                    try:
                        snapshot_data[current_cluster_id].append((current_mol_id, current_rod_state))
                    except KeyError:
                        snapshot_data[current_cluster_id] = [(current_mol_id, current_rod_state)]
                    current_mol_id = mol_id
                    current_rod_structure = []
                current_rod_structure.append(bead_type)
                if cluster_id > 0:
                    current_cluster_id = cluster_id
            
            if i == max_row:
                current_rod_state = state_struct_to_id(''.join(bead_type for bead_type in current_rod_structure))
                try:
                    snapshot_data[current_cluster_id].append((current_mol_id, current_rod_state))
                except KeyError:
                    snapshot_data[current_cluster_id] = [(current_mol_id, current_rod_state)]
                data.append([timestep, snapshot_data])
                snapshot_data = {}
                current_mol_id = None
                current_cluster_id = None
                current_rod_structure = []
                i = 0
    return data

#TODO provide more advanced analysis (like changes in clusters between snapshots)
def analyse_data(raw_data):
    '''
    Analyses the output of "parse_dump_file"...
    
    return : a list of [timestep, cluster_sizes] pairs, where "cluster_sizes" is a dictionary by
    cluster type (same as rod state ID if homogeneous, otherwise -1) whose values are
    lists of (cluster size, occurrences) pairs 
    '''
    cluster_data = []
    
    for timestep, snapshot_data in raw_data:
        cluster_sizes = {}
        for cluster in snapshot_data.values():
            cluster_type = cluster[0][1] # cluster_type is the same as state id of rods if they are all in the same state
            for rod_id, rod_state in cluster:
                if rod_state != cluster_type:
                    cluster_type = -1
            try:
                cluster_sizes[cluster_type].append(len(cluster))
            except KeyError:
                cluster_sizes[cluster_type] = [len(cluster)]
    
        for cluster_type, value in cluster_sizes.iteritems():
            cluster_sizes[cluster_type] = zip(*np.unique(value, return_counts=True))
        
        cluster_data.append([timestep, cluster_sizes])
        
    return cluster_data

def analyse_dump_file(dump_file, config_file):
    '''
    return : the output of the "analyse_data" method
    '''
    global model
    model = lammps_multistate_rods.Model(config_file)

    raw_data = parse_dump_file(dump_file)
    
    cluster_data = analyse_data(raw_data)
    
    return cluster_data

def output_cluster_data(cluster_data, output_path):
    
    with open(output_path, 'w') as out_file:
        for timestamp, cluster_sizes in cluster_data:
            out_str = []
            for key, value in cluster_sizes.iteritems():
                out_str.append('({})'.format(key))
                i = 1
                for size, occurrence in value:
                    while i < size:
                        out_str.append('0')
                        i += 1
                    out_str.append('{:d}'.format(occurrence))
                    i += 1
            out_file.write('{:^10d} | {:s}\n'.format(timestamp, ' '.join(out_str)))

def read_cluster_data(in_path):
    
    cluster_data = []
    
    with open(in_path, 'r') as data_file:
        for line in data_file:
            cluster_sizes = {}
            timestep, raw_data = line.split(' | ')
            chunks = re.compile(r'\(([-\d]+)\) ([ \d]+)').findall(raw_data)
            for chunk in chunks:
                key = int(chunk[0])
                value = []
                i = 1
                for elem in chunk[1].split(' '):
                    if elem != '':
                        elem = int(elem)
                        if elem > 0:
                            value.append((i,elem))
                        i += 1
                cluster_sizes[key] = value
            cluster_data.append([int(timestep), cluster_sizes])
            
    return cluster_data

# if called as a program does the cluster analysis of a single dump file
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='''
        Application for the analysis of clusters of lammps_multistate_rods as molecules
        from LAMMPS dump files with the output structure: id x y z type mol c_''' +
        lammps_multistate_rods.Simulation.cluster_compute)
    parser.add_argument('config_file', help='path to the "lammps_multistate_rods" model config file')
    parser.add_argument('in_file', help='path of the dump file to analyse')
    args = parser.parse_args()
    
    output_path = os.path.splitext(args.in_file)[0]+"_cluster_data"
    
    processed_data = analyse_dump_file(args.in_file, args.config_file)

    output_cluster_data(processed_data, output_path)
        