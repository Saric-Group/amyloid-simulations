# encoding: utf-8
'''
Created on 16 May 2018

@author: Eugen Rožić
'''

import re
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

import lammps_multistate_rods.model as rod_model
import argparse

def state_struct_to_name(state_struct):
    for i in range(rod_model.num_states):
        if state_struct == rod_model.state_structures[i].replace('|',''):
            return rod_model.rod_states[i]
    return None

def wrap_periodic((elem, L)):
    if elem > L/2.:
        elem -= L
    elif elem < -L/2.:
        elem += L
    return elem

def parse_dump_file(dump_file_path):
    
    with open(dump_file_path, 'r') as dump_file:
        lines = dump_file.readlines()
    
    box_size = None
    data = []
    i = 0
    while i < len(lines):
        timestep = int(lines[i+1])
        atoms_to_read = int(lines[i+3])
        pattern = re.compile(r'([-\+\d\.eE]+) ([-\+\d\.eE]+)')
        xmin, xmax = [float(match) for match in pattern.match(lines[i+5]).groups()]
        ymin, ymax = [float(match) for match in pattern.match(lines[i+6]).groups()]
        zmin, zmax = [float(match) for match in pattern.match(lines[i+7]).groups()]
        if box_size is None:
            box_size = [(xmax - xmin), (ymax - ymin), (zmax - zmin)]
            
        #TODO get info on which column is what, don't force exact format... just require mol, type, x, y & z ...
        i += 9
        active_particles = []
        mol_states = []
        pattern = re.compile(r'(\d+) (\d+) ([-\+\d\.eE]+) ([-\+\d\.eE]+) ([-\+\d\.eE]+)')
        mol_id, bead_type, x, y, z = pattern.match(lines[i]).groups()
        current_mol_id = int(mol_id)
        current_mol_data = []
        current_mol_structure = []
        for n in range(atoms_to_read):
            mol_id, bead_type, x, y, z = pattern.match(lines[i+n]).groups()
            mol_id = int(mol_id)
            bead_type = int(bead_type)
            if mol_id != current_mol_id:
                active_particles.append(current_mol_data)
                mol_states.append(''.join(str(bead_type) for bead_type in current_mol_structure))
                current_mol_id = mol_id
                current_mol_data = []
                current_mol_structure = []
            current_mol_structure.append(bead_type)
            if bead_type in active_types:
                current_mol_data.append(np.array([float(x), float(y), float(z)]))
        active_particles.append(current_mol_data)
        mol_states.append(''.join(str(bead_type) for bead_type in current_mol_structure))
        
        data.append( [timestep, (active_particles, map(state_struct_to_name, mol_states))] )
        
        i += atoms_to_read
    
    return box_size, data

def analyse_snapshot((active_particles, mol_states), box_size):
    
    num_rods = len(active_particles)
    adj_mat = np.zeros((num_rods, num_rods), dtype='int')
    for i in range(num_rods):
        if len(active_particles[i]) == 1:
            delta_i = r_body
        else:
            delta_i = r_int
        adj_mat[i][i] = 1
        for bead_i in active_particles[i]:
            for j in range(i, num_rods):
                if adj_mat[i][j] != 0:
                    continue
                if len(active_particles[j]) == 1:
                    delta_j = r_body 
                else:
                    delta_j = r_int
                delta = delta_i + delta_j + cutoff_dist
                for bead_j in active_particles[j]:
                    dist = np.linalg.norm(map(wrap_periodic, zip((bead_i-bead_j),box_size)))
                    if dist < delta:
                        adj_mat[i][j] = adj_mat[j][i] = 1
                        break
                    
    adj_mat = csr_matrix(adj_mat)
    num_clusters, clusters = connected_components(adj_mat, directed=False)
    cluster_labels, cluster_sizes = np.unique(clusters, return_counts=True)
    sizes, size_occurences = np.unique(cluster_sizes, return_counts=True)
    
    #TODO use mol_states for analysis of cluster types, e.g. fibril, micelle or mixed...
    
    ret = [0]*max(sizes)
    for size,n in zip(sizes,size_occurences):
        ret[size-1] = n
    return ret

def analyse_dump_file(dump_file, config_file):
    #get information from the config file
    rod_model.set_model_params(config_file)
    global r_body, r_int, cutoff_dist, active_types
    r_body = rod_model.rod_radius
    cutoff_dist = 0.5*rod_model.rod_radius #??
    r_int = rod_model.int_radius
    active_types = set()
    for key, value in rod_model.eps.iteritems():
        if value != rod_model.vx and value > 0:
            active_types.update(key)
    active_types = list(active_types)

    box_size, raw_data = parse_dump_file(dump_file)
    n_snapshots = len(raw_data)

    cluster_data = raw_data

    for i in range(n_snapshots):
        cluster_data[i][1] = analyse_snapshot(raw_data[i][1], box_size)
    
    return cluster_data #a list of [ timestep, [ cluster occurences ] ] {cluster size implicit as index+1}

def output_cluster_data(cluster_data, output_path):
    
    with open(output_path, 'w') as out_file:
        for timestamp, data in cluster_data:
            out_file.write('{:^10d} | {:s}\n'.format(timestamp, ','.join(str(elem) for elem in data)))

def read_cluster_data(in_path):
    
    cluster_data = []
    
    with open(in_path, 'r') as data_file:
        for line in data_file:
            timestep, data = line.split(' | ')
            cluster_data.append([ int(timestep), [int(elem) for elem in data.split(',')] ])
    
    return cluster_data

# if called as a program does the cluster analysis of a single dump file
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='''
        Application for the analysis of clusters of lammps_multistate_rods as molecules
        from LAMMPS dump files with the structure: mol type x y z
        ''')
    parser.add_argument('config_file', help='path to the "lammps_multistate_rods" model config file describing the molecules')
    parser.add_argument('in_file', help='path of the dump file to analyse')
    args = parser.parse_args()
    
    processed_data = analyse_dump_file(args.in_file, args.config_file)

    output_cluster_data(processed_data, os.path.splitext(args.in_file)[0]+"_cluster_data")
        