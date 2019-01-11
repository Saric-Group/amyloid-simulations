# encoding: utf-8
'''
TODO

Created on 16 May 2018

@author: Eugen Rožić
'''

import os
import argparse

import lammps_multistate_rods as rods
import lammps_multistate_rods.tools as rods_tools

parser = argparse.ArgumentParser(description='Application for the processing of LAMMPS'\
                                 'dump files generated with lammps_multistate_rods library',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('config_file',
                    help='path to the "lammps_multistate_rods" model config file')
parser.add_argument('in_files', nargs='+',
                    help='path(s) of the dump file(s) to analyse')
parser.add_argument('-t', '--type_offset', type=int, default=0,
                    help='the type offset for the rod models of the simulation')
parser.add_argument('-n', '--every', type=int, default=1,
                    help='every which snapshot to analyse')

args = parser.parse_args()
    
model = rods.Model(args.config_file)
    
for in_file in args.in_files:
    raw_data = rods_tools.parse_dump_file(in_file)
    
    cluster_output_path = os.path.splitext(in_file)[0]+"_cluster_data"
    box_size, timesteps, cluster_data = rods_tools.clusters.get_cluster_data(
        raw_data, args.every, model, args.type_offset)
    
    rods_tools.clusters.write_cluster_data(box_size, timesteps, cluster_data,
                                           cluster_output_path)
    
    last_dump_output_path = os.path.splitext(in_file)[0]+"_last_dump"
    for timestep, box_bounds, data_structure, data in rods_tools.parse_dump_file(in_file):
        pass
    rods_tools.write_dump_snapshot(timestep, box_bounds, data_structure, data,
                                       last_dump_output_path)
