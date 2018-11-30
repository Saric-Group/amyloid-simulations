# encoding: utf-8
'''
TODO

Created on 16 May 2018

@author: Eugen Rožić
'''

import os
import argparse

from lammps_multistate_rods import Model
from lammps_multistate_rods.tools import cluster_analysis

parser = argparse.ArgumentParser(description='Application for the analysis of clusters '\
                                 'of lammps_multistate_rods from LAMMPS dump files',
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
    
model = Model(args.config_file)
    
for in_file in args.in_files:
    output_path = os.path.splitext(in_file)[0]+"_cluster_data"
    box_size, timesteps, raw_data = cluster_analysis.parse_dump_file(
        in_file, args.every, model, args.type_offset)
    cluster_analysis.output_raw_data(box_size, timesteps, raw_data, output_path)
