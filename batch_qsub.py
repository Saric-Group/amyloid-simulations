# encoding: utf-8
'''
An application for generating and submitting simulation jobs with "qsub"...

Created on 16 Oct 2018

@author: Eugen Rožić
'''

import argparse

parser = argparse.ArgumentParser(description='TODO',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('job', type=str,
                    help='the job template to populate and run')
parser.add_argument('cell_size', type=float,
                    help='size of an SC cell (i.e. room for one rod)')
parser.add_argument('num_cells', type=float,
                    help='the number of cells per dimension')
parser.add_argument('sim_length', type=int,
                    help='the total number of MD steps to simulate')
parser.add_argument('cfgs', nargs='+', type=str,
                    help='paths to .cfg files to run as simulations')

parser.add_argument('-t', '--walltime', type=str, default='24:00:00',
                    help='walltime for the job (HH:MM:SS)')
parser.add_argument('-m', '--memory', type=str, default='1G',
                    help='memory for the job')

parser.add_argument('-n', '--repeat', type=int, default=1,
                    help='number of copies of each job to start')

parser.add_argument('--args', type=str,
                    help='arguments to pass to the job script as is')

args = parser.parse_args()

#TODO...