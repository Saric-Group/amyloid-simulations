# encoding: utf-8
'''
An application for generating and submitting simulation jobs with "qsub"...

Created on 16 Oct 2018

@author: Eugen Rožić
'''

import argparse
import os, re
import subprocess
import traceback

#----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='An application for generating and submitting \
simulation jobs with "qsub"',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('job_template', type=str,
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

parser.add_argument('--args', type=str, default='',
                    help='arguments to pass to the job script as is')

job_args = parser.parse_args()

#----------------------------------------------------------------------------------------

def match_replace(matchobj):
    return str(globals()[matchobj.group(1)])

def fill_template(line):
    return re.sub(r'<([a-zA-Z_][a-zA-Z_0-9\-]*)>', match_replace, line)

walltime = job_args.walltime
memory = job_args.memory
cell_size = job_args.cell_size
num_cells = job_args.num_cells
sim_length = job_args.sim_length
args = job_args.args

temp_script_path = 'temp_job_script'

for cfg_file in job_args.cfgs:
    out_folder = os.path.splitext(cfg_file)[0]
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    with open(temp_script_path, 'w') as job_script:
        with open(job_args.job_template, 'r') as job_template:
            for line in job_template:
                job_script.write(fill_template(line))
    
    for n in range(job_args.repeat):
        try:
            subprocess.call(['qsub', 'temp_job_script'])
        except:
            traceback.print_exc()
            os.remove(temp_script_path)
            quit()
    
os.remove(temp_script_path)

