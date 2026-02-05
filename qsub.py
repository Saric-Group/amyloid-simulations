# encoding: utf-8
'''
An application for generating and submitting simulation jobs with "qsub"...

Created on 16 Oct 2018

@author: Eugen Rožić
'''

import argparse
import sys, os, re, shutil
import subprocess
import traceback
import time

script_loc = os.path.abspath(os.path.dirname(sys.argv[0]))

#----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='An application for generating and submitting \
simulation jobs with "qsub"', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('job_template',
                    help='the job template to populate and run')
parser.add_argument('cfg_file',
                    help='path to the model configuration file')
parser.add_argument('run_file',
                    help='path to the run configuration file')
parser.add_argument('simlen', type=int,
                    help='the length of the simulation')

parser.add_argument('--home', type=str, default=script_loc,
                    help='the location of the home directory')
parser.add_argument('--model', type=str, default="5p_v3",
                    help='model version')
parser.add_argument('--out', type=str, default=None,
                    help='the location of the output folder')

parser.add_argument('-t', '--walltime', type=str, default='24:00:00',
                    help='walltime for the job (HH:MM:SS)')
parser.add_argument('-m', '--memory', type=str, default='1G',
                    help='memory for the job')
parser.add_argument('-p', '--num_proc', type=str, default='8',
                    help='number of processors/cores')
parser.add_argument('-e', '--e_mail', type=str, default='',
                    help='e-mail address to send notification about jobs to')
parser.add_argument('-n', '--repeat', type=int, default=1,
                    help='number of copies of each job to start')

parser.add_argument('--args', type=str, default='',
                    help='any additional args for the program (to be passed verbatim)')
parser.add_argument('--aargs', type=str, default='',
                    help='any additional args for analysis (to be passed verbatim)')

job_args = parser.parse_args()

#----------------------------------------------------------------------------------------

def match_replace(matchobj):
    return str(globals()[matchobj.group(1)])

def fill_template(line):
    return re.sub(r'<([a-zA-Z_][a-zA-Z_0-9\-]*)>', match_replace, line)

cfg_file = job_args.cfg_file
run_file = job_args.run_file
simlen = job_args.simlen
project_home = job_args.home
walltime = job_args.walltime
memory = job_args.memory
num_proc = job_args.num_proc
mail_address = job_args.e_mail
if job_args.e_mail != '':
    e_mail = ''
else:
    e_mail = '#NOMAIL ' # to comment-out the e-mail PBS instructions in the job script
args = job_args.args
aargs = job_args.aargs

temp_script_path = 'temp_job_script'

if job_args.out == None:
    job_name = os.path.basename(job_args.job_template)
    cfg_name = os.path.splitext(os.path.basename(cfg_file))[0]
    out_folder = os.path.join(project_home, "data", job_name, job_args.model, cfg_name)
else:
    out_folder = job_args.out
if not os.path.exists(out_folder):
    os.makedirs(out_folder) #has to exist because of PBS script header

# copy cfg file to output folder parent, for later processing
cfg_copy = os.path.join(out_folder, "..", os.path.basename(cfg_file))
if not os.path.exists(cfg_copy):
    shutil.copy(cfg_file, cfg_copy)
    
with open(temp_script_path, 'w') as job_script,\
     open(job_args.job_template, 'r') as job_template:
    for line in job_template:
        job_script.write(fill_template(line))

# TODO?: implement job copies with qsub -j N:M option (job array)  
for n in range(job_args.repeat):
    try:
        subprocess.call(['qsub', temp_script_path])
        time.sleep(1)
    except:
        traceback.print_exc()
        os.remove(temp_script_path)
        quit()
    
os.remove(temp_script_path)
