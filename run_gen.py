# encoding: utf-8
'''
A tool for generating run configuration files...

Created on 29 Mar 2019

@author: Eugen Rožić
'''
import sys, os, re
#import numpy as np

template_path = sys.argv[1]

out_folder = template_path + '_out'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

def match_replace(matchobj):
    return str(globals()[matchobj.group(1)])

def fill_template(line):
    return re.sub(r'<([a-zA-Z_][a-zA-Z_0-9\-]*)>', match_replace, line)

num_cells = 8.0
cell_sizes = [round(x,1) for x in [8.1 * 1.15**i for i in range(10)]]
for cell_size in cell_sizes:
    run_filename = "{:.1f}_{:.1f}.run".format(cell_size, num_cells)
    run_filepath = os.path.join(out_folder, run_filename)
    with open(run_filepath, 'w') as run_file:
        with open(template_path, 'r') as template_file:
            for line in template_file:
                run_file.write(fill_template(line))

# mem_wc = 1.4
# for mem_eps in (1/0.7, 1/0.8, 1/0.9, 1.0, 1/1.1, 1/1.2):
#     for int_eps in np.arange(0.5, 4.01, 0.25):
#         run_filename = "{:.2f}-{:.2f}-{:.2f}.run".format(mem_wc, mem_eps, int_eps)
#         run_filepath = os.path.join(out_folder, run_filename)
#         with open(run_filepath, 'w') as run_file:
#             with open(template_path, 'r') as template_file:
#                 for line in template_file:
#                     run_file.write(fill_template(line))
# 
# mem_eps = 1.0
# for mem_wc in (1.3, 1.5, 1.6, 1.7):
#     for int_eps in np.arange(0.5, 4.01, 0.25):
#         run_filename = "{:.2f}-{:.2f}-{:.2f}.run".format(mem_wc, mem_eps, int_eps)
#         run_filepath = os.path.join(out_folder, run_filename)
#         with open(run_filepath, 'w') as run_file:
#             with open(template_path, 'r') as template_file:
#                 for line in template_file:
#                     run_file.write(fill_template(line))

