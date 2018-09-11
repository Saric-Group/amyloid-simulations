# encoding: utf-8
'''
Draws energies and forces for different LAMMPS potentials (nm/cut, morse, gauss/cut)

Created on 30 Apr 2018

@author: Eugen Rožić
'''

import os
import matplotlib.pyplot as plt
import math

from lammps import PyLammps
import sys

py_lmp = PyLammps()
#py_lmp.log("none")

py_lmp.units("lj")
py_lmp.dimension(3)
py_lmp.atom_style("atomic")

py_lmp.boundary("p p p")
box_size = 10
py_lmp.region("box", "block", -box_size / 2, box_size / 2, -box_size / 2, box_size / 2, -box_size / 2, box_size / 2)
py_lmp.create_box(3, "box")

py_lmp.mass("*", 1)

# general potential parameters
eps = 4.5
try:
    r_body = float(sys.argv[1])
except:
    r_body = 1.0
r_int = 0.25*r_body
sigma = r_body + r_body
cutoff = sigma + 1.0*r_body

# output parameters
output_filename = "interactions.dat"
if os.path.isfile(output_filename):
    os.remove(output_filename)

num_points = 1000
min_r = 0.1*sigma
max_r = 5

#########################################################################################
### interaction definitions #############################################################
#########################################################################################

py_lmp.pair_style("lj/cut", max_r)
#py_lmp.pair_modify("shift yes")
py_lmp.pair_coeff("*", "*", eps, sigma/math.pow(2, 1./6), cutoff)
py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'LJ_12-6')

# for i in (9,):
#     py_lmp.pair_style("nm/cut", max_r)
#     #py_lmp.pair_modify("shift yes")
#     py_lmp.pair_coeff("*", "*", eps, sigma, 12, i, cutoff)
#     py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'LJ_12-'+str(i))
      
# for std_dev in (0.5, 0.6):
#     std_dev *= r_body
#     py_lmp.pair_style("gauss/cut", max_r)
#     #py_lmp.pair_modify("shift yes")
#     py_lmp.pair_coeff("*", "*", -math.sqrt(2*math.pi)*std_dev*eps, sigma, std_dev)
#     py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Gauss-{:0.2f}'.format(std_dev))
       
for a in (2.5,):
    a /= r_body
    py_lmp.pair_style("morse", max_r)
    #py_lmp.pair_modify("shift yes")
    py_lmp.pair_coeff("*", "*", eps, a, sigma, cutoff)
    py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Morse-'+str(a))

for reach in (1.0,):
    py_lmp.pair_style("cosine/squared", max_r, "wca yes")
    py_lmp.pair_coeff("*", "*", eps, sigma, sigma + reach*r_body)
    py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'cos_sq_'+str(reach))
    

#########################################################################################
#########################################################################################
#########################################################################################

data = []
with open(output_filename, 'r') as data_file:
    line = data_file.readline()
    while (line!=''):
        data_file.readline()
        name = data_file.readline().strip()
        rs = []
        Es = []
        Fs = []
        data.append([name, rs, Es, Fs])
        N = int(data_file.readline().split()[1])
        data_file.readline()
        for i in range(N):
            parts = data_file.readline().split()
            rs.append(float(parts[1]))
            Es.append(float(parts[2]))
            Fs.append(float(parts[3]))
        line = data_file.readline()    
n = len(data)

axis_font = {'size':12}
cmap = plt.cm.get_cmap('nipy_spectral')
colors = [cmap((i+1.0)/(n+1)) for i in range(n)]
#colors = ['b','g','r','c','m','y']

E_fig, E_axes = plt.subplots(num='LAMMPS interactions - energies', figsize=(10,8))
F_fig, F_axes = plt.subplots(num='LAMMPS interactions - forces', figsize=(10,8))

for i in range(n):
    E_axes.plot(data[i][1], data[i][2], label=data[i][0], color=colors[i], lw=1.5)
    F_axes.plot(data[i][1], data[i][3], label=data[i][0], color=colors[i], lw=1.5)

E_axes.legend(loc='lower right')
E_axes.set_xlabel(r'$r$', **axis_font)
E_axes.set_ylabel(r'$E(r)$', rotation='vertical', **axis_font)
E_axes.axis(xmin = 0, xmax = 4.5, ymin = -15, ymax = 5)
E_axes.axvline(sigma+1.0*r_body, color='black', linestyle='--', linewidth=0.5)
E_axes.axvline(sigma+1.5*r_body, color='black', linestyle='--', linewidth=0.5)
E_axes.grid(True)

F_axes.legend(loc='lower right')
F_axes.set_xlabel(r'$r$', **axis_font)
F_axes.set_ylabel(r'$F(r)$', **axis_font)
force_min = min(data[0][3])
F_axes.axis(xmin = 0, xmax = 4.5, ymin = 1.05*force_min, ymax = -0.2*force_min)
F_axes.axvline(sigma+1.0*r_body, color='black', linestyle='--', linewidth=0.5)
F_axes.axvline(sigma+1.5*r_body, color='black', linestyle='--', linewidth=0.5)
F_axes.grid(True)

plt.show()

