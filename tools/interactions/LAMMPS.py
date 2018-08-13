# encoding: utf-8
'''
Draws energies and forces for different LAMMPS potentials (nm/cut, morse, gauss/cut)

Created on 30 Apr 2018

@author: Eugen Rožić
'''

import os
import matplotlib.pyplot as plt
import math
import numpy as np

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
eps = 1.0
try:
    r_body = float(sys.argv[1])
except:
    r_body = 1.0
r_int = 0.5*r_body
sigma = r_int + r_int
cutoff = sigma + 1.5*r_body

# output parameters
output_filename = "interactions.dat"
if os.path.isfile(output_filename):
    os.remove(output_filename)

num_points = 1000
min_r = 0.5*sigma
max_r = cutoff

#########################################################################################
### interaction definitions #############################################################
#########################################################################################

for i in (6, 4):
    py_lmp.pair_style("nm/cut", cutoff)
    py_lmp.pair_modify("shift yes")
    py_lmp.pair_coeff("*", "*", eps, sigma, 12, i)
    py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'LJ_12-'+str(i))
    
# # Buckingham
# for alpha in (8, ):
#     py_lmp.pair_style("buck", cutoff)
#     py_lmp.pair_modify("shift yes")
#     py_lmp.pair_coeff("*", "*", 6*eps*math.exp(alpha)/(alpha-6), sigma/alpha, eps*alpha*sigma**6/(alpha-6))
#     py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Buckingham-'+str(alpha)) 

# Gaussian
for std_dev in (0.5*r_body, 0.6*r_body):
    py_lmp.pair_style("gauss/cut", cutoff)
    py_lmp.pair_modify("shift yes")
    py_lmp.pair_coeff("*", "*", -math.sqrt(2*math.pi)*std_dev*eps, sigma, std_dev)
    py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Gauss-{:0.2f}'.format(std_dev))
      
# Morse
for a in (2.5/r_body, 2.0/r_body):
    py_lmp.pair_style("morse", cutoff)
    py_lmp.pair_modify("shift yes")
    py_lmp.pair_coeff("*", "*", eps, a, sigma)
    py_lmp.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Morse-'+str(a))                          
    

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

axis_font = {'size':12}
fig = plt.figure('LAMMPS interactions - energies', figsize=(10,8))

#LJ
plt.plot(data[0][1], data[0][2], label=data[0][0], color='black', linewidth=2)
plt.plot(data[1][1], data[1][2], label=data[1][0], color='black', linewidth=1)
#gauss
plt.plot(data[2][1], data[2][2], label=data[2][0], color='green', linewidth=2)
plt.plot(data[3][1], data[3][2], label=data[3][0], color='green', linewidth=1)
#morse
plt.plot(data[4][1], data[4][2], label=data[4][0], color='red', linewidth=2)
plt.plot(data[5][1], data[5][2], label=data[5][0], color='red', linewidth=1)

# from potentials import *
# import numpy as np
# rs = np.linspace(min_r, max_r, num_points)
# plt.plot(rs, [lj_n_m(12, 4, r, sigma, eps) for r in rs], 'k-', label='LJ_12-4')
# plt.plot(rs, [gauss(0.5, r, sigma, eps) for r in rs], 'g-', label='Gauss_0.5')
# plt.plot(rs, [morse(3, r, sigma, eps) for r in rs], 'r-', label='Morse_3')

plt.legend(loc='lower right')
plt.xlabel(r'$r$', **axis_font)
plt.ylabel(r'$E(r)$', rotation='vertical', **axis_font)
plt.axis(xmin = min_r, xmax = max_r, ymin = -1.05*eps, ymax = 1.05*eps)
plt.axvline(sigma+r_body, color='black', linestyle='--', linewidth=0.5)


fig = plt.figure('LAMMPS interactions - forces', figsize=(10,8))

#LJ
plt.plot(data[0][1], data[0][3], label=data[0][0], color='black', linewidth=2)
plt.plot(data[1][1], data[1][3], label=data[1][0], color='black', linewidth=1)
#gauss
plt.plot(data[2][1], data[2][3], label=data[2][0], color='green', linewidth=2)
plt.plot(data[3][1], data[3][3], label=data[3][0], color='green', linewidth=1)
#morse
plt.plot(data[4][1], data[4][3], label=data[4][0], color='red', linewidth=2)
plt.plot(data[5][1], data[5][3], label=data[5][0], color='red', linewidth=1)


plt.legend(loc='lower right')
plt.xlabel(r'$r$', **axis_font)
plt.ylabel(r'$F(r)$', **axis_font)
force_min = min(data[0][3])
plt.axis(xmin = min_r, xmax = max_r, ymin = 1.05*force_min, ymax = -1.05*force_min)
plt.axvline(sigma+r_body, color='black', linestyle='--', linewidth=0.5)

#fig.savefig("./LAMMPS_potenatials.pdf", dpi=1000)
plt.show()

