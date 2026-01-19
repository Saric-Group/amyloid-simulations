# encoding: utf-8
'''
Draws energies and forces for different LAMMPS potentials (nm/cut, morse, gauss/cut)

Created on 30 Apr 2018

@author: Eugen Rožić
'''

import os
import matplotlib.pyplot as plt
import math

from lammps import lammps
import sys

lmp = lammps()
#lmp.cmd.log("none")

lmp.cmd.units("lj")
lmp.cmd.dimension(3)
lmp.cmd.atom_style("atomic")

lmp.cmd.boundary("p p p")
box_size = 10
lmp.cmd.region("box", "block", -box_size / 2, box_size / 2, -box_size / 2, box_size / 2, -box_size / 2, box_size / 2)
lmp.cmd.create_box(3, "box")

lmp.cmd.mass("*", 1)

# general potential parameters
eps = 1.0
try:
    r_body = float(sys.argv[1])
except:
    r_body = 1.0
r_int = 1.0*r_body
sigma = r_body + r_body

# output parameters
output_filename = "interactions.dat"
if os.path.isfile(output_filename):
    os.remove(output_filename)

num_points = 1000
min_r = 0.1*sigma
max_r = 3*sigma

#########################################################################################
### interaction definitions #############################################################
#########################################################################################

for lj_cutoff in (2.6, 3.0, 3.75):
    lmp.cmd.pair_style("lj/cut", max_r)
    #lmp.cmd.pair_modify("shift yes")
    lmp.cmd.pair_coeff("*", "*", eps, sigma/math.pow(2, 1./6), lj_cutoff*r_body)
    lmp.cmd.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'LJ_12-6_'+str(lj_cutoff))

# for i in (9,):
#     lmp.cmd.pair_style("nm/cut", max_r)
#     #lmp.cmd.pair_modify("shift yes")
#     lmp.cmd.pair_coeff("*", "*", eps, sigma, 12, i, cutoff)
#     lmp.cmd.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'LJ_12-'+str(i))
      
# for std_dev in (0.5, 0.6):
#     std_dev *= r_body
#     lmp.cmd.pair_style("gauss/cut", max_r)
#     #lmp.cmd.pair_modify("shift yes")
#     lmp.cmd.pair_coeff("*", "*", -math.sqrt(2*math.pi)*std_dev*eps, sigma, std_dev)
#     lmp.cmd.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Gauss-{:0.2f}'.format(std_dev))
       
# for a in (2.5,):
#     a /= r_body
#     lmp.cmd.pair_style("morse", max_r)
#     #lmp.cmd.pair_modify("shift yes")
#     lmp.cmd.pair_coeff("*", "*", eps, a, sigma, cutoff)
#     lmp.cmd.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'Morse-'+str(a))

for reach in (0.0, 0.1, 1.0):
    lmp.cmd.pair_style("cosine/squared", max_r)
    lmp.cmd.pair_coeff("*", "*", eps, sigma, sigma + reach*r_body, "wca")
    lmp.cmd.pair_write(1, 1, num_points, 'r', min_r, max_r, output_filename, 'cos-sq_'+str(reach))
    

#########################################################################################
#########################################################################################
#########################################################################################

data = []
with open(output_filename, 'r') as data_file:
    line = data_file.readline() # pair_write comment line
    data_file.readline() # 1st interaction comment line
    while (line!=''):
        data_file.readline() # empty line
        name = data_file.readline().strip()
        rs = []
        Es = []
        Fs = []
        data.append([name, rs, Es, Fs])
        N = int(data_file.readline().split()[1])
        data_file.readline() # empty line
        for i in range(N):
            parts = data_file.readline().split()
            rs.append(float(parts[1]))
            Es.append(float(parts[2]))
            Fs.append(float(parts[3]))
        line = data_file.readline() #next interaction comment line, or empty
n = len(data)

axis_font = {'size':14}
cmap = plt.cm.get_cmap('nipy_spectral')
colors = [cmap((i+1.0)/(n+1)) for i in range(n)]
#colors = ['b','g','r','c','m','y']

E_fig, E_axes = plt.subplots(num='LAMMPS interactions - energies', figsize=(10,8))
F_fig, F_axes = plt.subplots(num='LAMMPS interactions - forces', figsize=(10,8))

for i in range(n):
    E_axes.plot(data[i][1], data[i][2], label=data[i][0], color=colors[i], lw=1.5)
    F_axes.plot(data[i][1], data[i][3], label=data[i][0], color=colors[i], lw=1.5)

E_axes.legend(loc='lower right', prop = axis_font)
E_axes.set_xlabel(r'$r$', **axis_font)
E_axes.set_ylabel(r'$E(r)$', rotation='vertical', **axis_font)
E_axes.axis(xmin = 0, xmax = 4.5, ymin = -15, ymax = 5)
E_axes.axvline(sigma+1.0*r_body, color='black', linestyle='--', linewidth=0.5)
E_axes.axvline(sigma+1.5*r_body, color='black', linestyle='--', linewidth=0.5)
E_axes.grid(True)

F_axes.legend(loc='lower right', prop = axis_font)
F_axes.set_xlabel(r'$r$', **axis_font)
F_axes.set_ylabel(r'$F(r)$', **axis_font)
force_min = min(data[0][3])
F_axes.axis(xmin = 0, xmax = 4.5, ymin = 1.05*force_min, ymax = -0.2*force_min)
F_axes.axvline(sigma+1.0*r_body, color='black', linestyle='--', linewidth=0.5)
F_axes.axvline(sigma+1.5*r_body, color='black', linestyle='--', linewidth=0.5)
F_axes.grid(True)

plt.show()