# encoding: utf-8
'''
TODO

Created on 15 Aug 2018

@author: Eugen Rožić
'''
import os
import argparse

import numpy as np
from matplotlib import pyplot as plt

from lammps_multistate_rods.tools.clusters import read_cluster_data, free_rods

def concentration(number, volume):
    '''
    Concentration in M for volume in nm^3
    '''
    return (10.0*number)/(6.022*volume)

parser = argparse.ArgumentParser(description='''
        Application for plotting the free monomer concentrations''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('in_files', nargs='+', help='paths of the *_cluster_data files')
parser.add_argument('-n', '--avg_over', type=int, default=10, help='last how many timesteps to average over')
parser.add_argument('-s', '--save', type=str, default=None, help='saves the figure to the given location')
args = parser.parse_args()

data_dir = os.path.dirname(args.in_files[0])
data_len = len(args.in_files)
avg_over = args.avg_over

labels = [os.path.basename(os.path.dirname(args.in_files[0]))]
data = [[]]

for n in range(data_len):
    label = os.path.basename(os.path.dirname(args.in_files[n]))
    if (label != labels[-1]):
        labels.append(label)
        data.append([])
    
    box_size, timesteps, raw_data = read_cluster_data(args.in_files[n])
    volume = reduce(lambda x,y: x*y, box_size)
    
    analyse_from = len(raw_data)-avg_over
    results = free_rods(raw_data[analyse_from:])
    
    temp_x = temp_y = 0.0
    for free_monomers, total_monomers in results:
        if temp_x == 0.0:
            temp_x = total_monomers 
        elif total_monomers != temp_x:
            raise Exception('Total number of monomers is not constant!') #??
        temp_y += free_monomers
    temp_y /= avg_over

    data[-1].append((concentration(temp_x, volume),concentration(temp_y, volume)))

n = len(labels)

fig = plt.figure()

cmap = plt.cm.get_cmap('nipy_spectral')
colors = [cmap((i+1.0)/(n+1)) for i in range(n)]
#colors = ['b','g','r','c','m','y']

xmax = -np.infty
xmin = np.infty
for i in range(n):
    
    data[i].sort(key=lambda i: i[0])
    
    xs, ys = zip(*data[i])
    xs = np.log(np.array(xs))
    ys = np.log(np.array(ys))

    temp_max = xs.max()
    if temp_max > xmax:
        xmax = temp_max
    temp_min = xs.min()
    if temp_min < xmin:
        xmin = temp_min

    plt.plot(xs, ys, label=labels[i], c=colors[i], marker='o', lw=1.0)

plt.plot((xmin, xmax), (xmin, xmax), 'k--', lw=1.0)

plt.legend(loc='upper left')
plt.xlabel(r'$\ln\left(c_{total}\;/\;M\right)$')
plt.ylabel(r'$\ln\left(c_{free}\;/\;M\right)$', rotation='vertical')
plt.axis(ymin=xmin-0.5, ymax=xmax+0.5, xmin=xmin-0.5, xmax=xmax+0.5)
plt.grid(True)

if args.save:
    fig.savefig(args.save, dpi=1000, facecolor='white')

plt.show()

