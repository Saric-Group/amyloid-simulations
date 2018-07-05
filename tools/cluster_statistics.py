# encoding: utf-8
'''
Created on 21 Jun 2018

@author: Eugen Rožić
'''
import os
import numpy as np
from math import sqrt
from cluster_analysis import read_cluster_data

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

import argparse

parser = argparse.ArgumentParser(description='''
        Application for calculating the statistics of cluster data over multiple simulations
        ''')
parser.add_argument('in_files', nargs='+', help='paths of the *_cluster_data files')
args = parser.parse_args()

data_dir = os.path.dirname(args.in_files[0])
all_data = []
for in_file in args.in_files:
    if data_dir != os.path.dirname(in_file):
        raise Exception('All files should be from the same directory!')
    all_data.append(read_cluster_data(in_file))

n_snapshots = len(all_data[0])
n_simulations = len(all_data)

timesteps = [0]*n_snapshots
avgs = [0]*n_snapshots
avg_devs = [0]*n_snapshots
maxs = [0]*n_snapshots
max_devs = [0]*n_snapshots
mcs = [0]*n_snapshots
mc_devs = [0]*n_snapshots
fms = [0]*n_snapshots
fm_devs = [0]*n_snapshots
distributions = [None]*n_snapshots

for i in range(n_snapshots):
    timesteps[i] = all_data[0][i][0]
    
    avg_sizes = np.zeros(n_simulations)
    sizes_devs = np.zeros(n_simulations)
    max_sizes = np.zeros(n_simulations, int)
    mc_sizes = np.zeros(n_simulations, int)
    free_mons = np.zeros(n_simulations, int)
    for n in range(n_simulations):
        free_mons[n] = all_data[n][i][1][0]
        n_size_occurrences = np.array(all_data[n][i][1][1:])
        max_sizes[n] = len(n_size_occurrences) + 1
        if max_sizes[n] > 1: # else it should be 0, which it is by initialisation
            n_sizes = np.arange(2, max_sizes[n]+1)
            avg_sizes[n] = np.average(n_sizes, weights=n_size_occurrences)
            sizes_devs[n] = sqrt(np.average((n_sizes-avg_sizes[n])**2, weights=n_size_occurrences))
            mc_sizes[n] = np.argmax(n_size_occurrences) + 2
    temp = np.sum(1. / sizes_devs**2)
    avgs[i] = np.sum(avg_sizes / sizes_devs**2) / temp
    avg_devs[i] = 1. / sqrt(temp)
    maxs[i] = np.average(max_sizes)
    max_devs[i] = sqrt(np.average((max_sizes-maxs[i])**2))
    mcs[i] = np.average(mc_sizes)
    mc_devs[i] = sqrt(np.average((mc_sizes-mcs[i])**2))
    fms[i] = np.average(free_mons)
    fm_devs[i] = sqrt(np.average((free_mons-fms[i])**2))
            
    biggest_cluster = max(max_sizes)
    aggregate_occurrences = [None]*biggest_cluster
    for j in range(biggest_cluster):
        aggregate_occurrences[j] = np.zeros(n_simulations)
        for n in range(n_simulations):
            try:
                aggregate_occurrences[j][n] = all_data[n][i][1][j]
            except IndexError:
                pass # if not there should be 0, which it is by initialisation
    avg_occurrences = [np.average(elem) for elem in aggregate_occurrences]
    dev_occurrences = [sqrt(np.average((aggregate_occurrences[j] - avg_occurrences[j])**2)) for j in range(biggest_cluster)]
    distributions[i] = (avg_occurrences, dev_occurrences)    

# drawing the static general statistics figure
stat_fig = plt.figure('Cluster size statistics ({})'.format(data_dir)) #TODO take name of the containing folder?

plt.errorbar(timesteps, avgs, yerr=avg_devs, fmt='ko-', label='average size', capsize=4, linewidth=1)
plt.errorbar(timesteps, maxs, yerr=max_devs, fmt='ro-', label='maximum size', capsize=4, linewidth=1)
plt.errorbar(timesteps, mcs, yerr=mc_devs, fmt='go-', label='most common size', capsize=4, linewidth=1)
plt.errorbar(timesteps, fms, yerr=fm_devs, fmt='bo-', label='free monomers', capsize=4, linewidth=1)

plt.legend(loc='upper left')
plt.xlabel('timestep')
plt.ylabel('number of rods', rotation='vertical')
plt.axis(ymin=0, ymax = max(maxs)*1.2)
plt.grid(axis='y')

# drawing the dynamic clutser distribution figure
dist_fig, ax = plt.subplots(num='Cluster size distribution ({})'.format(data_dir))
dist_fig.subplots_adjust(bottom=0.20)

plot_data = distributions[-1]
dist_plot = ax.errorbar(range(1, len(plot_data[0])+1), plot_data[0], yerr=plot_data[1],
                            fmt='bo-', capsize=4, linewidth=1)
avg_lines = (ax.axvline(avgs[-1], color='black', linestyle='-', linewidth=0.5),
             ax.axvline(avgs[-1]-avg_devs[-1], color='black', linestyle='--', linewidth=0.5),
             ax.axvline(avgs[-1]+avg_devs[-1], color='black', linestyle='--', linewidth=0.5))
    
ax.set_xlabel('cluster size')
ax.set_ylabel('occurrences', rotation='vertical')

timestep_slider_axes = dist_fig.add_axes([0.35, 0.02, 0.32, 0.03])
timestep_slider = Slider(timestep_slider_axes, 'timestep', 0, n_snapshots-1,
                            valinit=n_snapshots-1, valfmt="%d", valstep=1)
timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[-1])

def update_dist_plot(new_timestep_index):
    new_timestep_index = int(new_timestep_index)
    timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[new_timestep_index])
    
    plot_data = distributions[new_timestep_index]
    global dist_plot
    dist_plot.remove()
    dist_plot = ax.errorbar(range(1, len(plot_data[0])+1), plot_data[0], yerr=plot_data[1],
                                fmt='bo-', capsize=4, linewidth=1)
    ax.axis(ymax=max(plot_data[0])+max(plot_data[1]))
    
    avg = avgs[new_timestep_index]
    dev = avg_devs[new_timestep_index]
    avg_lines[0].set_xdata([avg, avg])
    avg_lines[1].set_xdata([avg-dev, avg-dev])
    avg_lines[2].set_xdata([avg+dev, avg+dev])
    dist_fig.canvas.draw_idle()
    
timestep_slider.on_changed(update_dist_plot)

plt.show()
