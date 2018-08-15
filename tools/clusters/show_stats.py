# encoding: utf-8
'''
Created on 21 Jun 2018

@author: Eugen Rožić

TODO
'''
import os
import argparse

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

import numpy as np
from math import sqrt, ceil
import analysis

parser = argparse.ArgumentParser(description='''
        Application for calculating the statistics of cluster data over multiple simulations''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('in_files', nargs='+', help='paths of the *_cluster_data files')
parser.add_argument('--bin_size', default=1, type=int, help='size of bins for doing cluster size distributions')
parser.add_argument('-s', '--save', default=None, type=str, help='if set saves the figures to given location with extensions "_dist.pdf" & "_stat.pdf"')
args = parser.parse_args()
data_dir = os.path.dirname(args.in_files[0])
bin_size = args.bin_size

#get data...

box_size = None
timesteps = None
n_snapshots = None
cluster_sizes_data = []
for in_file in args.in_files:
    if data_dir != os.path.dirname(in_file):
        raise Exception('All files should be from the same directory!')
    
    _box_size, _timesteps, raw_data = analysis.read_raw_data(in_file)
    
    if box_size == None:
        box_size = _box_size
    elif box_size != _box_size:
        raise Exception('All files should have the same box size!')
    if timesteps == None:
        timesteps = _timesteps
        n_snapshots = len(timesteps)
    elif len(_timesteps) != n_snapshots: #quick, non thorough check
        raise Exception('All files should have the same timesteps!')
    
    cluster_sizes_data.append(analysis.cluster_sizes_by_type(raw_data))

n_simulations = len(cluster_sizes_data)

avgs = {}
avg_devs = {}
maxs = {}
max_devs = {}
mcs = {}
mc_devs = {}
fms = {}
fm_devs = {}
distributions = {}

#this is a "hack", a fair and strong assumption really so I can skip exhaustive search...
cluster_types = set()
for sim_data in cluster_sizes_data:
    cluster_types.update(sim_data[0].keys())
    cluster_types.update(sim_data[n_snapshots/2].keys())
    cluster_types.update(sim_data[-1].keys())
    
# do statistics...

for cluster_type in cluster_types:

    avgs[cluster_type] = [0]*n_snapshots
    avg_devs[cluster_type] = [0]*n_snapshots
    maxs[cluster_type] = [0]*n_snapshots
    max_devs[cluster_type] = [0]*n_snapshots
    mcs[cluster_type] = [0]*n_snapshots
    mc_devs[cluster_type] = [0]*n_snapshots
    fms[cluster_type] = [0]*n_snapshots
    fm_devs[cluster_type] = [0]*n_snapshots
    distributions[cluster_type] = [None]*n_snapshots
    
    for i in range(n_snapshots):
        avg_sizes = np.zeros(n_simulations)
        sizes_devs = np.zeros(n_simulations)
        max_sizes = np.zeros(n_simulations, int)
        mc_sizes = np.zeros(n_simulations, int)
        free_mons = np.zeros(n_simulations, int)
        for n in range(n_simulations):
            try:
                n_sizes = np.array(cluster_sizes_data[n][i][cluster_type][0])
            except KeyError:
                continue
            n_size_occurrences = np.array(cluster_sizes_data[n][i][cluster_type][1])
            max_sizes[n] = n_sizes[-1]
            if n_sizes[0] == 1: #if smallest size is 1 (if there are single monomers)
                free_mons[n] = n_size_occurrences[0]
                n_sizes = n_sizes[1:]
                n_size_occurrences = n_size_occurrences[1:]
            if max_sizes[n] > 1: # else it should be 0, which it is by initialisation
                avg_sizes[n] = np.average(n_sizes, weights=n_size_occurrences)
                sizes_devs[n] = sqrt(np.average((n_sizes-avg_sizes[n])**2, weights=n_size_occurrences))
                mc_sizes[n] = n_sizes[np.argmax(n_size_occurrences)]
        temp = np.sum(1. / sizes_devs**2)
        avgs[cluster_type][i] = np.sum(avg_sizes / sizes_devs**2) / temp
        avg_devs[cluster_type][i] = 1. / sqrt(temp)
        maxs[cluster_type][i] = np.average(max_sizes)
        max_devs[cluster_type][i] = sqrt(np.average((max_sizes-maxs[cluster_type][i])**2))
        mcs[cluster_type][i] = np.average(mc_sizes)
        mc_devs[cluster_type][i] = sqrt(np.average((mc_sizes-mcs[cluster_type][i])**2))
        fms[cluster_type][i] = np.average(free_mons)
        fm_devs[cluster_type][i] = sqrt(np.average((free_mons-fms[cluster_type][i])**2))
        
        biggest_cluster = max(max_sizes)
        num_bins = int(ceil(float(biggest_cluster-1)/bin_size))+1
        aggregate_occurrences = [None]*num_bins
        for j in range(num_bins):
            aggregate_occurrences[j] = np.zeros(n_simulations)
        for n in range(n_simulations):
            try:
                for size, occurrences in zip(*cluster_sizes_data[n][i][cluster_type]):
                    bin_no = int((size-2)/bin_size)+1
                    if bin_no == 0:
                        aggregate_occurrences[bin_no][n] = occurrences
                    else:
                        aggregate_occurrences[bin_no][n] += occurrences/float(bin_size)
                    #all other are 0 by initialisation
            except KeyError:
                pass
        avg_occurrences = [np.average(elem) for elem in aggregate_occurrences]
        dev_occurrences = [sqrt(np.average((aggregate_occurrences[j] - avg_occurrences[j])**2)) for j in range(num_bins)]
        distributions[cluster_type][i] = (avg_occurrences, dev_occurrences)

def draw_cluster_stats(cluster_type):
    
    fig = plt.figure('Cluster type {} statistics ({})'.format(cluster_type, data_dir))

    plt.errorbar(timesteps, avgs[cluster_type], yerr=avg_devs[cluster_type], fmt='kx-', label='average size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, maxs[cluster_type], yerr=max_devs[cluster_type], fmt='rx-', label='maximum size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, mcs[cluster_type], yerr=mc_devs[cluster_type], fmt='gx-', label='most common size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, fms[cluster_type], yerr=fm_devs[cluster_type], fmt='bx-', label='free monomers', capsize=4, linewidth=0.7)

    plt.legend(loc='upper left')
    plt.xlabel('timestep')
    plt.ylabel('number of rods', rotation='vertical')
    plt.axis(ymin=0, ymax = max(maxs[cluster_type])*1.2)
    plt.grid(axis='y')
    
    return fig

def fixed_size_range(size, start=0.0, step=1.0):
    ret = [start]*size
    i = 1
    while i < size:
        ret[i] = ret[i-1] + step
        i += 1
    return ret

def draw_cluster_distributions(cluster_type):

    fig, ax = plt.subplots(num='Cluster type {} size distributions ({})'.format(cluster_type, data_dir))
    fig.subplots_adjust(bottom=0.20)

    plot_data = distributions[cluster_type][-1]
    sizes = [1]
    sizes.extend(fixed_size_range(len(plot_data[0])-1, np.average(range(2, bin_size+2)), bin_size))
    ax.errorbar(sizes, plot_data[0], yerr=plot_data[1], fmt='bo-', capsize=4, linewidth=1)
    avg = avgs[cluster_type][-1]
    dev = avg_devs[cluster_type][-1]
    avg_lines = (ax.axvline(avg, color='black', linestyle='-', linewidth=0.5),
                 ax.axvline(avg-dev, color='black', linestyle='--', linewidth=0.5),
                 ax.axvline(avg+dev, color='black', linestyle='--', linewidth=0.5))
    
    ax.set_xlabel('cluster size')
    ax.set_ylabel('occurrences', rotation='vertical')

    timestep_slider_axes = fig.add_axes([0.35, 0.02, 0.32, 0.03])
    timestep_slider = Slider(timestep_slider_axes, 'timestep', 0, n_snapshots-1,
                            valinit=n_snapshots-1, valfmt="%d", valstep=1)
    timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[-1])

    def update_dist_plot(new_timestep_index):
        new_timestep_index = int(new_timestep_index)
        timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[new_timestep_index])
    
        plot_data = distributions[cluster_type][new_timestep_index]
        sizes = fixed_size_range(len(plot_data[0]), np.average(range(1, bin_size+1)), bin_size)
        ax.containers[-1].remove()
        ax.errorbar(sizes, plot_data[0], yerr=plot_data[1], fmt='bo-', capsize=4, linewidth=1)
        ax.axis(ymax=max(plot_data[0])+max(plot_data[1]))
    
        avg = avgs[cluster_type][new_timestep_index]
        dev = avg_devs[cluster_type][new_timestep_index]
        avg_lines[0].set_xdata([avg, avg])
        avg_lines[1].set_xdata([avg-dev, avg-dev])
        avg_lines[2].set_xdata([avg+dev, avg+dev])
        fig.canvas.draw_idle()
    
    timestep_slider.on_changed(update_dist_plot)
    
    return fig, timestep_slider

widgets = []
for cluster_type in cluster_types:
    stat_fig = draw_cluster_stats(cluster_type)
    dist_fig, dist_slider = draw_cluster_distributions(cluster_type)
    widgets.append(dist_slider)
    if args.save:
        print args.save+"_"+str(cluster_type)+"_stat.pdf"
        stat_fig.savefig(args.save+"_"+str(cluster_type)+"_stat.pdf", dpi=1000, facecolor='white')
        dist_fig.savefig(args.save+"_"+str(cluster_type)+"_dist.pdf", dpi=1000, facecolor='white')

plt.show()
