# encoding: utf-8
'''
TODO

Created on 21 Jun 2018

@author: Eugen Rožić
'''
import os
import argparse

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

import numpy as np
from math import sqrt
from lammps_multistate_rods.tools import cluster_analysis

def fetch_data(in_files):
    box_size = None
    timesteps = None
    n_snapshots = None
    cluster_sizes_data = []
    max_cluster_size = 0
    for in_file in in_files:
        if data_dir != os.path.dirname(in_file):
            raise Exception('All files should be from the same directory!')
    
        _box_size, _timesteps, raw_data = cluster_analysis.read_raw_data(in_file)
    
        if box_size == None:
            box_size = _box_size
        elif box_size != _box_size:
            raise Exception('All files should have the same box size!')
        if timesteps == None:
            timesteps = _timesteps
            n_snapshots = len(timesteps)
        elif len(_timesteps) != n_snapshots: #quick, non thorough check
            raise Exception('All files should have the same timesteps!')
    
        file_cluster_data, _max_cluster_size = cluster_analysis.cluster_sizes_by_type(raw_data) 
        cluster_sizes_data.append(file_cluster_data)
        if _max_cluster_size > max_cluster_size:
            max_cluster_size = _max_cluster_size
    
    return box_size, timesteps, cluster_sizes_data, max_cluster_size

def do_overall_stats():
    
    avgs = [0]*n_snapshots
    avg_devs = [0]*n_snapshots
    maxs = [0]*n_snapshots
    max_devs = [0]*n_snapshots
    mcs = [0]*n_snapshots
    mc_devs = [0]*n_snapshots
    fms = [0]*n_snapshots
    fm_devs = [0]*n_snapshots
    
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
        avgs[i] = np.sum(avg_sizes / sizes_devs**2) / temp
        avg_devs[i] = 1. / sqrt(temp)
        maxs[i] = np.average(max_sizes)
        max_devs[i] = sqrt(np.average((max_sizes-maxs[i])**2))
        mcs[i] = np.average(mc_sizes)
        mc_devs[i] = sqrt(np.average((mc_sizes-mcs[i])**2))
        fms[i] = np.average(free_mons)
        fm_devs[i] = sqrt(np.average((free_mons-fms[i])**2))
    
    return (avgs, avg_devs, maxs, max_devs, mcs, mc_devs, fms, fm_devs)

def generate_pdfs(last_n):
    
    aggregate_occurrences = [None]*last_n
    
    offset = n_snapshots - last_n
    for i in range(last_n):
        aggregate_occurrences[i] = [None]*max_cluster_size
        for j in range(max_cluster_size):
            aggregate_occurrences[i][j] = np.zeros(n_simulations)
        for n in range(n_simulations):
            try:
                data = cluster_sizes_data[n][offset+i][cluster_type]
            except KeyError:
                continue
            for size, occurrences in zip(*data):
                aggregate_occurrences[i][size-1][n] = occurrences        
                
    return aggregate_occurrences

def plot_overall_stats(plot_data):
    
    #TODO make this a subplot of a figure...
    fig = plt.figure('Cluster type {} overall statistics ({})'.format(cluster_type, data_dir))

    plt.errorbar(timesteps, plot_data[0], yerr=plot_data[1], fmt='kx-',
                 label='average size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, plot_data[2], yerr=plot_data[3], fmt='rx-',
                 label='maximum size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, plot_data[4], yerr=plot_data[5], fmt='gx-',
                 label='most common size', capsize=4, linewidth=0.7)
    plt.errorbar(timesteps, plot_data[6], yerr=plot_data[7], fmt='bx-',
                 label='free monomers', capsize=4, linewidth=0.7)

    plt.legend(loc='upper left')
    plt.xlabel('timestep')
    plt.ylabel('number of rods', rotation='vertical')
    plt.axis(ymin=0, ymax = max(plot_data[2])*1.2)
    plt.grid(axis='y')
    
    return fig

def plot_interactive_pdf(pdf_data):
    
    distributions = [None]*len(pdf_data)
    for i in range(len(pdf_data)):
        avg_occurrences = [np.average(occurrences) for occurrences in pdf_data[i]]
        dev_occurrences = [sqrt(np.average((pdf_data[i][j] - avg_occurrences[j])**2))
                           for j in range(max_cluster_size)]
        distributions[i] = (avg_occurrences, dev_occurrences)
        
    #TODO make this just a subplot of a figure... (how, with the slider?!?!)
    fig, ax = plt.subplots(num='Cluster type {} size distributions ({})'.format(cluster_type, data_dir))
    fig.subplots_adjust(bottom=0.20)

    plot_data = distributions[-1]
    sizes = range(1, max_cluster_size+1)
    ax.errorbar(sizes, plot_data[0], yerr=plot_data[1], fmt='bo-', capsize=4, linewidth=1)
    
    ax.set_xlabel('cluster size')
    ax.set_ylabel('occurrences', rotation='vertical')

    timestep_slider_axes = fig.add_axes([0.35, 0.02, 0.32, 0.03])
    timestep_slider = Slider(timestep_slider_axes, 'timestep', 0, n_snapshots-1,
                            valinit=n_snapshots-1, valfmt="%d", valstep=1)
    timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[-1])

    def update_dist_plot(new_timestep_index):
        new_timestep_index = int(new_timestep_index)
        timestep_slider.valtext.set_text(timestep_slider.valfmt % timesteps[new_timestep_index])
    
        plot_data = distributions[new_timestep_index]
        ax.containers[-1].remove()
        ax.errorbar(sizes, plot_data[0], yerr=plot_data[1], fmt='bo-', capsize=4, linewidth=1)

        fig.canvas.draw_idle()
    
    timestep_slider.on_changed(update_dist_plot)
    
    return fig, timestep_slider

def plot_aggregate_pdf(pdf_data):
    
    pdf_data = [np.array([item for sublist in snapshot for item in sublist])
                for snapshot in zip(*pdf_data)]
    avg_occurrences = [np.average(occurrences) for occurrences in pdf_data]
    dev_occurrences = [sqrt(np.average((pdf_data[j] - avg_occurrences[j])**2))
                       for j in range(max_cluster_size)]
        
    #TODO make this a subplot of a figure...
    fig = plt.figure('Cluster type {} aggregate size distributions ({})'.format(cluster_type, data_dir))

    sizes = range(1, max_cluster_size+1)
    plt.errorbar(sizes, avg_occurrences, yerr=dev_occurrences, fmt='bo-', capsize=4, linewidth=1)
    
    plt.xlabel('cluster size')
    plt.ylabel('occurrences', rotation='vertical')
    
    return fig


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                description='Application for calculating the statistics of cluster data over \
multiple simulations',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_files', nargs='+', 
                        help='paths of the *_cluster_data files')
    parser.add_argument('--overall', action='store_true',
                        help='whether to plot the overall cluster statistics')
    parser.add_argument('--pdf', type=str, 
                        help='"i" to plot an interactive cluster PDF, or an int <N> to plot a \
static cluster PDF averaged over the last N timesteps')
    parser.add_argument('-s', '--save', type=str,
                        help='saves the figures to the given directory')
    args = parser.parse_args()

    data_dir = os.path.dirname(args.in_files[0])
    box_size, timesteps, cluster_sizes_data, max_cluster_size = fetch_data(args.in_files)
    n_snapshots = len(timesteps)
    n_simulations = len(cluster_sizes_data)
    
    #this is a "hack", but based on a fair and strong assumption really, so I can skip exhaustive search...
    cluster_types = set()
    for sim_data in cluster_sizes_data:
        cluster_types.update(sim_data[0].keys())
        cluster_types.update(sim_data[n_snapshots/2].keys())
        cluster_types.update(sim_data[-1].keys())        
    
    overall_figure = None
    pdf_figure = None
    widgets = []
    for cluster_type in cluster_types:
        
        if args.overall:
            if not overall_figure:
                overall_figure = None #TODO create the figure with len(cluster_types) subplots
            
            overall_stats = do_overall_stats()
            overall_plot = plot_overall_stats(overall_stats)
            
            #TODO put the "overall_plot" in the "overall_figure" in the right spot...
        
        if args.pdf:
            if not pdf_figure:
                pdf_figure = None #TODO create the figure with len(cluster_types) subplots
            
            if args.pdf == 'i':
                pdf_data = generate_pdfs(n_snapshots)
                pdf_plot, slider = plot_interactive_pdf(pdf_data)
                widgets.append(slider)
            else:
                try:
                    last_n = int(args.pdf)
                    if last_n > n_snapshots:
                        raise Exception()
                except:
                    print 'ERROR: nonsupported value ({}) for the "--pdf" option! (only integers smaller \
than the number of snapshots, or "i" for interactive mode, accepted)'.format(args.pdf)
                    quit()
                pdf_data = generate_pdfs(last_n)
                pdf_plot = plot_aggregate_pdf(pdf_data)
            
            #TODO put the "pdf_plot" in the "pdf_figure" in the right spot...
    
    if args.save:
        if overall_figure:
            overall_figure.savefig(os.path.join(args.save,"overall_cluster_stats.pdf"), dpi=1000, facecolor='white')
        if pdf_figure:
            pdf_figure.savefig(os.path.join(args.save,"cluster_distribution.pdf"), dpi=1000, facecolor='white')
            
    plt.show()

