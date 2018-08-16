# encoding: utf-8
'''
Created on 15 Aug 2018

@author: Eugen Rožić

TODO
'''
import os
import argparse

import numpy as np
from matplotlib import pyplot as plt

import analysis

def concentration(number, volume):
    '''
    Concentration in M for volume in nm^3
    '''
    return (10.0*number)/(6.022*volume)

parser = argparse.ArgumentParser(description='''
        Application for plotting the free monomer concentrations''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('in_files', nargs='+', help='paths of the *_cluster_data files')
parser.add_argument('-s', '--save', default=None, type=str, help='saves the figure to the given location')
args = parser.parse_args()

data_dir = os.path.dirname(args.in_files[0])
data_len = len(args.in_files)
avg_over = 10 #TODO make this a parameter to set...

xs = np.zeros(data_len)
ys = np.zeros(data_len)

for n in range(data_len):
    if data_dir != os.path.dirname(args.in_files[n]):
        raise Exception('All files should be from the same directory!') #??
    
    box_size, timesteps, raw_data = analysis.read_raw_data(args.in_files[n])
    volume = reduce(lambda x,y: x*y, box_size)
    
    for i in range(avg_over):
        free_monomers, total_monomers = analysis.free_monomers(raw_data[-1-i])
        if xs[n] == 0:
            xs[n] = total_monomers 
        elif total_monomers != xs[n]:
            raise Exception('Total number of monomers is not constant!') #??
        ys[n] += free_monomers
    ys[n] = ys[n] / avg_over
    
    xs[n] = concentration(xs[n], volume)
    ys[n] = concentration(ys[n], volume)
    
fig = plt.figure(data_dir)

plt.plot(np.log(xs), np.log(ys), 'bo')
    
plt.xlabel(r'$\ln\left(c_{total}\;/\;M\right)$')
plt.ylabel(r'$\ln\left(c_{free}\;/\;M\right)$', rotation='vertical')
plt.axis(ymin=-11, ymax=-5.5, xmin=-11, xmax=-5.5)
plt.grid(True)

if args.save:
    fig.savefig(args.save, dpi=1000, facecolor='white')

plt.show()
    
    