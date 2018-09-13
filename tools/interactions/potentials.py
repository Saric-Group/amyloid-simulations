# encoding: utf-8
'''
Contains various potentials to be (possibly) used for MD.

All accept a "sigma" parameter which is the minimum of the potential,
and a "cutoff" parameter where the potential is made 0 by shifting.

Created on 1 May 2018

@author: Eugen Rožić
'''

from math import pi, exp, sqrt, cos

def lj_n_m(n, m, r, sigma, cutoff, eps, shift = False):
    '''
    A generalised Lennard-Jones potential with repulsive exponent "n" and
    attractive exponent "m", and a minimum at (sigma, -eps)
    '''
    if r >= cutoff:
        return 0.0
    if shift: 
        co_val = m*(sigma/cutoff)**n - n*(sigma/cutoff)**m
        return eps*(m*(sigma/r)**n - n*(sigma/r)**m - co_val)/(n-m)
    else:
        return eps*(m*(sigma/r)**n - n*(sigma/r)**m)/(n-m)

def cos_sq(r, sigma, cutoff, eps):
    '''
    A potential with a cosine-squared attractive part connecting points (sigma, -eps) and (cutoff, 0) smoothly. 
    '''
    if r >= cutoff:
        return 0.0
    elif r <= sigma:
        return -eps
    else:
        return -eps*(cos(pi*(r-sigma)/(2*(cutoff-sigma))))**2

def gauss(std_dev, r, sigma, cutoff, eps, shift = False):
    '''
    A Gaussian potential with standard deviation "std_dev" centered at (sigma, -eps)
    '''
    if r >= cutoff:
        return 0.0
    if shift:
        co_val = exp(-((cutoff-sigma)/(sqrt(2)*std_dev))**2)
        return -eps*(exp(-((r-sigma)/(sqrt(2)*std_dev))**2) - co_val)
    else:
        return -eps*exp(-((r-sigma)/(sqrt(2)*std_dev))**2)

def morse(a, r, sigma, cutoff, eps, shift = False):
    '''
    A Morse potential with a minimum at (sigma, -eps) with "width" given by "a".
    '''
    if r >= cutoff:
        return 0.0
    if shift:
        co_val = exp(-2*a*(cutoff-sigma)) - 2*exp(-a*(cutoff-sigma))
        return eps*(exp(-2*a*(r-sigma)) - 2*exp(-a*(r-sigma)) - co_val)
    else:
        return eps*(exp(-2*a*(r-sigma)) - 2*exp(-a*(r-sigma)))
