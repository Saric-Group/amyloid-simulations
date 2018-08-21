# encoding: utf-8
'''
Contains methods that calculate the interaction energy between rods in
different states in Eugen Rožić's spherocylinder MD/bead model
for amyloid aggregation.

Created on 12 Apr 2018

@author: Eugen Rožić
'''

import numpy as np

N = 7
r_body = 1.0
delta_body = 1.0 #in terms of r_body

r_body_sq = r_body**2

# these have to be set from outside
M = None
r_int = None
delta_int = None #in terms of r_int
bulge_out = None
int_range = None

def set_parameters(num_int_spots, int_bead_r, int_bead_overlap, int_bead_bulge, range_of_int):
    '''
    Sets the "M", "r_int" and "delta_int" global parameters
    '''
    global M, r_int, delta_int, bulge_out, int_range
    M = num_int_spots
    r_int = int_bead_r
    delta_int = int_bead_overlap
    bulge_out = int_bead_bulge
    int_range = range_of_int

def sb_interaction(int_f, r, z, phi):
    '''
    Interaction between a rod centered at (0,0) and extending along the
    z-axis and an interaction site at (r,z), which represents the center
    of the attractive end of a soluble rod.
    
    int_f : the point-point interaction function with signature (r, sigma), where
        "r" is the distance between points and "sigma" is the minimum of the interaction
    phi : the orientation of the interaction sites of the rod at (0,0) 
    '''
    # check for body overlap
    r_sq = r**2
    for i in range(N):
        z_i = (i - (N-1)/2.)*(2 - delta_body)*r_body
        if (r_sq + (z - z_i)**2 < 4*r_body_sq):
            return 0
    
    U = 0
    r_i = r_body - r_int + bulge_out
    R = r_body + r_int
    cutoff = R + int_range
    for i in range(M):
        z_i = (i - (M-1)/2.)*(2 - delta_int)*r_int
        dist = np.sqrt(r**2 + r_i**2 - 2*r_i*r*np.cos(phi) + (z-z_i)**2)
        U += int_f(dist, R, cutoff)
    return U

def bb_interaction(int_f, r, z, theta, phi, psi1 = 0, psi2 = None):
    '''
    Inter-rod beta-beta interaction. The interaction is relative to a rod
    centered at (0,0) and extending along the z-axis.
    
    int_f : the point-point interaction function with signature (r, sigma), where
        "r" is the distance between points and "sigma" is the minimum of the interaction
    theta, phi : the direction of the rod at point (r.z)
    psi1 : the direction of the patch vector of the rod at (0,0)
    psi2 : the direction of the patch vector of the rod at (r,z)
    
    NOTE: if psi2 is not given it will be set so that the "patch" on the rod at
        (r-z) points as much as possible toward the rod at (0,0)
    '''
    U = 0
    if psi2 is None:
        psi2 = np.pi - phi # turn maximally towards the patch at (0,0)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    c_psi1 = np.cos(psi1)
    s_psi1 = np.sin(psi1)
    c_psi2 = np.cos(psi2)
    s_psi2 = np.sin(psi2)
    
    # check for body overlap
    for i in range(N):
        z_i = (i - (N-1)/2.)*(2 - delta_body)*r_body
        for j in range(N):
            z_j = (j - (N-1)/2.)*(2 - delta_body)*r_body
            x_j = z_j*c_phi*s_theta
            y_j = z_j*s_phi*s_theta
            z_j = z_j*c_theta
            if ((r + x_j)**2 + (y_j)**2 + (z + z_j - z_i)**2 < 4*r_body_sq):
                return 0
    
    R = 2*r_int
    cutoff = R + int_range
    r_i = r_j = r_body - r_int + bulge_out
    x_i = c_psi1*r_i
    y_i = s_psi1*r_i
    for i in range(M):
        z_i = (i - (M-1)/2.)*(2 - delta_int)*r_int
        for j in range(M):
            z_j = (j - (M-1)/2.)*(2 - delta_int)*r_int
            x_j = r_j*(c_psi2*c_phi*c_theta - s_psi2*s_phi) + z_j*c_phi*s_theta
            y_j = r_j*(c_psi2*s_phi*c_theta + s_psi2*c_phi) + z_j*s_phi*s_theta
            z_j = -r_j*c_psi2*s_theta + z_j*c_theta
            dist = np.sqrt((r + x_j - x_i)**2 + (y_j - y_i)**2 + (z + z_j - z_i)**2)
            U += int_f(dist, R, cutoff)
    return U
