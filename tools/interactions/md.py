# encoding: utf-8
'''
Contains methods that calculate the interaction energy between rods in
different states in Eugen Rožić's spherocylinder MD/bead model
for amyloid aggregation.

Created on 12 Apr 2018

@author: Eugen Rožić
'''

import potentials
import numpy as np

N = None
r_body = None
r_body_sq = None
delta_body = None #in terms of r_body

M = None
r_int = None
delta_int = None #in terms of r_int
bulge_out = None

sol_active = None
beta_active = None

int_potentials = {}

def int_f(int_type, R, eps):
    if int_type[0] == 'lj/cut':
        return lambda r: potentials.lj_n_m(12, 6, r, R, R + int_type[1], eps)
    elif int_type[0] == 'cosine/squared':
        return lambda r: potentials.cos_sq(r, R, R + int_type[1], eps)
    elif int_type[0] == 'nm/cut':
        return lambda r: potentials.lj_n_m(int_type[1], int_type[2], r, R, R + int_type[3], eps)
    elif int_type[0] == 'morse':
        return lambda r: potentials.morse(int_type[1], r, R, R + int_type[2], eps)
    elif int_type[0] == 'gauss/cut':
        return lambda r: potentials.gauss(int_type[1], r, R, R + int_type[2], eps)
    else:
        raise Exception('Unknown/invalid int_type parameter: '+ str(int_type))

def setup(model):
    '''
    model : a lammps_multistate_rods.Model instance
    sb : a list of "model.eps" entries/keys that constitute soluble-beta interaction
    bb : a list of "model.eps" entries/keys that constitute beta-beta interaction
    '''
    global sol_active, beta_active
    sol_active = filter(lambda x: x in model.active_bead_types,
           [int(bead_type) for bead_type in model.state_structures[0].replace('|','')])[0]
    beta_active= filter(lambda x: x in model.active_bead_types,
           [int(bead_type) for bead_type in model.state_structures[1].replace('|','')])
    
    global N, r_body, r_body_sq, delta_body
    N = model.body_beads
    r_body = model.rod_radius
    r_body_sq = r_body**2
    delta_body = model.body_bead_overlap / r_body
    
    global M, r_int, delta_int, bulge_out
    M = model.int_sites
    r_int = model.int_radius
    delta_int = model.int_bead_overlap / r_int
    bulge_out = model.int_bulge_out
    
    beta_types = sorted(set(beta_active))
    for i in range(len(beta_types)):
        R = r_body + r_int
        (eps, int_type_key) = model.eps[(sol_active, beta_types[i])]
        int_type = model.int_types[int_type_key]
        
        int_potentials[(sol_active, beta_types[i])] = int_f(int_type, R, eps)
        
        for j in range(i, len(beta_types)):
            R = 2*r_int
            (eps, int_type_key) = model.eps[(beta_types[i], beta_types[j])]
            int_type = model.int_types[int_type_key]
        
            int_potentials[(beta_types[i], beta_types[j])] = \
                int_potentials[(beta_types[j], beta_types[i])] = int_f(int_type, R, eps)
    

def sb_interaction(r, z, phi):
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
    for i in range(M):
        z_i = (i - (M-1)/2.)*(2 - delta_int)*r_int
        dist = np.sqrt(r**2 + r_i**2 - 2*r_i*r*np.cos(phi) + (z-z_i)**2)
        U += int_potentials[(sol_active, beta_active[i])](dist)
    return U

def bb_interaction(r, z, theta, phi, psi1 = 0, psi2 = None):
    '''
    Inter-rod beta-beta interaction. The interaction is relative to a rod
    centered at (0,0) and extending along the z-axis.
    
    int_f : the point-point interaction function with signature (r, sigma), where
        "r" is the distance between points and "sigma" is the minimum of the interaction
    theta, phi : the direction of the rod at point (r,z)
    psi1 : the direction of the patch vector of the rod at (0,0)
    psi2 : the direction of the patch vector of the rod at (r,z)
    
    NOTE: if psi2 is not given it will be set so that the "patch" of the rod at
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
            U += int_potentials[(beta_active[i], beta_active[j])](dist)
    return U
