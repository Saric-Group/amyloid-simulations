# encoding: utf-8
'''
Contains methods that calculate the interaction energy between rods in
different states in Anđela Šarić's spherocylinder Monte-Carlo model
for amyloid aggregation.

Created on 23 Apr 2018

@author: Eugen Rožić
'''

import numpy as np

nugget = 10**(-4)

sigma = None
L = L_half = None
L_patch = L_patch_half = None
L_patch_half_sq = None

def model_setup(model):
    '''
    model : a lammps_multistate_rods.Rod_model instance
    '''
    global sigma, sigma_sq, L, L_half, L_patch, L_patch_half, L_patch_half_sq
    sigma = 2*model.rod_radius
    sigma_sq = sigma*sigma
    L = 3.0*sigma
    L_half = L/2
    L_patch = 0.7*L
    L_patch_half = L_patch/2
    L_patch_half_sq = L_patch_half*L_patch_half

def min_patch_dist_sq(r, z, theta, phi):
    e1e2 = np.cos(theta)
    Re1 = z
    if (e1e2**2 > 1-nugget):
        e1e2 = 1 #sign not important anymore
        if (Re1 < nugget):
            return r*r
        Re2 = z
        lambda0 = Re1/2
        mu0 = -Re2/2
        d_perp_sq = r**2
    else :
        denom = (1 - e1e2**2)
        Re2 = r*np.sin(theta)*np.cos(phi) + z*e1e2
        lambda0 = (Re1 - Re2*e1e2)/denom
        mu0 = -(Re2 - Re1*e1e2)/denom
        d_perp_sq = (r*np.sin(phi))**2
    #d_perp_sq = r**2 + z**2 + mu0**2 + lambda0**2 + 2*(mu0*Re2 - lambda0*Re1 - mu0*lambda0*e1e2)
    
    if (lambda0**2 <= L_patch_half_sq and mu0**2 <= L_patch_half_sq):
        return d_perp_sq #if in-plane patches overlap
    
    gamma1 = -lambda0 - L_patch_half
    gamma2 = -lambda0 + L_patch_half
    gamma_min = gamma1 if gamma1**2 < gamma2**2 else gamma2
    delta1 = -mu0 - L_patch_half
    delta2 = -mu0 + L_patch_half
    delta_min = delta1 if delta1**2 < delta2**2 else delta2
    
    if (e1e2 == 1):
        return d_perp_sq + (delta_min - gamma_min)**2 #if parallel then it's simpler
    
    if gamma_min**2 > delta_min**2 :
        gamma = gamma_min
        delta = gamma*e1e2
        if (delta + mu0)**2 >= L_patch_half_sq :
            delta = delta1 if (delta - delta1)**2 < (delta - delta2)**2 else delta2
    else :
        delta = delta_min
        gamma = delta*e1e2
        if (gamma + lambda0)**2 >= L_patch_half_sq :
            gamma = gamma1 if (gamma - gamma1)**2 < (gamma - gamma2)**2 else gamma2
            
    #TODO do both and take smaller/bigger(?) ...
    
    d_in_plane_sq = gamma**2 + delta**2 - 2*gamma*delta*e1e2
    
    return d_perp_sq + d_in_plane_sq


def sb_interaction(r, z, phi, eps = 1.0):
    '''
    Interaction between a rod centered at (0,0) and extending along the
    z-axis and an interaction site at (r,z), which represents the center
    of the attractive cap of a soluble rod. 
    '''    
    if (phi > np.pi/2 and phi < 3*np.pi/2):
        return 0.0
    
    dist_sq = r*r
    
    z_temp = 0.0
    if z > L_half:
        z_temp = z - L_half
    elif z < -L_half:
        z_temp = z + L_half
    
    if dist_sq + z_temp**2 < sigma_sq:
        return 0.0
    
    if z > L_patch_half:
        z -= L_patch_half
        dist_sq += z*z
    elif z < -L_patch_half:
        z += L_patch_half
        dist_sq += z*z
    
    if (dist_sq > 2.25*sigma_sq or dist_sq < sigma_sq):
        return 0.0
    else:
        return -eps


def bb_interaction(r, z, theta, phi, psi1 = 0, psi2 = 0, eps = 1.0):
    '''
    Inter-rod beta-beta interaction. The interaction is relative to a rod
    centered at (0,0) and extending along the z-axis.
    
    theta, phi : the direction of the rod at point (r.z)
    psi1 : the direction of the patch vector of the rod at (0,0)
    psi2 : the direction of the patch vector of the rod at (r,z)
    
    NOTE: it depends on phi and psis only for on/off, not for strength
    
    TODO: doesn't consider the overlap of spherocylinder bodies (hard core repulsion)
    TODO: directions of patches (psi1 & psi2) are not considered
    '''
    min_dist_sq = min_patch_dist_sq(r, z, theta, phi)
    
    if (min_dist_sq > 2.25*sigma_sq or min_dist_sq < sigma_sq):
        return 0.0
    dist = np.sqrt(r*r + z*z)
    return -eps * (np.cos(theta)**2 + sigma/(dist**3))
