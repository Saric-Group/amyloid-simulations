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
L_patch = None
aux1 = None
aux2 = None

def setup(model):
    '''
    model : a lammps_multistate_rods.Model instance
    '''
    global sigma, L_patch, aux1, aux2
    sigma = 2*model.rod_radius
    L_patch = 0.7*3*sigma # 4*sigma is the length of the rod
    aux1 = L_patch/2
    aux2 = aux1**2

def min_patch_dist(r, z, theta, phi):
    e1e2 = np.cos(theta)
    Re1 = z
    if (e1e2**2 > 1-nugget):
        e1e2 = 1 #sign not important anymore
        if (Re1 < nugget):
            return np.sqrt(r**2)
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
    
    if (lambda0**2 <= aux2 and mu0**2 <= aux2):
        return np.sqrt(d_perp_sq) #if in-plane patches overlap
    
    gamma1 = -lambda0 - aux1
    gamma2 = -lambda0 + aux1
    gamma_min = gamma1 if gamma1**2 < gamma2**2 else gamma2
    delta1 = -mu0 - aux1
    delta2 = -mu0 + aux1
    delta_min = delta1 if delta1**2 < delta2**2 else delta2
    
    if (e1e2 == 1):
        return np.sqrt(d_perp_sq + (delta_min - gamma_min)**2) #if parallel then it's simpler
    
    if gamma_min**2 > delta_min**2 :
        gamma = gamma_min
        delta = gamma*e1e2
        if (delta + mu0)**2 >= aux2 :
            delta = delta1 if (delta - delta1)**2 < (delta - delta2)**2 else delta2
    else :
        delta = delta_min
        gamma = delta*e1e2
        if (gamma + lambda0)**2 >= aux2 :
            gamma = gamma1 if (gamma - gamma1)**2 < (gamma - gamma2)**2 else gamma2
            
    #TODO do both and take smaller/bigger(?) ...
    
    d_in_plane_sq = gamma**2 + delta**2 - 2*gamma*delta*e1e2
    
    return np.sqrt(d_perp_sq + d_in_plane_sq)


def sb_interaction(r, z, phi, eps = 1):
    '''
    Interaction between a rod centered at (0,0) and extending along the
    z-axis and an interaction site at (r,z), which represents the center
    of the attractive cap of a soluble rod.
    
    TODO: doesn't consider the overlap of spherocylinder bodies (hard core repulsion) 
    '''    
    if (phi > np.pi/2 and phi < 3*np.pi/2):
        return 0
    
    dist = r*r
    if z > L_patch/2 :
        z = z - L_patch/2
        dist += z*z
    elif z < -L_patch/2 :
        z = z + L_patch/2
        dist += z*z
        
    dist = np.sqrt(dist)
    
    if (dist > 1.5*sigma or dist < sigma):
        return 0
    else:
        return -eps


def bb_interaction(r, z, theta, phi, psi1 = 0, psi2 = 0, eps = 1):
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
    min_dist = min_patch_dist(r, z, theta, phi)
    
    if (min_dist > 1.5*sigma or min_dist < sigma):
        return 0
    dist = np.sqrt(r*r + z*z)
    return -eps * (np.cos(theta)**2 + sigma/(dist**3))
