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
import lammps_multistate_rods

model = None

N = None # num of body beads
K = None # num of patches
M = None # nums of patch beads by patch

r_rod_sq = None

body_z = None
patch_z = None
patch_r = None
patch_phi = None

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

def setup(rod_model):
    '''
    rod_model : a lammps_multistate_rods.Model instance
    '''
    global model
    model = rod_model
    
    global N, K, M, r_rod_sq
    N = model.body_beads
    K = model.num_patches
    M = model.patch_beads
    r_rod_sq = model.rod_radius**2
    
    global body_z, patch_z, patch_r, patch_phi
    body_z = [(i - (N-1)/2.)*(2*model.rod_radius - model.body_bead_overlap)
              for i in range(N)]
    patch_z = [[(i - (M[k]-1)/2.)*(2*model.patch_bead_radii[k] + model.patch_bead_sep[k])
                for i in range(M[k])] for k in range(K)]
    patch_r = [model.rod_radius - model.patch_bead_radii[k] + model.patch_bulge_out[k]
               for k in range(K)]
    patch_phi = [np.deg2rad(model.patch_angles[k])
                 for k in range(K)]
    
    global sol_active, beta_active
    # assumption that there is only one active body bead type (the tip)
    sol_active = filter(lambda x: x in model.active_bead_types, model.body_bead_types)[0]
    # assumption that only patch beads are active, and all of them
    beta_active = [filter(lambda x: x in model.active_bead_types, patch)
                   for patch in model.state_structures[1][1:]]
    
    for k1 in range(K):
        k1_patch_types = sorted(set(beta_active[k1]))
        for i in range(len(k1_patch_types)):
            R = model.rod_radius + model.patch_bead_radii[k1]
            try:
                (eps, int_type_key) = model.eps[tuple(sorted([sol_active, k1_patch_types[i]]))]
            except:
                (eps, int_type_key) = (0.0, lammps_multistate_rods.model.vx)
            int_type = model.int_types[int_type_key]
        
            int_potentials[(sol_active, k1_patch_types[i])] = \
            int_potentials[(k1_patch_types[i], sol_active)] = int_f(int_type, R, eps)
        
            for j in range(i, len(k1_patch_types)):
                R = 2*model.patch_bead_radii[k1]
                try:
                    (eps, int_type_key) = model.eps[(k1_patch_types[i], k1_patch_types[j])]
                except:
                    (eps, int_type_key) = (0.0, lammps_multistate_rods.model.vx)
                int_type = model.int_types[int_type_key]
        
                int_potentials[(k1_patch_types[i], k1_patch_types[j])] = \
                int_potentials[(k1_patch_types[j], k1_patch_types[i])] = int_f(int_type, R, eps)
            
            for k2 in range(k1+1, K):
                k2_patch_types = list(set(beta_active[k2]))
                for j in range(len(k2_patch_types)):
                    R = model.patch_bead_radii[k1] + model.patch_bead_radii[k2]
                    try:
                        (eps, int_type_key) = model.eps[tuple(sorted([k1_patch_types[i], k2_patch_types[j]]))]
                    except:
                        (eps, int_type_key) = (0.0, lammps_multistate_rods.model.vx)
                    int_type = model.int_types[int_type_key]
        
                    int_potentials[(k1_patch_types[i], k2_patch_types[j])] = \
                    int_potentials[(k2_patch_types[j], k1_patch_types[i])] = int_f(int_type, R, eps)


def point_rod_int(bead_type, rod_state, r, z, phi):
    '''
    Interaction between a rod centered at (0,0) and extending along the
    z-axis and an interaction site at (r,z), which represents the center
    of the attractive end of a soluble rod.
    
    int_f : the point-point interaction function with signature (r, sigma), where
        "r" is the distance between points and "sigma" is the minimum of the interaction
    phi : the orientation of the rod at (0,0) 
    '''
    # check for body overlap
    r_sq = r**2
    for i in range(N):
        if (r_sq + (z - body_z[i])**2 < 4*r_rod_sq):
            return 0
    
    U = 0
    for k in range(K):
        r_i = patch_r[k]
        phi_i = phi + patch_phi[k]
        ort_part = r**2 + r_i**2 - 2*r_i*r*np.cos(phi_i)
        for i in range(M[k]):
            z_i = patch_z[k][i]
            dist = np.sqrt(ort_part + (z - z_i)**2)
            U += int_potentials[(sol_active, beta_active[k][i])](dist)
    return U

def rod_rod_int(rod1_state, rod2_state, r, z, theta, phi, psi1, psi2):
    '''
    Inter-rod beta-beta interaction. The interaction is relative to a rod
    centered at (0,0) and extending along the z-axis.
    
    int_f : the point-point interaction function with signature (r, sigma), where
        "r" is the distance between points and "sigma" is the minimum of the interaction
    theta, phi : the direction of the rod at point (r,z)
    psi1 : the direction of the patch vector of the rod at (0,0)
    psi2 : the direction of the patch vector of the rod at (r,z); the patch starts facing
        down and turns in the opposite direction to psi1 (as if the rods are mirror images)
    '''
    if theta < 0 or theta > np.pi:
        raise Exception("Theta has to be in interval [0, pi]")
    elif theta < np.pi/2:
        psi2 = np.pi - phi - psi2
    else:
        psi2 = phi + psi2
    
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    
    # check for body overlap
    for i in range(N):
        for j in range(N):
            temp = body_z[j]*s_theta
            x_j = c_phi*temp
            y_j = s_phi*temp
            z_j = body_z[j]*c_theta
            if ((r + x_j)**2 + (y_j)**2 + (z + z_j - body_z[i])**2 < 4*r_rod_sq):
                return 0
    
    x1 = [0.0]*K
    y1 = [0.0]*K
    z1 = patch_z
    x2 = [None]*K
    y2 = [None]*K
    z2 = [None]*K
    for k in range(K):
        psi1_k = psi1 + patch_phi[k]
        c_psi1 = np.cos(psi1_k)
        s_psi1 = np.sin(psi1_k)
        psi2_k = psi2 + patch_phi[k]
        c_psi2 = np.cos(psi2_k)
        s_psi2 = np.sin(psi2_k)
        r1 = r2 = patch_r[k]
        x1[k] = r1*c_psi1
        y1[k] = r1*s_psi1
        x2[k] = [0.0]*M[k]
        y2[k] = [0.0]*M[k]
        z2[k] = [0.0]*M[k]
        for i in range(M[k]):
            x2[k][i] = r2*(c_psi2*c_phi*c_theta - s_psi2*s_phi) + patch_z[k][i]*c_phi*s_theta
            y2[k][i] = r2*(c_psi2*s_phi*c_theta + s_psi2*c_phi) + patch_z[k][i]*s_phi*s_theta
            z2[k][i] = -r2*c_psi2*s_theta + patch_z[k][i]*c_theta
    U = 0
    for k1 in range(K):
        for i in range(M[k1]):        
            for k2 in range(K):
                for j in range(M[k2]):
                    dist = np.sqrt((r + x2[k2][j] - x1[k1])**2 + 
                                   (y2[k2][j] - y1[k1])**2 +
                                   (z + z2[k2][j] - z1[k1][i])**2)
                    U += int_potentials[(beta_active[k1][i], beta_active[k2][j])](dist)
    return U
