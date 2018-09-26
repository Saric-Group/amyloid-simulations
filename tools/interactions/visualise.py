# encoding: utf-8
'''
This can draw interactive plots of MC and MD potentials (and their difference),
defined in "mc_interactions" and "md_interactions", with different possible
interaction potentials for MD.

Created on 23 Apr 2018

@author: Eugen Rožić
'''
import os

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

import numpy as np

from lammps_multistate_rods import Model
import mc
import md

def draw_models(axes = None, r = 0, z = 0, phi = 0, old = []):
    
    if axes is None :
        axes = plt.gca()
    r_body = model.rod_radius
    r_int = model.patch_bead_radii
    
    for i in range(len(old),0,-1):
        old[i-1].remove()
        del old[i-1]
            
    new = old
    
    # the MD model body
    for i in range(md.N): 
        new.append(
            axes.add_patch(plt.Circle((z + md.body_z[i], r), 2*r_body, fill=True, color='#DDDDDD', alpha=0.7)))
    for i in range(md.N):
        new.append(
            axes.add_patch(plt.Circle((z + md.body_z[i], r), r_body, fill=True, color='#777777')))
        
    # the MC model patch axis
    new.append(
        axes.add_patch(plt.Rectangle((r - mc.L_patch/2, z - 0.02), mc.L_patch, 0.04, color='black')))
    
    # the MD model interaction centers (patches)
    for k in range(md.K):
        r_k = r + md.patch_r[k]*np.cos(md.patch_phi[k] + phi)
        for i in range(md.M[k]):
            new.append(
                axes.add_patch(plt.Circle((z + md.patch_z[k][i], r_k), r_int[k], fill=True, color='red')))
    
    return new

def draw_sb_2D(sb_vals, md_sb_min, prefix = ''):
    
    fig_2D, ax_2D = plt.subplots(num=prefix+" soluble-beta interaction (side 2D)", figsize=(11,7))
    fig_2D.subplots_adjust(left=0.13, bottom=0.13, right=0.99, top=0.99) #0.86 x 0.86
    
    img = ax_2D.imshow(sb_vals['md'][0], extent=[zmin, zmax, rmin, rmax],
                       vmin=md_sb_min, vmax=-md_sb_min,
                       origin="lower", interpolation='bilinear', cmap=cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    phi_slider_axes = fig_2D.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    data_chooser_axes = fig_2D.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_img(_):
        curr_theta_index = theta_index(phi_slider.val)
        curr_data = data_chooser.value_selected
        draw_models(ax_2D, phi = phi_slider.val, old = model_patches)
        img.set_data(sb_vals[curr_data][curr_theta_index])
        fig_2D.canvas.draw_idle()
    
    phi_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
    
    return phi_slider, data_chooser

def draw_sb_z_slice(sb_vals, md_sb_min, prefix = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=prefix+" soluble-beta interaction (z-slice)", figsize=(9,5))
    fig_r.subplots_adjust(left=0.18, bottom=0.2, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = ax_r.plot(rs, sb_vals['md'][0].T[zero_z], 'r-', lw=1.0)
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.axis([rmin, rmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    phi_slider_axes = fig_r.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    z_slider_axes = fig_r.add_axes([0.41, 0.06, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, 'z', zs[0], zs[-1], valinit=zs[zero_z])
    
    data_chooser_axes = fig_r.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_plot(_):
        curr_z_index = z_index(z_slider.val)
        curr_theta_index = theta_index(phi_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    z_slider.on_changed(update_plot)
    phi_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    
    return z_slider, phi_slider, data_chooser

def draw_sb_r_slice(sb_vals, md_sb_min, prefix = ''):
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=prefix+" soluble-beta interaction (r-slice)", figsize=(11,5))
    fig_z.subplots_adjust(left=0.18, bottom=0.2, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = ax_z.plot(zs, sb_vals['md'][0][zero_r], 'r-', lw=1.0)
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.axis([zmin, zmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    phi_slider_axes = fig_z.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    r_slider_axes = fig_z.add_axes([0.41, 0.06, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, 'r', rs[0], rs[-1], valinit=rs[zero_r])
    
    data_chooser_axes = fig_z.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_plot(_):
        curr_r_index = r_index(r_slider.val)
        curr_theta_index = theta_index(phi_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index][curr_r_index])
        fig_z.canvas.draw_idle()
    
    r_slider.on_changed(update_plot)
    phi_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    
    return r_slider, phi_slider, data_chooser


def draw_bb_2D(bb_vals, md_bb_min, prefix = ''):
    
    fig_2D, ax_2D = plt.subplots(num=prefix+" beta-beta interaction", figsize=(10,7))
    fig_2D.subplots_adjust(left=0.14, bottom=0.24, right=0.99, top=0.99) #0.88 x 0.88
    
    img = ax_2D.imshow(bb_vals['md'][0][0][0],
                    extent=[zmin, zmax, rmin, rmax], vmin=md_bb_min, vmax=-md_bb_min,
                    origin="lower", interpolation='bilinear', cmap=cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    theta_slider = phi_slider = psi_slider = None
    
    if len(bb_vals['md']) > 1:
        psi_slider_axes = fig_2D.add_axes([0.38, 0.10, 0.35, 0.03], facecolor=widget_color)
        psi_slider = Slider(psi_slider_axes, 'Psi', thetas[0], thetas[-1], valinit=thetas[0])
        
    if len(bb_vals['md'][0]) > 1:
        phi_slider_axes = fig_2D.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
        phi_slider = Slider(phi_slider_axes, 'Phi', phis[0], phis[-1], valinit=phis[0])
    
    if len(bb_vals['md'][0][0]) > 1:
        theta_slider_axes = fig_2D.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        theta_slider = Slider(theta_slider_axes, 'Theta', thetas[0], thetas[-1], valinit=thetas[0])
    
    data_chooser_axes = fig_2D.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_img(_):
        curr_phi_index = curr_psi_index = curr_theta_index = 0
        if psi_slider != None:
            curr_psi_index = theta_index(psi_slider.val)
            draw_models(ax_2D, phi = psi_slider.val, old = model_patches)
        if phi_slider != None:
            curr_phi_index = phi_index(phi_slider.val)
        if theta_slider != None:
            curr_theta_index = theta_index(theta_slider.val)
        curr_data = data_chooser.value_selected
        
        img.set_data(bb_vals[curr_data][curr_psi_index][curr_phi_index][curr_theta_index])
        fig_2D.canvas.draw_idle()
    
    if theta_slider != None:
        theta_slider.on_changed(update_img)
    if phi_slider != None:
        phi_slider.on_changed(update_img)
    if psi_slider != None:
        psi_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
        
    return theta_slider, phi_slider, psi_slider, data_chooser

def draw_bb_z_slice(bb_vals, md_bb_min, prefix = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=prefix+" beta-beta interaction (z-slice)", figsize=(9,6))
    fig_r.subplots_adjust(left=0.14, bottom=0.28, right=0.99, top=0.99) #0.86 x 0.86
    
    lines = ax_r.plot(rs, bb_vals['md'][0][0][0].T[zero_z], 'r-', lw=1.0)
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.axis([rmin, rmax, 1.05*md_bb_min, -1.05*md_bb_min])
    
    theta_slider = phi_slider = psi_slider = None
    
    if len(bb_vals['md']) > 1:
        psi_slider_axes = fig_r.add_axes([0.38, 0.10, 0.35, 0.03], facecolor=widget_color)
        psi_slider = Slider(psi_slider_axes, 'Psi', thetas[0], thetas[-1], valinit=thetas[0])
        
    if len(bb_vals['md'][0]) > 1:
        phi_slider_axes = fig_r.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
        phi_slider = Slider(phi_slider_axes, 'Phi', phis[0], phis[-1], valinit=phis[0])
    
    if len(bb_vals['md'][0][0]) > 1:
        theta_slider_axes = fig_r.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        theta_slider = Slider(theta_slider_axes, 'Theta', thetas[0], thetas[-1], valinit=thetas[0])
    
    z_slider_axes = fig_r.add_axes([0.38, 0.14, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, 'z', zs[0], zs[-1], valinit=zs[zero_z])
    
    data_chooser_axes = fig_r.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_plot(_):
        curr_phi_index = curr_psi_index = curr_theta_index = 0
        if psi_slider != None:
            curr_psi_index = theta_index(psi_slider.val)
        if phi_slider != None:
            curr_phi_index = phi_index(phi_slider.val)
        if theta_slider != None:
            curr_theta_index = theta_index(theta_slider.val)
        curr_z_index = z_index(z_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(bb_vals[curr_data][curr_psi_index][curr_phi_index][curr_theta_index].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    if psi_slider != None:
        psi_slider.on_changed(update_plot)
    if phi_slider != None:
        phi_slider.on_changed(update_plot)
    if theta_slider != None:
        theta_slider.on_changed(update_plot)
    z_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
        
    return psi_slider, phi_slider, theta_slider, z_slider, data_chooser

def draw_bb_r_slice(bb_vals, md_bb_min, prefix = ''):
    # 1D r-E plot
    fig_z, ax_z = plt.subplots(num=prefix+" beta-beta interaction (r-slice)", figsize=(9,6))
    fig_z.subplots_adjust(left=0.14, bottom=0.28, right=0.99, top=0.99) #0.86 x 0.86
    
    lines = ax_z.plot(zs, bb_vals['md'][0][0][0][zero_r], 'r-', lw=1.0)
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.axis([zmin, zmax, 1.05*md_bb_min, -1.05*md_bb_min])
    
    theta_slider = phi_slider = psi_slider = None
    
    if len(bb_vals['md']) > 1:
        psi_slider_axes = fig_z.add_axes([0.38, 0.10, 0.35, 0.03], facecolor=widget_color)
        psi_slider = Slider(psi_slider_axes, 'Psi', thetas[0], thetas[-1], valinit=thetas[0])
        
    if len(bb_vals['md'][0]) > 1:
        phi_slider_axes = fig_z.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
        phi_slider = Slider(phi_slider_axes, 'Phi', phis[0], phis[-1], valinit=phis[0])
    
    if len(bb_vals['md'][0][0]) > 1:
        theta_slider_axes = fig_z.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        theta_slider = Slider(theta_slider_axes, 'Theta', thetas[0], thetas[-1], valinit=thetas[0])
    
    r_slider_axes = fig_z.add_axes([0.38, 0.14, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, 'r', rs[0], rs[-1], valinit=rs[zero_z])
    
    data_chooser_axes = fig_z.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_plot(_):
        curr_phi_index = curr_psi_index = curr_theta_index = 0
        if psi_slider != None:
            curr_psi_index = theta_index(psi_slider.val)
        if phi_slider != None:
            curr_phi_index = phi_index(phi_slider.val)
        if theta_slider != None:
            curr_theta_index = theta_index(theta_slider.val)
        curr_r_index = r_index(r_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(bb_vals[curr_data][curr_psi_index][curr_phi_index][curr_theta_index][curr_r_index])
        fig_z.canvas.draw_idle()
    
    if psi_slider != None:
        psi_slider.on_changed(update_plot)
    if phi_slider != None:
        phi_slider.on_changed(update_plot)
    if theta_slider != None:
        theta_slider.on_changed(update_plot)
    r_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
        
    return psi_slider, phi_slider, theta_slider, r_slider, data_chooser

#======================================================================================

def calculate_sb():
    # calculate soluble-beta interaction values
    sb_vals = {}
    sb_vals['md'] = np.array([md.sb_interaction(r, z, phi)
                              for phi in thetas for r in rs for z in zs])
    sb_vals['md'] = sb_vals['md'].reshape(theta_points, r_points, z_points)
    md_sb_min = sb_vals['md'].min()
    print "md_sb_vals.min() =", md_sb_min
    sb_vals['mc'] = np.array([mc.sb_interaction(r, z, phi, eps=-md_sb_min*sb_factor)
                              for phi in thetas for r in rs for z in zs])
    sb_vals['mc'] = sb_vals['mc'].reshape(theta_points, r_points, z_points)
    sb_vals['diff'] = sb_vals['md'] - sb_vals['mc']
    
    return sb_vals, md_sb_min
    
def calculate_bb(theta = None, phi = None, psi = None):
    '''
    Calculates beta-beta interaction values...
    theta : polar angle of rod at (r,z); if None calculated for all thetas
    phi : azimuthal angle of rod at (r,z); if None calculated for all phis
    psi : angle of patch of rod as (0,0); if None calculated for all psis
    '''
    bb_vals = {}
    if theta != None:
        if phi != None:
            bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, 0)
                                          for phi in [phi] for theta in [theta] for r in rs for z in zs])
            bb_vals['mc'] = bb_vals['mc'].reshape(1, 1, 1, r_points, z_points)
            if psi != None:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in [psi] for phi in [phi] for theta in [theta] for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(1, 1, 1, r_points, z_points)
            else:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in thetas for phi in [phi] for theta in [theta] for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(theta_points, 1, 1, r_points, z_points)

                bb_vals['mc'] = np.array([bb_vals['mc'][0] for psi in thetas]) # so MC data has the same "shape"
        else:
            bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, 0)
                                          for phi in phis for theta in [theta] for r in rs for z in zs])
            bb_vals['mc'] = bb_vals['mc'].reshape(1, phi_points, 1, r_points, z_points)
            if psi != None:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in [psi] for phi in phis for theta in [theta] for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(1, phi_points, 1, r_points, z_points)
            else:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in thetas for phi in phis for theta in [theta] for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(theta_points, phi_points, 1, r_points, z_points)
                
                bb_vals['mc'] = np.array([bb_vals['mc'][0] for psi in thetas]) # so MC data has the same "shape"
    else:
        if phi != None:
            bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, 0)
                                          for phi in [phi] for theta in thetas for r in rs for z in zs])
            bb_vals['mc'] = bb_vals['mc'].reshape(1, 1, theta_points, r_points, z_points)
            if psi != None:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in [psi] for phi in [phi] for theta in thetas for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(1, 1, theta_points, r_points, z_points)
            else:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in thetas for phi in [phi] for theta in thetas for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(theta_points, 1, theta_points, r_points, z_points)
                
                bb_vals['mc'] = np.array([bb_vals['mc'][0] for psi in thetas]) # so MC data has the same "shape"
        else:
            bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, 0)
                                          for phi in phis for theta in thetas for r in rs for z in zs])
            bb_vals['mc'] = bb_vals['mc'].reshape(1, phi_points, theta_points, r_points, z_points)
            if psi != None:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in [psi] for phi in phis for theta in thetas for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(1, phi_points, theta_points, r_points, z_points)
            else:
                bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi)
                                          for psi in thetas for phi in phis for theta in thetas for r in rs for z in zs])
                bb_vals['md'] = bb_vals['md'].reshape(theta_points, phi_points, theta_points, r_points, z_points)
                
                bb_vals['mc'] = np.array([bb_vals['mc'][0] for psi in thetas]) # so MC data has the same "shape"
    
    md_bb_min = bb_vals['md'].min()
    print "md_bb_vals.min() = ", md_bb_min
    
    bb_vals['mc'] *= -md_bb_min*bb_factor
    
    bb_vals['diff'] = bb_vals['md'] - bb_vals['mc']
    
    return bb_vals, md_bb_min

#=======================================================================================

cfg_filename = '9-4.cfg'
model = Model(os.path.join('./test cases/',cfg_filename))
mc.setup(model)
md.setup(model)

# plot parameters 
rmin = -0*model.rod_radius; zero_r = 0
rmax = 5*model.rod_radius
zmin = -0*model.rod_radius; zero_z = 0
zmax = 7.5*model.rod_radius
point_density = 0.1
z_points = int((zmax-zmin)/point_density) + 1
r_points = int((rmax-rmin)/point_density) + 1
phi_points = 12 + 1 # 30 deg slices
theta_points = 6 + 1 # 30 deg slices

# plotting points
zs = np.linspace(zmin, zmax, z_points)
rs = np.linspace(rmin, rmax, r_points)
phis = np.linspace(0, 2*np.pi, phi_points)
thetas = np.linspace(0, np.pi, theta_points)

def z_index(z):
    return int((z/zs[-1])*(z_points-1))
    
def r_index(r):
    return int((r/rs[-1])*(r_points-1))

def phi_index(phi):
    return int((phi/phis[-1])*(phi_points-1))
    
def theta_index(theta):
    return int((theta/thetas[-1])*(theta_points-1))

# figure parameters
axis_font = {'size':12}
cmap = plt.get_cmap("RdBu") #coolwarm, seismic, viridis
widget_color = 'lightgoldenrodyellow'

sb_factor = 1/2.0
bb_factor = 1/3.0

widgets = []

draw = ['sb', 'bb']

if 'sb' in draw:
    sb_vals, md_sb_min = calculate_sb()
    widgets.append(draw_sb_2D(sb_vals, md_sb_min, prefix=cfg_filename))
    widgets.append(draw_sb_z_slice(sb_vals, md_sb_min, prefix=cfg_filename))
    widgets.append(draw_sb_r_slice(sb_vals, md_sb_min, prefix=cfg_filename))
if 'bb' in draw:
    bb_vals, md_bb_min = calculate_bb(phi=0, theta=0)
    widgets.append(draw_bb_2D(bb_vals, md_bb_min, prefix=cfg_filename))
    widgets.append(draw_bb_z_slice(bb_vals, md_bb_min, prefix=cfg_filename))
    widgets.append(draw_bb_r_slice(bb_vals, md_bb_min, prefix=cfg_filename))

plt.show()
