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

def draw_models(axes = None, r = 0, z = 0):
    if axes is None :
        axes = plt.gca()
    # the "halo" of the interacting body (closest possible approach)
    for i in range(md.N):
        z_i = z + (i - (md.N-1)/2.)*(2 - md.delta_body)*md.r_body 
        axes.add_patch(plt.Circle((z_i,r), 2*md.r_body, fill=True, color='#DDDDDD', alpha=0.7))
    # the soft-sphere interaction shell...
    for i in range(md.N):
        z_i = z + (i - (md.N-1)/2.)*(2 - md.delta_body)*md.r_body 
        axes.add_patch(plt.Circle((z_i,r), md.r_body, fill=True, color='#777777'))
    # the interaction centers
    r = r + (md.r_body - md.r_int + md.bulge_out)
    for i in range(md.M):
        z_i = z + (i - (md.M-1)/2.)*(2 - md.delta_int)*md.r_int 
        axes.add_patch(plt.Circle((z_i,r), md.r_int, fill=True, color='red'))
    # the patch axis
    axes.plot([-mc.L_patch/2, mc.L_patch/2], [0, 0], color='red', linestyle='-', linewidth=5)

def draw_sb_2D(sb_vals, md_sb_min, prefix = ''):
    
    fig_2D, ax_2D = plt.subplots(num=prefix+" soluble-beta interactions (side 2D)", figsize=(11,7))
    fig_2D.subplots_adjust(left=0.13, bottom=0.13, right=0.99, top=0.99) #0.86 x 0.86
    
    img = ax_2D.imshow(sb_vals['md'][0], extent=[zmin, zmax, rmin, rmax],
                       vmin=md_sb_min, vmax=-md_sb_min,
                       origin="lower", interpolation='bilinear', cmap=cmap)
    draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    theta_slider_axes = fig_2D.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    data_chooser_axes = fig_2D.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_theta(curr_theta):
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        curr_data = data_chooser.value_selected
        img.set_data(sb_vals[curr_data][theta_index])
        fig_2D.canvas.draw_idle()
        
    def update_data(curr_data):
        curr_theta = theta_slider.val
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        img.set_data(sb_vals[curr_data][theta_index])
        fig_2D.canvas.draw_idle()
    
    theta_slider.on_changed(update_theta)    
    data_chooser.on_clicked(update_data)
    
    return theta_slider, data_chooser

def draw_sb_z_slice(sb_vals, md_sb_min, prefix = ''):
    
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=prefix+" soluble-beta interactions (z-slice)", figsize=(9,5))
    fig_r.subplots_adjust(left=0.18, bottom=0.2, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = ax_r.plot(rs, sb_vals['md'][0].T[zero_z], 'r-', lw=1.0)
    ax_r.axvline(2.0*md.r_body, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*md.r_body, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.axis([rmin, rmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    theta_slider_axes = fig_r.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    z_slider_axes = fig_r.add_axes([0.41, 0.06, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, 'z', zs[0], zs[-1], valinit=zs[zero_z])
    
    data_chooser_axes = fig_r.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_z(curr_z):
        curr_z_index = z_index(curr_z)
        curr_theta_index = theta_index(theta_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    def update_theta(curr_theta):
        curr_z_index = z_index(z_slider.val)
        curr_theta_index = theta_index(curr_theta)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index].T[curr_z_index])
        fig_r.canvas.draw_idle()
        
    def update_data(curr_data):
        curr_z_index = z_index(z_slider.val)
        curr_theta_index = theta_index(theta_slider.val)
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    z_slider.on_changed(update_z)
    theta_slider.on_changed(update_theta)
    data_chooser.on_clicked(update_data)
    
    return z_slider, theta_slider, data_chooser

def draw_sb_r_slice(sb_vals, md_sb_min, prefix = ''):
    
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=prefix+" soluble-beta interactions (r-slice)", figsize=(11,5))
    fig_z.subplots_adjust(left=0.18, bottom=0.2, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = ax_z.plot(zs, sb_vals['md'][0][zero_r], 'r-', lw=1.0)
    ax_z.axvline(3.0*md.r_body, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*md.r_body, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.axis([zmin, zmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    theta_slider_axes = fig_z.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'phi', thetas[0], thetas[-1], valinit=thetas[0])
    
    r_slider_axes = fig_z.add_axes([0.41, 0.06, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, 'r', rs[0], rs[-1], valinit=rs[zero_r])
    
    data_chooser_axes = fig_z.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_r(curr_r):
        curr_r_index = r_index(curr_r)
        curr_theta_index = theta_index(theta_slider.val)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index][curr_r_index])
        fig_z.canvas.draw_idle()
    
    def update_theta(curr_theta):
        curr_r_index = r_index(r_slider.val)
        curr_theta_index = theta_index(curr_theta)
        curr_data = data_chooser.value_selected
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index][curr_r_index])
        fig_z.canvas.draw_idle()
        
    def update_data(curr_data):
        curr_r_index = r_index(r_slider.val)
        curr_theta_index = theta_index(theta_slider.val)
        lines[0].set_ydata(sb_vals[curr_data][curr_theta_index][curr_r_index])
        fig_z.canvas.draw_idle()
    
    r_slider.on_changed(update_r)
    theta_slider.on_changed(update_theta)
    data_chooser.on_clicked(update_data)
    
    return r_slider, theta_slider, data_chooser

def draw_bb(bb_vals, md_bb_min, prefix = ''):
    
    fig, ax = plt.subplots(num=prefix+" beta-beta interactions", figsize=(11,7))
    fig.subplots_adjust(left=0.131, bottom=0.19, right=0.99, top=0.99) #0.88 x 0.88
    
    img = ax.imshow(bb_vals['md'][0][0], extent=[zmin, zmax, rmin, rmax],
                    vmin=md_bb_min, vmax=-md_bb_min, origin="lower", interpolation='bilinear', cmap=cmap)
    draw_models(ax)
    
    ax.set_xlabel(r'$z$', **axis_font)
    ax.set_ylabel(r'$r$', **axis_font)
    ax.axis([zmin, zmax, rmin, rmax])
    
    theta_slider_axes = fig.add_axes([0.41, 0.07, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'Theta', thetas[0], thetas[-1], valinit=thetas[0])
    
    phi_slider_axes = fig.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, 'Phi', phis[0], phis[-1], valinit=phis[0])
    
    data_chooser_axes = fig.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_theta(curr_theta):
        curr_phi = phi_slider.val
        curr_data = data_chooser.value_selected
        
        phi_index = int((curr_phi/phis[-1])*(phi_points-1))
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        
        img.set_data(bb_vals[curr_data][phi_index][theta_index])
        fig.canvas.draw_idle()
        
    def update_phi(curr_phi):
        curr_theta = theta_slider.val
        curr_data = data_chooser.value_selected
        
        phi_index = int((curr_phi/phis[-1])*(phi_points-1))
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        
        img.set_data(bb_vals[curr_data][phi_index][theta_index])
        fig.canvas.draw_idle()
        
    def update_data(curr_data):
        curr_phi = phi_slider.val
        curr_theta = theta_slider.val
        
        phi_index = int((curr_phi/phis[-1])*(phi_points-1))
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        
        img.set_data(bb_vals[curr_data][phi_index][theta_index])
        fig.canvas.draw_idle()
    
    theta_slider.on_changed(update_theta)
    phi_slider.on_changed(update_phi)
    data_chooser.on_clicked(update_data)
        
    return theta_slider, phi_slider, data_chooser

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
    
def calculate_bb():
    # calculate beta-beta interaction values
    bb_vals = {}
    bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi) for phi in phis for theta in thetas for r in rs for z in zs])
    bb_vals['md'] = bb_vals['md'].reshape(phi_points, theta_points, r_points, z_points)
    md_bb_min = bb_vals['md'].min()
    print "md_bb_vals.min() = ", md_bb_min
    bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, eps=-md_bb_min*bb_factor)
                              for phi in phis for theta in thetas for r in rs for z in zs])
    bb_vals['mc'] = bb_vals['mc'].reshape(phi_points, theta_points, r_points, z_points)
    bb_vals['diff'] = bb_vals['md'] - bb_vals['mc']
    return bb_vals, md_bb_min

#=======================================================================================

cfg_filename = '6-3_cos-sq_0.25.cfg'
model = Model(os.path.join('./test cases/',cfg_filename))
mc.setup(model)
md.setup(model)

# plot parameters 
rmin = -0*md.r_body; zero_r = 0
rmax = 5*md.r_body
zmin = -0*md.r_body; zero_z = 0
zmax = 7.5*md.r_body
point_density = 0.1
z_points = int((zmax-zmin)/point_density) + 1
r_points = int((rmax-rmin)/point_density) + 1
phi_points = 64
theta_points = 32

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

#calculating & drawing soluble-beta
sb_vals, md_sb_min = calculate_sb()
sb_widgets = []
sb_widgets.append(draw_sb_2D(sb_vals, md_sb_min, cfg_filename))
sb_widgets.append(draw_sb_z_slice(sb_vals, md_sb_min, cfg_filename))
sb_widgets.append(draw_sb_r_slice(sb_vals, md_sb_min, cfg_filename))

#calculating & drawing beta-beta
# bb_vals, md_bb_min = calculate_bb()
# bb_widgets = draw_bb(bb_vals, md_bb_min, cfg_filename)

plt.show()
