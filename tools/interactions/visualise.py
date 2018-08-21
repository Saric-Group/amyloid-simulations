# encoding: utf-8
'''
This can draw interactive plots of MC and MD potentials (and their difference),
defined in "mc_interactions" and "md_interactions", with different possible
interaction potentials for MD.

Created on 23 Apr 2018

@author: Eugen Rožić
'''

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

import numpy as np

from tools.interactions import mc
from tools.interactions import md
from potentials import lj_n_m, lj_cos_sq, gauss, morse

# MD interaction parameters
M = 6
r_int = 0.25*md.r_body
#delta_int = 2 - ((md.N - 2)*(2 - md.delta_body)*md.r_body)/((M - 1)*r_int)
delta_int = -2.0*md.r_body
int_range = 1.75*md.r_body
bulge_out = 0.0*md.r_body
md.set_parameters(M, r_int, delta_int, bulge_out, int_range)

# interactions
lj_12_6 = lambda r, sigma, cutoff : lj_n_m(12, 6, r, sigma, cutoff)
lj_12_4 = lambda r, sigma, cutoff : lj_n_m(12, 4, r, sigma, cutoff)
gauss_05 = lambda r, sigma, cutoff : gauss(0.5*md.r_body, r, sigma, cutoff)
gauss_06 = lambda r, sigma, cutoff : gauss(0.6*md.r_body, r, sigma, cutoff)
morse_275 = lambda r, sigma, cutoff : morse(2.75/md.r_body, r, sigma, cutoff)
morse_25 = lambda r, sigma, cutoff : morse(2.5/md.r_body, r, sigma, cutoff)
morse_20 = lambda r, sigma, cutoff : morse(2.0/md.r_body, r, sigma, cutoff)

# plot parameters 
rmin = -0*md.r_body; zero_r = 0
rmax = 5*md.r_body
zmin = -0*md.r_body; zero_z = 0
zmax = 7.5*md.r_body
point_density = 0.2
z_points = int((zmax-zmin)/point_density) + 1
r_points = int((rmax-rmin)/point_density) + 1
phi_points = 16#32
theta_points = 8

# plotting points
zs = np.linspace(zmin, zmax, z_points)
rs = np.linspace(rmin, rmax, r_points)
phis = np.linspace(0, 2*np.pi, phi_points)
thetas = np.linspace(0, np.pi, theta_points)

# figure parameters
axis_font = {'size':12}
cmap = plt.get_cmap("RdBu") #coolwarm, seismic, viridis
widget_color = 'lightgoldenrodyellow'


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


def draw_sb(prefix, int_f, polar=True):
    
    # calculate soluble-beta interaction values
    sb_vals = {}
    sb_vals['md'] = np.array([md.sb_interaction(int_f, r, z, phi)  for phi in thetas for r in rs for z in zs])
    sb_vals['md'] = sb_vals['md'].reshape(theta_points, r_points, z_points)
    md_sb_min = sb_vals['md'].min()
    print "md_sb_vals.min() =", md_sb_min
    sb_vals['mc'] = np.array([mc.sb_interaction(r, z, phi, eps=-md_sb_min/2) for phi in thetas for r in rs for z in zs])
    sb_vals['mc'] = sb_vals['mc'].reshape(theta_points, r_points, z_points)
    sb_vals['diff'] = sb_vals['md'] - sb_vals['mc']
    
    init_theta_index = 0
    init_data = 'md'
    
    # do drawing
    fig, ax = plt.subplots(num=prefix+" soluble-beta interactions", figsize=(11,7))
    fig.subplots_adjust(left=0.13, bottom=0.13, right=0.99, top=0.99) #0.86 x 0.86
    
    img = ax.imshow(sb_vals[init_data][init_theta_index], extent=[zmin, zmax, rmin, rmax], vmin=md_sb_min, vmax=-md_sb_min,
                    origin="lower", interpolation='bilinear', cmap=cmap)
    draw_models(ax)
    
    ax.set_xlabel(r'$z$', **axis_font)
    ax.set_ylabel(r'$r$', **axis_font)
    ax.axis([zmin, zmax, rmin, rmax])
    
    if polar:
        fig_polar, ax_polar = plt.subplots(subplot_kw={'projection':'polar'}, num=fig.canvas.get_window_title() + ' - side view (z=0)')
     
        angle, rad = np.meshgrid(thetas, rs)
        polar_mesh = ax_polar.pcolormesh(angle, rad, sb_vals[init_data].T[zero_z], vmin=md_sb_min, vmax=-md_sb_min,
                                         cmap=cmap)
     
        ax_polar.add_artist(plt.Circle((0,0), 2*md.r_body, transform=plt.gca().transData._b, fill=True, color='#DDDDDD', alpha=0.7))
     
        ax_polar.axis(rmax=rmax)
    
    theta_slider_axes = fig.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'Phi', thetas[0], thetas[-1], valinit=thetas[init_theta_index])
    
    data_chooser_axes = fig.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_theta(curr_theta):
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        curr_data = data_chooser.value_selected
        img.set_data(sb_vals[curr_data][theta_index])
        fig.canvas.draw_idle()
        
    def update_data(curr_data):
        curr_theta = theta_slider.val
        theta_index = int((curr_theta/thetas[-1])*(theta_points-1))
        img.set_data(sb_vals[curr_data][theta_index])
        if (polar):
            polar_mesh.set_array(sb_vals[curr_data].T[zero_z][:-1, :-1].ravel())
        fig.canvas.draw_idle()
        fig_polar.canvas.draw_idle()
    
    theta_slider.on_changed(update_theta)    
    data_chooser.on_clicked(update_data)
        
    return theta_slider, data_chooser


def draw_bb(prefix, inf_f):
    
    # calculate beta-beta interaction values
    bb_vals = {}
    bb_vals['md'] = np.array([md.bb_interaction(inf_f, r, z, theta, phi) for phi in phis for theta in thetas for r in rs for z in zs])
    bb_vals['md'] = bb_vals['md'].reshape(phi_points, theta_points, r_points, z_points)
    md_bb_min = bb_vals['md'].min()
    print "md_bb_vals.min() = ", md_bb_min
    bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi, eps=-md_bb_min/2) for phi in phis for theta in thetas for r in rs for z in zs])
    bb_vals['mc'] = bb_vals['mc'].reshape(phi_points, theta_points, r_points, z_points)
    bb_vals['diff'] = bb_vals['md'] - bb_vals['mc']
    
    init_theta_index = 0
    init_phi_index = 0
    init_data = 'md'
    
    # do drawing
    fig, ax = plt.subplots(num=prefix+" beta-beta interactions", figsize=(11,7))
    fig.subplots_adjust(left=0.131, bottom=0.19, right=0.99, top=0.99) #0.88 x 0.88
    
    img = ax.imshow(bb_vals[init_data][init_phi_index][init_theta_index], extent=[zmin, zmax, rmin, rmax],
                    vmin=md_bb_min, vmax=-md_bb_min, origin="lower", interpolation='bilinear', cmap=cmap)
    draw_models(ax)
    
    ax.set_xlabel(r'$z$', **axis_font)
    ax.set_ylabel(r'$r$', **axis_font)
    ax.axis([zmin, zmax, rmin, rmax])
    
    theta_slider_axes = fig.add_axes([0.41, 0.07, 0.35, 0.03], facecolor=widget_color)
    theta_slider = Slider(theta_slider_axes, 'Theta', thetas[0], thetas[-1], valinit=thetas[init_theta_index])
    
    phi_slider_axes = fig.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, 'Phi', phis[0], phis[-1], valinit=phis[init_phi_index])
    
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

# sb_widgets = draw_sb("LJ_12-6", lj_12_6, polar=True)
# sb_widgets = draw_sb("LJ_12-4", lj_12_4, polar=True)
# sb_widgets = draw_sb("morse_2.5", morse_25, polar=True)
# sb_widgets = draw_sb("morse_2.0", morse_20, polar=True)
# sb_widgets = draw_sb("gauss_0.5", gauss_05, polar=True)
# sb_widgets = draw_sb("gauss_0.6", gauss_06, polar=True)
sb_widgets = draw_sb("lj/cos_sq", lj_cos_sq, polar=True)

# bb_widgets = draw_bb("LJ_12-6", lj_12_6)
# bb_widgets = draw_bb("LJ_12-4", lj_12_4)
# bb_widgets = draw_bb("morse_2.5", morse_25)
# bb_widgets = draw_bb("morse_2.0", morse_20)
# bb_widgets = draw_bb("gauss_0.5", gauss_05)
# bb_widgets = draw_bb("gauss_0.6", gauss_06)
bb_widgets = draw_bb("LJ/cos_sq", lj_cos_sq)

plt.show()
