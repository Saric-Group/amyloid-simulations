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
    fig_2D.subplots_adjust(left=0.13, bottom=0.15, right=0.99, top=0.99) #0.86 x 0.86
    
    img = ax_2D.imshow(sb_vals['md'][0], extent=[zmin, zmax, rmin, rmax],
                       vmin=md_sb_min, vmax=-md_sb_min,
                       origin="lower", interpolation='bilinear', cmap=img_cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    phi_slider_axes = fig_2D.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, r'$\phi$', 0, theta_points-1, valinit=zero_theta,
                        valstep=1, valfmt = '%1.1f deg')
    
    data_chooser_axes = fig_2D.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_img(_):
        curr_phi_index = int(phi_slider.val)
        curr_phi_val = thetas[curr_phi_index]
        phi_slider.valtext.set_text(phi_slider.valfmt % (curr_phi_val*180/np.pi))
        curr_data = data_chooser.value_selected
        draw_models(ax_2D, phi = curr_phi_val, old = model_patches)
        img.set_data(sb_vals[curr_data][curr_phi_index])
        fig_2D.canvas.draw_idle()
    
    phi_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
    
    return phi_slider, data_chooser

def draw_sb_z_slice(sb_vals, md_sb_min, prefix = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=prefix+" soluble-beta interaction (z-slice)", figsize=(9,5))
    fig_r.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = []
    for i in range(theta_points):
        lines.append(ax_r.plot(rs, sb_vals['md'][i].T[zero_z], 'r-', lw=1.0,
                               label=r'$\phi = {:1.1f}^\circ$'.format(thetas[i]*180/np.pi),
                               color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.legend(loc = 'lower left', prop = axis_font)
    ax_r.axis([rmin, rmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    z_slider_axes = fig_r.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, r'$z$', 0, z_points-1, valinit=zero_z,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_r.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_plot(_):
        curr_z_index = int(z_slider.val)
        curr_z_val = zs[curr_z_index]
        z_slider.valtext.set_text(z_slider.valfmt % curr_z_val)
        curr_data = data_chooser.value_selected
        
        for i in range(theta_points):
            lines[i][0].set_ydata(sb_vals[curr_data][i].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    z_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    
    return z_slider, data_chooser

def draw_sb_r_slice(sb_vals, md_sb_min, prefix = ''):
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=prefix+" soluble-beta interaction (r-slice)", figsize=(11,5))
    fig_z.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.97) #0.86 x 0.86
    
    lines = []
    for i in range(theta_points):
        lines.append(ax_z.plot(zs, sb_vals['md'][i][zero_r], 'r-', lw=1.0,
                               label=r'$\phi = {:1.1f}^\circ$'.format(thetas[i]*180/np.pi),
                               color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.legend(loc = 'lower right', prop = axis_font)
    ax_z.axis([zmin, zmax, 1.05*md_sb_min, -1.05*md_sb_min])
    
    r_slider_axes = fig_z.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, r'$r$', 0, r_points-1, valinit=zero_r,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_z.add_axes([0.01, 0.45, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, sb_vals.keys())
    
    def update_plot(_):
        curr_r_index = int(r_slider.val)
        curr_r_val = rs[curr_r_index]
        r_slider.valtext.set_text(r_slider.valfmt % curr_r_val)
        curr_data = data_chooser.value_selected
        for i in range(theta_points):
            lines[i][0].set_ydata(sb_vals[curr_data][i][curr_r_index])
        fig_z.canvas.draw_idle()
    
    r_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    
    return r_slider, data_chooser

#===================================================================

def draw_bb_2D(bb_vals, md_bb_min, prefix = ''):
    
    fig_2D, ax_2D = plt.subplots(num=prefix+" beta-beta interaction", figsize=(10,7))
    fig_2D.subplots_adjust(left=0.14, bottom=0.20, right=0.99, top=0.99) #0.88 x 0.88
    
    img = ax_2D.imshow(bb_vals['md'][zero_phi][zero_theta],
                    extent=[zmin, zmax, rmin, rmax], vmin=md_bb_min, vmax=-md_bb_min,
                    origin="lower", interpolation='bilinear', cmap=img_cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    psi1_slider = psi2_slider = None
    
    if len(bb_vals['md']) > 1:
        psi2_slider_axes = fig_2D.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
        
    if len(bb_vals['md'][0]) > 1:
        psi1_slider_axes = fig_2D.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
        psi1_slider = Slider(psi1_slider_axes, r'$\psi_1$', 0, theta_points-1, valinit=zero_theta,
                        valstep=1, valfmt = '%1.1f deg')
    
    data_chooser_axes = fig_2D.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_img(_):
        curr_psi1_index = curr_psi2_index = 0
        if psi1_slider != None:
            curr_psi1_index = int(psi1_slider.val)
            curr_psi1_val = thetas[curr_psi1_index]
            psi1_slider.valtext.set_text(psi1_slider.valfmt % (curr_psi1_val*180/np.pi))
            draw_models(ax_2D, phi = curr_psi1_val, old = model_patches)
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (curr_psi2_val*180/np.pi))
        curr_data = data_chooser.value_selected
        
        img.set_data(bb_vals[curr_data][curr_psi2_index][curr_psi1_index])
        fig_2D.canvas.draw_idle()
    
    if psi1_slider != None:
        psi1_slider.on_changed(update_img)
    if psi2_slider != None:
        psi2_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
        
    return psi1_slider, psi2_slider, data_chooser

def draw_bb_z_slice(bb_vals, md_bb_min, prefix = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=prefix+" beta-beta interaction (z-slice)", figsize=(9,6))
    fig_r.subplots_adjust(left=0.14, bottom=0.20, right=0.99, top=0.99) #0.86 x 0.86
    
    lines = []
    for i in range(theta_points):
        lines.append(ax_r.plot(rs, bb_vals['md'][zero_phi][i].T[zero_z], 'r-', lw=1.0,
                               label=r'$\psi_1 = {:1.1f}^\circ$'.format(thetas[i]*180/np.pi),
                                   color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.legend(loc = 'lower left', prop = axis_font)
    ax_r.axis([rmin, rmax, 1.05*md_bb_min, -1.05*md_bb_min])
    
    psi2_slider = None
    
    if len(bb_vals['md']) > 1:
        psi2_slider_axes = fig_r.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
    
    z_slider_axes = fig_r.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, r'$z$', 0, z_points-1, valinit=zero_z,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_r.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_plot(_):
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (curr_psi2_val*180/np.pi))
        curr_z_index = int(z_slider.val)
        curr_z_val = zs[curr_z_index]
        z_slider.valtext.set_text(z_slider.valfmt % curr_z_val)
        curr_data = data_chooser.value_selected
        
        for i in range(theta_points):
            lines[i][0].set_ydata(bb_vals[curr_data][curr_psi2_index][i].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    if psi2_slider != None:
        psi2_slider.on_changed(update_plot)
    z_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
        
    return psi2_slider, z_slider, data_chooser

def draw_bb_r_slice(bb_vals, md_bb_min, prefix = ''):
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=prefix+" beta-beta interaction (r-slice)", figsize=(9,6))
    fig_z.subplots_adjust(left=0.14, bottom=0.20, right=0.99, top=0.99) #0.86 x 0.86
    
    lines = []
    for i in range(theta_points):
        lines.append(ax_z.plot(zs, bb_vals['md'][zero_phi][i][zero_r], 'r-', lw=1.0,
                               label=r'$\psi_1 = {:1.1f}^\circ$'.format(thetas[i]*180/np.pi),
                                   color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.legend(loc = 'lower right', prop = axis_font)
    ax_z.axis([zmin, zmax, 1.05*md_bb_min, -1.05*md_bb_min])
    
    psi2_slider = None
    
    if len(bb_vals['md']) > 1:
        psi2_slider_axes = fig_z.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
    
    r_slider_axes = fig_z.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, r'$r$', 0, r_points-1, valinit=zero_r,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_z.add_axes([0.01, 0.45, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, bb_vals.keys())
    
    def update_plot(_):
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (curr_psi2_val*180/np.pi))
        curr_r_index = int(r_slider.val)
        curr_r_val = zs[curr_r_index]
        r_slider.valtext.set_text(r_slider.valfmt % curr_r_val)
        curr_data = data_chooser.value_selected
        
        for i in range(theta_points):
            lines[i][0].set_ydata(bb_vals[curr_data][curr_psi2_index][i][curr_r_index])
        fig_z.canvas.draw_idle()
    
    if psi2_slider != None:
        psi2_slider.on_changed(update_plot)
    r_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
        
    return psi2_slider, r_slider, data_chooser

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
    
def calculate_bb(theta = 0, phi = 0):
    '''
    Calculates beta-beta interaction values...
    theta : polar angle of rod at (r,z)
    phi : azimuthal angle of rod at (r,z)
    psi1 : angle of patch of rod at (0,0); if "None" calculated for psi1 in [0,2pi>
    psi1 : angle of patch of rod at (r,z); if "None" calculated for psi2 in [0,2pi>
    '''
    bb_vals = {}
    
    bb_vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi) for r in rs for z in zs])
    bb_vals['mc'] = bb_vals['mc'].reshape(r_points, z_points)
    bb_vals['mc'] = np.array([[bb_vals['mc'] for psi1 in thetas] for psi2 in phis]) # so MC data has the same "shape"
    
    bb_vals['md'] = np.array([md.bb_interaction(r, z, theta, phi, psi1, psi2)
                              for psi2 in phis for psi1 in thetas for r in rs for z in zs])
    bb_vals['md'] = bb_vals['md'].reshape(phi_points, theta_points, r_points, z_points)
    
    md_bb_min = bb_vals['md'].min()
    print "md_bb_vals.min() = ", md_bb_min
    
    bb_vals['mc'] *= -md_bb_min*bb_factor
    
    bb_vals['diff'] = bb_vals['md'] - bb_vals['mc']
    
    return bb_vals, md_bb_min

#=======================================================================================

cfg_filename = '7-545.cfg'
model = Model(os.path.join('./test cases/',cfg_filename))
mc.setup(model)
md.setup(model)

# plot parameters 
rmin = -0*model.rod_radius; rmax = 5.0*model.rod_radius
zmin = -0*model.rod_radius; zmax = 7.5*model.rod_radius
point_density = 0.1
z_points = int((zmax-zmin)/point_density) + 1; zero_z = 0
r_points = int((rmax-rmin)/point_density) + 1; zero_r = 0
phimin = -np.pi; phimax = np.pi
phi_points = 8; zero_phi = 4 # 45 deg slices
thetamin = 0; thetamax = np.pi
theta_points = 4 + 1; zero_theta = 0 # 45 deg slices

# plotting points
zs = np.linspace(zmin, zmax, z_points)
rs = np.linspace(rmin, rmax, r_points)
phis = np.linspace(phimin, phimax, phi_points, endpoint=False)
thetas = np.linspace(thetamin, thetamax, theta_points)

# figure parameters
axis_font = {'size':13}
img_cmap = plt.get_cmap("RdBu")
plot_cmap = plt.cm.get_cmap('nipy_spectral')
widget_color = 'lightgoldenrodyellow'

sb_factor = 1/2.0
bb_factor = 1/3.0

widgets = []

draw = ['sb']

if 'sb' in draw:
    sb_vals, md_sb_min = calculate_sb()
    widgets.append(draw_sb_2D(sb_vals, md_sb_min, prefix=cfg_filename))
    widgets.append(draw_sb_z_slice(sb_vals, md_sb_min, prefix=cfg_filename))
    widgets.append(draw_sb_r_slice(sb_vals, md_sb_min, prefix=cfg_filename))
if 'bb' in draw:
    bb_vals, md_bb_min = calculate_bb(theta = 0, phi = 0)
    widgets.append(draw_bb_2D(bb_vals, md_bb_min, prefix=cfg_filename))
    widgets.append(draw_bb_z_slice(bb_vals, md_bb_min, prefix=cfg_filename))
    widgets.append(draw_bb_r_slice(bb_vals, md_bb_min, prefix=cfg_filename))

plt.show()
