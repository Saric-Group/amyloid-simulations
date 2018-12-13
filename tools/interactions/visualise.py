# encoding: utf-8
'''
An interactive console application for drawing interactive plots of interactions
between rods from a lammps_multistate_rods model, defined in a config file, and
the corresponding interactions in the MC model by Anđela Šarić along with 
an interactive comparison between the two.

Created on 23 Apr 2018

@author: Eugen Rožić
'''

import numpy as np

from lammps_multistate_rods import Model
import lammps_multistate_rods.tools.interactions as md
import mc

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, AxesWidget
import six

class VertSlider(AxesWidget):
    """
    A slider representing a floating point range.

    For the slider to remain responsive you must maintain a
    reference to it.

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*.

        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...).

        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in

        label : str
            Slider label

        valmin : float
            The minimum value of the slider

        valmax : float
            The maximum value of the slider

        valinit : float
            The slider initial position

        label : str
            The slider label

        valfmt : str
            Used to format the slider value, fprint format string

        dragging : bool
            if the slider can be dragged by the mouse

        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.hline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}
        
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            val = self.valmin
        elif val >= self.valmax:
            val = self.valmax

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for _, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)
            
#========================================================================================

def draw_models(axes = None, r = 0, z = 0, phi = 0, old = []):
    
    if axes is None :
        axes = plt.gca()
    r_body = md.model.rod_radius
    r_int = md.model.patch_bead_radii
    
    for i in range(len(old),0,-1):
        old[i-1].remove()
        del old[i-1]
            
    new = old
    
    # the MD model body
    for i in range(md.M[0]): 
        new.append(
            axes.add_patch(plt.Circle((z + md.bead_z[0][i], r), 2*r_body, fill=True, color='#DDDDDD', alpha=0.7)))
    for i in range(md.M[0]):
        new.append(
            axes.add_patch(plt.Circle((z + md.bead_z[0][i], r), r_body, fill=True, color='#777777')))
        
    # the MC model patch axis
    new.append(
        axes.add_patch(plt.Rectangle((z - mc.L_patch_half, r - 0.02), mc.L_patch, 0.04, color='black')))
    
    # the MD model interaction centers (patches)
    for k in range(1, md.K):
        r_k = r + md.patch_r[k]*np.cos(md.patch_phi[k] + phi)
        for i in range(md.M[k]):
            new.append(
                axes.add_patch(plt.Circle((z + md.bead_z[k][i], r_k), r_int[k-1], fill=True, color='red')))
    
    return new

def draw_point_rod_2D(vals, vals_min, fig_title = ''):
    
    fig_2D, ax_2D = plt.subplots(num=fig_title, figsize=(11,7))
    fig_2D.subplots_adjust(left=0.13, bottom=0.15, right=0.99, top=0.99) #0.86 x 0.86
    
    img = ax_2D.imshow(vals['md'][0], extent=[zmin, zmax, rmin, rmax],
                       vmin=vals_min, vmax=-vals_min,
                       origin="lower", interpolation='bilinear', cmap=img_cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    phi_slider_axes = fig_2D.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    phi_slider = Slider(phi_slider_axes, r'$\phi$', 0, theta_points-1, valinit=zero_theta,
                        valstep=1, valfmt = '%1.1f deg')
    
    data_chooser_axes = fig_2D.add_axes([0.03, 0.55, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_2D.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/2)
    
    def update_img(_):
        curr_phi_index = int(phi_slider.val)
        curr_phi_val = thetas[curr_phi_index]
        phi_slider.valtext.set_text(phi_slider.valfmt % (np.rad2deg(curr_phi_val)))
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(2.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        draw_models(ax_2D, phi = curr_phi_val, old = model_patches)
        img.set_data(to_plot[curr_phi_index])
        fig_2D.canvas.draw_idle()
    
    phi_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
    mc_scale_slider.on_changed(update_img)
    
    return phi_slider, data_chooser, mc_scale_slider

def draw_point_rod_z_slice(vals, vals_min, fig_title = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=fig_title, figsize=(9,5))
    fig_r.subplots_adjust(left=0.18, bottom=0.15, right=0.99, top=0.99) #0.86 x 0.86
    
    lines = []
    for i in range(theta_points):
        lines.append(ax_r.plot(rs, vals['md'][i].T[zero_z], 'r-', lw=1.0,
                               label=r'$\phi = {:1.1f}^\circ$'.format(np.rad2deg(thetas[i])),
                               color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.legend(loc = 'lower left', prop = axis_font)
    ax_r.axis([rmin, rmax, 1.05*vals_min, -1.05*vals_min])
    
    z_slider_axes = fig_r.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, r'$z$', 0, z_points-1, valinit=zero_z,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_r.add_axes([0.02, 0.55, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_r.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/2)
    
    def update_plot(_):
        curr_z_index = int(z_slider.val)
        curr_z_val = zs[curr_z_index]
        z_slider.valtext.set_text(z_slider.valfmt % curr_z_val)
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(2.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        for i in range(theta_points):
            lines[i][0].set_ydata(to_plot[i].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    z_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    mc_scale_slider.on_changed(update_plot)
    
    return z_slider, data_chooser, mc_scale_slider

def draw_point_rod_r_slice(vals, vals_min, fig_title = ''):
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=fig_title, figsize=(11,5))
    fig_z.subplots_adjust(left=0.18, bottom=0.15, right=0.99, top=0.99) #0.86 x 0.86
    
    r_init = int(model.rod_radius*2/dx)
    lines = []
    for i in range(theta_points):
        lines.append(ax_z.plot(zs, vals['md'][i][r_init], 'r-', lw=1.0,
                               label=r'$\phi = {:1.1f}^\circ$'.format(np.rad2deg(thetas[i])),
                               color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.legend(loc = 'lower right', prop = axis_font)
    ax_z.axis([zmin, zmax, 1.05*vals_min, -1.05*vals_min])
    
    r_slider_axes = fig_z.add_axes([0.41, 0.02, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, r'$r$', 0, r_points-1, valinit=r_init,
                      valstep=1, valfmt = '%1.2f')
    r_slider.valtext.set_text(r_slider.valfmt % (r_init*dx))
    
    data_chooser_axes = fig_z.add_axes([0.02, 0.55, 0.07, 0.15], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_z.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/2)
    
    def update_plot(_):
        curr_r_index = int(r_slider.val)
        curr_r_val = rs[curr_r_index]
        r_slider.valtext.set_text(r_slider.valfmt % curr_r_val)
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(2.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        for i in range(theta_points):
            lines[i][0].set_ydata(to_plot[i][curr_r_index])
        fig_z.canvas.draw_idle()
    
    r_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    mc_scale_slider.on_changed(update_plot)
    
    return r_slider, data_chooser, mc_scale_slider

#===================================================================

def draw_rod_rod_2D(vals, vals_min, fig_title = ''):
    
    fig_2D, ax_2D = plt.subplots(num=fig_title, figsize=(10,7))
    fig_2D.subplots_adjust(left=0.14, bottom=0.20, right=0.99, top=0.99) #0.88 x 0.88
    
    if len(vals['md']) > 1:
        psi2_init = zero_phi
    else:
        psi2_init = 0
    img = ax_2D.imshow(vals['md'][psi2_init][zero_theta],
                    extent=[zmin, zmax, rmin, rmax], vmin=vals_min, vmax=-vals_min,
                    origin="lower", interpolation='bilinear', cmap=img_cmap)
    model_patches = draw_models(ax_2D)
    
    ax_2D.set_xlabel(r'$z$', **axis_font)
    ax_2D.set_ylabel(r'$r$', **axis_font)
    ax_2D.axis([zmin, zmax, rmin, rmax])
    
    psi1_slider = psi2_slider = None
    
    if len(vals['md']) > 1:
        psi2_slider_axes = fig_2D.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
        
    if len(vals['md'][0]) > 1:
        psi1_slider_axes = fig_2D.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
        psi1_slider = Slider(psi1_slider_axes, r'$\psi_1$', 0, theta_points-1, valinit=zero_theta,
                        valstep=1, valfmt = '%1.1f deg')
    
    data_chooser_axes = fig_2D.add_axes([0.03, 0.55, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_2D.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/3)
    
    def update_img(_):
        curr_psi1_index = curr_psi2_index = 0
        if psi1_slider != None:
            curr_psi1_index = int(psi1_slider.val)
            curr_psi1_val = thetas[curr_psi1_index]
            psi1_slider.valtext.set_text(psi1_slider.valfmt % (np.rad2deg(curr_psi1_val)))
            draw_models(ax_2D, phi = curr_psi1_val, old = model_patches)
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (np.rad2deg(curr_psi2_val)))
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(3.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        img.set_data(to_plot[curr_psi2_index][curr_psi1_index])
        fig_2D.canvas.draw_idle()
    
    if psi1_slider != None:
        psi1_slider.on_changed(update_img)
    if psi2_slider != None:
        psi2_slider.on_changed(update_img)
    data_chooser.on_clicked(update_img)
    mc_scale_slider.on_changed(update_img)
        
    return psi1_slider, psi2_slider, data_chooser, mc_scale_slider

def draw_rod_rod_z_slice(vals, vals_min, fig_title = ''):
    # 1D r-E plot
    fig_r, ax_r = plt.subplots(num=fig_title, figsize=(9,6))
    fig_r.subplots_adjust(left=0.18, bottom=0.20, right=0.99, top=0.99) #0.86 x 0.86
    
    if len(vals['md']) > 1:
        psi2_init = zero_phi
    else:
        psi2_init = 0
    lines = []
    for i in range(theta_points):
        lines.append(ax_r.plot(rs, vals['md'][psi2_init][i].T[zero_z], 'r-', lw=1.0,
                               label=r'$\psi_1 = {:1.1f}^\circ$'.format(np.rad2deg(thetas[i])),
                                   color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_r.axvline(2.0*model.rod_radius, color='black', linestyle='-', lw=1.0)
    ax_r.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_r.grid()
    
    ax_r.set_xlabel(r'$r$', **axis_font)
    ax_r.set_ylabel(r'$E$', **axis_font)
    ax_r.legend(loc = 'lower left', prop = axis_font)
    ax_r.axis([rmin, rmax, 1.05*vals_min, -1.05*vals_min])
    
    psi2_slider = None
    
    if len(vals['md']) > 1:
        psi2_slider_axes = fig_r.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
    
    z_slider_axes = fig_r.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
    z_slider = Slider(z_slider_axes, r'$z$', 0, z_points-1, valinit=zero_z,
                      valstep=1, valfmt = '%1.2f')
    
    data_chooser_axes = fig_r.add_axes([0.03, 0.55, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_r.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/3)
    
    def update_plot(_):
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (np.rad2deg(curr_psi2_val)))
        else:
            curr_psi2_index = 0
        curr_z_index = int(z_slider.val)
        curr_z_val = zs[curr_z_index]
        z_slider.valtext.set_text(z_slider.valfmt % curr_z_val)
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(3.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        for i in range(theta_points):
            lines[i][0].set_ydata(to_plot[curr_psi2_index][i].T[curr_z_index])
        fig_r.canvas.draw_idle()
    
    if psi2_slider != None:
        psi2_slider.on_changed(update_plot)
    z_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    mc_scale_slider.on_changed(update_plot)
        
    return psi2_slider, z_slider, data_chooser, mc_scale_slider

def draw_rod_rod_r_slice(vals, vals_min, fig_title = ''):
    # 1D z-E plot
    fig_z, ax_z = plt.subplots(num=fig_title, figsize=(9,6))
    fig_z.subplots_adjust(left=0.18, bottom=0.20, right=0.99, top=0.99) #0.86 x 0.86
    
    if len(vals['md']) > 1:
        psi2_init = zero_phi
    else:
        psi2_init = 0
    r_init = int(model.rod_radius*2/dx)
    lines = []
    for i in range(theta_points):
        lines.append(ax_z.plot(zs, vals['md'][psi2_init][i][r_init], 'r-', lw=1.0,
                               label=r'$\psi_1 = {:1.1f}^\circ$'.format(np.rad2deg(thetas[i])),
                                   color = plot_cmap((i+1.0)/(theta_points+1))))
    ax_z.axvline(3.0*model.rod_radius, color='black', linestyle='--', lw=1.0)
    ax_z.axvline(4.0*model.rod_radius, color='black', linestyle='--', lw=0.5)
    ax_z.grid()
    
    ax_z.set_xlabel(r'$z$', **axis_font)
    ax_z.set_ylabel(r'$E$', **axis_font)
    ax_z.legend(loc = 'lower right', prop = axis_font)
    ax_z.axis([zmin, zmax, 1.05*vals_min, -1.05*vals_min])
    
    psi2_slider = None
    
    if len(vals['md']) > 1:
        psi2_slider_axes = fig_z.add_axes([0.38, 0.06, 0.35, 0.03], facecolor=widget_color)
        psi2_slider = Slider(psi2_slider_axes, r'$\psi_2$', 0, phi_points-1, valinit=zero_phi,
                        valstep=1, valfmt = '%1.1f deg')
        psi2_slider.valtext.set_text(psi2_slider.valfmt % 0.0)
    
    r_slider_axes = fig_z.add_axes([0.38, 0.02, 0.35, 0.03], facecolor=widget_color)
    r_slider = Slider(r_slider_axes, r'$r$', 0, r_points-1, valinit=r_init,
                      valstep=1, valfmt = '%1.2f')
    r_slider.valtext.set_text(r_slider.valfmt % (r_init*dx))
    
    data_chooser_axes = fig_z.add_axes([0.03, 0.55, 0.05, 0.1], facecolor=widget_color)
    data_chooser = RadioButtons(data_chooser_axes, vals.keys())
    
    mc_scale_slider_axes = fig_z.add_axes([0.04, 0.10, 0.03, 0.35], facecolor=widget_color)
    mc_scale_slider = VertSlider(mc_scale_slider_axes, r'MC scaling', 0.0, -vals_min,
                                 valinit=-vals_min/3)
    
    def update_plot(_):
        if psi2_slider != None:
            curr_psi2_index = int(psi2_slider.val)
            curr_psi2_val = phis[curr_psi2_index]
            psi2_slider.valtext.set_text(psi2_slider.valfmt % (np.rad2deg(curr_psi2_val)))
        else:
            curr_psi2_index = 0
        curr_r_index = int(r_slider.val)
        curr_r_val = rs[curr_r_index]
        r_slider.valtext.set_text(r_slider.valfmt % curr_r_val)
        curr_data = data_chooser.value_selected
        if curr_data == 'md':
            to_plot = vals['md']
        else:
            to_plot = vals['mc']*(3.0/-vals_min)*mc_scale_slider.val 
            if curr_data == 'diff':
                to_plot = vals['md'] - to_plot
        for i in range(theta_points):
            lines[i][0].set_ydata(to_plot[curr_psi2_index][i][curr_r_index])
        fig_z.canvas.draw_idle()
    
    if psi2_slider != None:
        psi2_slider.on_changed(update_plot)
    r_slider.on_changed(update_plot)
    data_chooser.on_clicked(update_plot)
    mc_scale_slider.on_changed(update_plot)
        
    return psi2_slider, r_slider, data_chooser, mc_scale_slider

#======================================================================================

def calculate_point_rod(bead_type, rod_state):
    '''
    Calculates the interaction between a bead of a given type at (z,r) and a rod
    in a given state at (0,0), for all zs and rs and internal angles psi1.
    '''
    vals = {}
    
    vals['md'] = np.array([md.point_rod(bead_type, rod_state, r, z, psi1)
                              for psi1 in thetas for r in rs for z in zs])
    vals['md'] = vals['md'].reshape(theta_points, r_points, z_points)
    vals_min = vals['md'].min()
    print "min MD val = {}".format(vals_min)
    
    vals['mc'] = np.array([mc.sb_interaction(r, z, psi1)
                              for psi1 in thetas for r in rs for z in zs])
    vals['mc'] = vals['mc'].reshape(theta_points, r_points, z_points)
    vals['mc'] = vals['mc']*(-vals_min)/2
    
    vals['diff'] = vals['md'] - vals['mc']
    
    return vals, vals_min
    
def calculate_rod_rod(rod1_state, rod2_state, theta = 0, phi = 0, psi2 = None):
    '''
    Calculates the interaction between a rod in one state at (z,r) with orientation
    (theta,phi) and internal angle rotation psi2, and a rod in another state at (0,0),
    for all zs and rs and internal angles psi1.
    
    theta : polar angle of rod at (r,z)
    phi : azimuthal angle of rod at (r,z)
    psi2 : angle of patch of rod at (r,z); if "None" calculated for psi2 in [-pi,pi>
    '''
    if psi2 != None:
        psi2s = [psi2]
        psi2_points = 1
    else:
        psi2s = phis
        psi2_points = phi_points
    
    vals = {}
    
    vals['md'] = np.array([md.rod_rod(rod1_state, rod2_state, r, z, theta, phi, psi1, psi2)
                              for psi2 in psi2s for psi1 in thetas for r in rs for z in zs])
    vals['md'] = vals['md'].reshape(psi2_points, theta_points, r_points, z_points)
    vals_min = vals['md'].min()
    print "min MD val = {}".format(vals_min)
    
    vals['mc'] = np.array([mc.bb_interaction(r, z, theta, phi)
                              for r in rs for z in zs])
    vals['mc'] = vals['mc'].reshape(r_points, z_points)
    vals['mc'] = np.array([[vals['mc'] for psi1 in thetas] for psi2 in psi2s]) # so MC data has the same "shape"
    vals['mc'] = vals['mc']*(-vals_min)/3
    
    vals['diff'] = vals['md'] - vals['mc']
    
    return vals, vals_min

#=======================================================================================

import argparse

parser = argparse.ArgumentParser(description=
'''An application for interactive visualisation of interaction potentials for
rod models of the "lammps_multistate_rods" library''',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('config_file', type=str,
                    help='path to the "lammps_multistate_rods" model config file')
parser.add_argument('--rmin', type=float, default=0.0,
                        help='lower bound for the "r" variable')
parser.add_argument('--rmax', type=float, default=5.0,
                        help='upper bound for the "r" variable')
parser.add_argument('--zmin', type=float, default=0.0,
                        help='lower bound for the "z" variable')
parser.add_argument('--zmax', type=float, default=7.5,
                        help='upper bound for the "z" variable')
parser.add_argument('--dx', type=float, default=0.1,
                        help='the grid spacing (in both "r" and "z")')
parser.add_argument('--da', type=float, default=30,
                        help='the angle step (for all angles; in deg)')
args = parser.parse_args()

if __name__ != '__main__':
    print "(visualise.py) ERROR: This module should only be called as a stand-alone application!"
    quit()

cfg_filename = args.config_file
model = Model(cfg_filename)
mc.setup(model)
md.setup(model)

# plot parameters and grid points
dx = args.dx*model.rod_radius
da = np.deg2rad(args.da)

rmin = args.rmin*model.rod_radius; rmax = args.rmax*model.rod_radius
r_points = int((rmax-rmin)/dx) + 1; zero_r = 0
rs = np.linspace(rmin, rmax, r_points)

zmin = args.zmin*model.rod_radius; zmax = args.zmax*model.rod_radius
z_points = int((zmax-zmin)/dx) + 1; zero_z = 0
zs = np.linspace(zmin, zmax, z_points)

phimin = -np.pi; phimax = np.pi
phi_points = int((phimax-phimin)/da); zero_phi = phi_points/2
phis = np.linspace(phimin, phimax, phi_points, endpoint=False)

thetamin = 0; thetamax = np.pi
theta_points = int((thetamax-thetamin)/da) + 1; zero_theta = 0
thetas = np.linspace(thetamin, thetamax, theta_points)

# figure parameters
axis_font = {'size':13}
img_cmap = plt.get_cmap("RdBu")
plot_cmap = plt.cm.get_cmap('nipy_spectral')
widget_color = 'lightgoldenrodyellow'

# interactive console - choice of interaction
while True:
    print
    int_type = raw_input("Enter '1' for point-rod or '2' for rod-rod interaction: ")
    int_type = int_type.strip()
    if int_type == '1':
        while True:
            bead_type = raw_input("Enter type of bead at (z,r): ")
            try:
                bead_type = int(bead_type)
                if bead_type not in model.all_bead_types:
                    raise Exception("")
                break
            except:
                print "Unexpected input ({:}), please enter one of the following: ".format(bead_type),
                print model.all_bead_types
        for i in range(model.num_states):
                print "{:2d} : {:s}".format(i, model.rod_states[i])
        while True:
            rod_state = raw_input("Enter state ID of rod at (0,0): ")
            try:
                rod_state = int(rod_state)
                if rod_state < 0 or rod_state >= model.num_states:
                    raise Exception("")
                break
            except:
                print "Unexpected input ({:}), please try again...".format(rod_state)
        print "Calculating..."
        vals, vals_min = calculate_point_rod(bead_type, rod_state)
        break
    elif int_type == '2':
        for i in range(model.num_states):
            print "{:2d} : {:s}".format(i, model.rod_states[i])
        while True:
            rod1_state = raw_input("Enter state ID of rod at (0,0): ")
            try:
                rod1_state = int(rod1_state)
                if rod1_state < 0 or rod1_state >= model.num_states:
                    raise Exception("")
                break
            except:
                print "Unexpected input ({:}), please try again...".format(rod1_state)
        while True:
            rod2_state = raw_input("Enter state ID of rod at (z,r): ")
            try:
                rod2_state = int(rod2_state)
                if rod2_state < 0 or rod2_state >= model.num_states:
                    raise Exception("")
                break
            except:
                print "Unexpected input ({:}), please try again...".format(rod2_state)
        while True:
            theta = raw_input("Enter angle from z axis (theta; in deg) of rod at (z,r): ")
            if theta.strip() == '':
                theta = 0
                break
            else:
                try:
                    theta = np.deg2rad(float(theta))
                    break
                except:
                    print "Unexpected input ({:}), please try again...".format(theta)
        while True:
            phi = raw_input("Enter angle around z axis (phi; in deg) of rod at (z,r): ")
            if phi.strip() == '':
                phi = 0
                break
            else:
                try:
                    phi = np.deg2rad(float(phi))
                    break
                except:
                    print "Unexpected input ({:}), please try again...".format(phi)
        while True:
            psi2 = raw_input("Enter internal rotation (psi2; in deg) of rod at (z,r): ")
            if psi2.strip() == '':
                psi2 = 0
                break
            else:
                try:
                    psi2 = np.deg2rad(float(psi2))
                    break
                except:
                    print "Unexpected input ({:}), please try again...".format(psi2)
        print "Calculating..."
        vals, vals_min = calculate_rod_rod(rod1_state, rod2_state, theta, phi, psi2)
        break
    else:
        print "Unexpected input ({:s}), please try again...".format(int_type)

# interactive console - choice of plot
while True:
    print
    print "1 : 2D plot"
    print "2 : z-slice plot"
    print "3 : r-slice plot"
    print "q : Quit"
    plot_type = raw_input("What to plot? (enter 1, 2 or 3) ")
    plot_type = plot_type.strip()
    if plot_type == '1':
        if int_type == '1':
            widgets = draw_point_rod_2D(vals, vals_min,
                                        "point({:d})-rod({:s}) interaction  ({:s})".format(
                                            bead_type, model.rod_states[rod_state],
                                            cfg_filename))
        else:
            widgets = draw_rod_rod_2D(vals, vals_min,
                                      "rod({:s})-rod({:s}) interaction  ({:s})".format(
                                            model.rod_states[rod1_state],
                                            model.rod_states[rod2_state],
                                            cfg_filename))
    elif plot_type == '2':
        if int_type == '1':
            widgets = draw_point_rod_z_slice(vals, vals_min,
                                        "point({:d})-rod({:s}) interaction  ({:s})".format(
                                            bead_type, model.rod_states[rod_state],
                                            cfg_filename))
        else:
            widgets = draw_rod_rod_z_slice(vals, vals_min,
                                      "rod({:s})-rod({:s}) interaction  ({:s})".format(
                                            model.rod_states[rod1_state],
                                            model.rod_states[rod2_state],
                                            cfg_filename))
    elif plot_type == '3':
        if int_type == '1':
            widgets = draw_point_rod_r_slice(vals, vals_min,
                                        "point({:d})-rod({:s}) interaction  ({:s})".format(
                                            bead_type, model.rod_states[rod_state],
                                            cfg_filename))
        else:
            widgets = draw_rod_rod_r_slice(vals, vals_min,
                                      "rod({:s})-rod({:s}) interaction  ({:s})".format(
                                            model.rod_states[rod1_state],
                                            model.rod_states[rod2_state],
                                            cfg_filename))
    elif plot_type.lower() == 'q':
        quit()
    else:
        print "Unexpected input ({:s}), please try again...".format(plot_type)
        continue
    plt.show()
