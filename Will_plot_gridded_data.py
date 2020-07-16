# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 13:13:17 2020

@author: wrk10_home

functions to plot data in gridded arrays

"""

import numpy as np
import matplotlib.pyplot as plt
import coordinates_2D as coord_2D

from figure_tools.cmap_tools import convert_to_colormap_object

def imshow_rotated(data, x, y, axes = None, cmap = 'viridis', **imshow_kwargs):
    
    if axes is not None:
        plt.sca(axes)
    
    # test if mgrid used (if true transpose data)
    if coord_2D.used_mgrid(x):
        data = data.T
        
    cmap = convert_to_colormap_object(cmap)
    
    x_step = coord_2D.get_spacing(x)
    y_step = coord_2D.get_spacing(y)
    
    axes_limits = [x.min() - x_step/2, x.max() + x_step/2, y.min() - y_step/2, y.max() + y_step/2]   # steps correct for pixel offset from edges
    
    image = plt.imshow(data, origin = 'lower', extent = axes_limits, cmap = cmap, **imshow_kwargs)
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    return image



def imsave_rotated(filename, data, x, y, cmap = 'inferno', vmin = None, vmax = None, **imsave_kwargs):
    # vmin/vmax - min/max values to apply to colormap

    # test if mgrid used (if true transpose data)
    if coord_2D.used_mgrid(x):
        data = data.T
        
    cmap = convert_to_colormap_object(cmap)
        
    plt.imsave(filename, data, origin = 'lower', cmap = cmap, vmin = vmin, vmax = vmax, **imsave_kwargs)
        
    return




def imshow_rotated_complex(data, x, y, axes = None, cmap_angle = 'hsv', magnitude_exponent = 1, **imshow_kwargs):
    # plot complex data where phase angle maps to colour, and magnitude maps to transparency
    # data - complex values to plot
    # x - x coordinates
    # y - y coordinates
    # axes - axes handle to plot into
    # cmap_angle - colormap to use for phase component
    # magnitude exponent - map magnitude**exponent to transparency, i.e. 2 would be squared magnitude (i.e. intensity)
    
    if axes is not None:
        plt.sca(axes)
    
    # test if mgrid used (if true transpose data)
    if coord_2D.used_mgrid(x):
        data = data.T
    
    angle = (np.angle(data) + np.pi) / (2*np.pi)    # normalised angle 0 -> 1
    
    magnitude = np.sqrt(np.real(data*data.conj()))           # normalised magnitude 0 -> 1
    magnitude /= np.max(magnitude)
    magnitude = magnitude**magnitude_exponent
    
    axes_limits = [x.min(), x.max(), y.min(), y.max()]
    
    colormap_object = convert_to_colormap_object(cmap_angle)
    data_RGBA = colormap_object(angle)
        
    for i in range(0,3):
        data_RGBA[:,:,i] *= magnitude
    
    
    image = plt.imshow(data_RGBA, origin='lower', extent = axes_limits, **imshow_kwargs)
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    return image

def set_phase_mapping(image, positive_phase = False, add_colorbar = True):
    # scale colormap between 0->2pi, or -pi -> pi
    # image = handle of image object
    # positive_phase = True for 0->2pi scale, False for -pi->pi scale
    # add_colorbar = add colorbar with formatted ticks or not
    
    if positive_phase:
        clim_min,clim_max = [0, 2*np.pi]
        cticks = np.linspace(clim_min, clim_max, 5, endpoint = True)
        ctick_labels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$\pi/2$', r'$2\pi$']
    else:
        clim_min,clim_max = [-np.pi, np.pi]
        cticks = np.linspace(clim_min, clim_max, 5, endpoint = True)                         
        ctick_labels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']                         

    image.set_clim([clim_min, clim_max])
        
    if add_colorbar:
        cbar = plt.colorbar()
        cbar.set_ticks(cticks)
        cbar.set_ticklabels(ctick_labels)
    
    return
