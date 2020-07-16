# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:36:29 2020

@author: wrk10_home

fucntions to generate and manipulate 2D coordinate systems

"""

import numpy as np
    
def reduce_to_broadcastable_arrays(x,y):
    
    axis = find_axis_of_change(x)
    
    # select appropriate axis of change
    if axis == 0:
        x_reduced = x[:,0]
        x_reduced = x_reduced[:,np.newaxis]   # slice return 1D array, ensure 2D in correct dimension
        y_reduced = y[0,:]
        y_reduced = y_reduced[np.newaxis,:]
    elif axis == 1:
        x_reduced = x[0,:]
        x_reduced = x_reduced[np.newaxis,:]
        y_reduced = y[:,0]
        y_reduced = y_reduced[:,np.newaxis]
    else:
        raise ValueError('Cannot reduce coordinate grid, invalid axis of change')
    
    return [x_reduced, y_reduced]       

def get_spacing(x):
    axis = find_axis_of_change(x)
    
    if axis == 0:
        dx = np.abs(x[1,0] - x[0,0])
    elif axis == 1:
        dx = np.abs(x[0,1] - x[0,0])
    else:
        raise ValueError('Invalid axis of change, must be 0 or 1')
    
    return dx

def find_axis_of_change(x):
    #   finds the axis along which a coordinate matrix/vector changes
    
    if x.ndim != 2:
        raise ValueError('coordinate system must be at least 2D')
    
    x_shape = np.array(x.shape)
    
    if 1 in x_shape:
        axis = np.argmax(x_shape)
    elif x[0,0] != x[-1,0]:
        axis = 0
    else:
        axis = 1    
    
    return axis

def get_grid_xy(x_width, y_width, n_pixels_x, n_pixels_y):
    
    x_min = -x_width / 2 
    x_max = x_width / 2
    y_min = -y_width / 2
    y_max = y_width / 2

    x, y = np.mgrid[x_min:x_max:(n_pixels_x*1j), y_min:y_max:(n_pixels_y*1j)]
    
    return x, y

def get_grid_xy_square(side_length, n_pixels = 100):
    return get_grid_xy(side_length, side_length, n_pixels, n_pixels)

def get_grid_open_xy(x_width, y_width, n_pixels_x, n_pixels_y):
    
    x_min = -x_width / 2 
    x_max = x_width / 2
    y_min = -y_width / 2
    y_max = y_width / 2

    x,y = np.ogrid[x_min:x_max:(n_pixels_x*1j), y_min:y_max:(n_pixels_y*1j)]
    
    return x, y

def get_grid_open_xy_square(side_length, n_pixels = 100):
    return get_grid_open_xy(side_length, side_length, n_pixels, n_pixels)

def used_mgrid(x):
    
    axis = find_axis_of_change(x)
    
    if axis == 0:
        return True
    elif axis == 1:
        return False
    else:
        raise ValueError('Cannot determine if grid used mgrid')
        
def find_closest_value_index(value, array):
    return np.argmin(np.abs(array-value))

def find_coordinate_expected_index(value, coordinate):
    # find expected index of a coordinate value in coordinate array, will
    # extrapolate outside of range
    # value - find the index of this value in the coordinate array
    # coordinate - coordinate array to find value in
    #
    # will return negative indices
    
    if coordinate.shape[0] == 1 and coordinate.shape[1] == 1:
        raise ValueError(f'To use pad_coordinate array must be in broadcastable form, currently: rows {coordinate.shape[0]}, columns{coordinate.shape[1]}')
        return
    
    if coordinate[-1, -1] < coordinate[0, 0]:
        raise ValueError('To use pad_coordinate array must be increasing monatonically')
        
    dx = get_spacing(coordinate)
    
    index = np.round( (value - coordinate.min()) / dx )
    
    return np.int(index)
    

def pad_coordinate(x, pad_number, pad_location = 'both'):
    # pad coordinate array by number pad_number. Assume monatonically increasing 
    # and single row/column, i.e. used in broadcasting mode
    # x - coordinate to pad, must be in broadcastable form
    # pad_number - number of elements to pad coordinate on each end
    # pad_location - 'lower' - pad lower range of coordinate
    #                'upper' - pad upper range of coordinate
    #                'both' - pad both ends of range
    
    
    if x.shape[0] == 1 and x.shape[1] == 1:
        raise ValueError(f'To use pad_coordinate array must be in broadcastable form, currently: rows {x.shape[0]}, columns{x.shape[1]}')
        return
    
    if x[-1, -1] < x[0, 0]:
        raise ValueError('To use pad_coordinate array must be increasing monatonically')
    
    pad_number_lower = []
    pad_number_upper = []
    
    if pad_location == 'both':
        pad_number_lower = pad_number
        pad_number_upper = pad_number
    elif pad_location == 'upper':
        pad_number_lower = 0
        pad_number_upper = pad_number
    elif pad_location == 'lower':
        pad_number_lower = pad_number
        pad_number_upper = 0
    else:
        raise ValueError(f"pad_location must be 'lower', 'upper', or 'both', not {pad_location}")
        
    dx = get_spacing(x)
    
    pad_array_lower = dx * np.arange(-pad_number_lower, 0) + x.min()                              # better to use np.arange with integers
    pad_array_upper = dx * (np.flip(np.arange(0, -pad_number_upper, -1)) + pad_number) + x.max()      # do this to ensure max coordinate value not repeated
    
    pad_array_lower = np.atleast_2d(pad_array_lower)
    pad_array_upper = np.atleast_2d(pad_array_upper)
    
    axis = find_axis_of_change(x)
    
    if axis == 0:
        pad_array_lower = pad_array_lower.T
        pad_array_upper = pad_array_upper.T
    
    x = np.concatenate((pad_array_lower, x, pad_array_upper), axis = axis)
    
    return x

def find_coordinate_pad_number(x, x_range):
    # returns number of elements to add to coordinate given a required total range
    # will add elements symettrically on each side
    # x - coordinate to pad
    # x_range - required range of x
    
    if x.shape[0] == 1 and x.shape[1] == 1:
        raise ValueError(f'To use pad_coordinate array must be in broadcastable form, currently: rows {x.shape[0]}, columns{x.shape[1]}')
        return
    
    if x[-1, -1] < x[0, 0]:
        raise ValueError('To use pad_coordinate array must be increasing monatonically')
    
    range_to_add = x_range - x.ptp()
    
    if range_to_add < 0:
        return 0
    
    dx = get_spacing(x)
    
    return np.int(np.ceil(range_to_add / 2 / dx))

def check_coordinates_consistent(x, y, array):
    # check that array sizes of x and y coordinates match intended array
    
    x, y = reduce_to_broadcastable_arrays(x, y)
    
    if array.shape == (x.size, y.size) or array.shape == (y.size, x.size):
        return True
    else:
        raise ValueError(f'shape of array does not match coordinates, x: {x.shape}, y: {y.shape}, array: {array.shape}')
    