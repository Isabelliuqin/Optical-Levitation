# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 17:59:54 2020

@author: wrk10_home
module to aid in plotting figures in sensible SI units
"""

import numpy as np
from matplotlib.ticker import FuncFormatter

def rescale_axis(axis, unit_seperator = '/'):
    # rescales single axis into sensible SI units, modifies ticks and label
    # unit_seperator - str of method to separate label and units, e.g. x / m, used to seperate label from units
    # notes: currently MUST have a seprator and units, otherwise will get confused finding label and units
    
    axis_label_original = axis.get_label().get_text()      # get original label string, first get label text object, then label within that
    
    axis_name, units_original = axis_label_original.split(unit_seperator)
    axis_name = axis_name.strip()
    units_original = units_original.strip(' )')  # remove whitespace and trailing parenthesis if () unit_separator used
    
    tick_labels_original = axis.get_majorticklocs()
    
    multiplier, new_units = get_scale_units(np.abs(tick_labels_original).max(), units_original)
    
    tick_label_formatter = FuncFormatter(lambda value, position: '{:.0f}'.format(multiplier*value))
    
    axis.set_major_formatter(tick_label_formatter)
    axis.set_minor_formatter(tick_label_formatter)
    
    axis.set_label_text(axis_name + ' / ' + new_units)
    
def rescale_all_axes(axes, unit_seperator = '/'):
    # rescales both axes on plot into sensible SI units, modifies ticks and label, see rescale_axis
    # unit_seperator - str of method to separate label and units, e.g. x / m, used to seperate label from units
    # notes: currently MUST have a seprator and units, otherwise will get confused finding label and units
    
    rescale_axis(axes.xaxis)
    rescale_axis(axes.yaxis)

def get_scale_units(max_value, units):
    # get modifiers to make plot scale  'nice' multiples of SI units
    #   max_value = maximum (or desired max) of plot scale
    #   units = str of axes units
    # returns:
    #   multiplier = multiplier for axes values
    #   units = new unit value
    
    SI_table = {12: 'T',
                 9: 'G',
                 6: 'M',
                 3: 'k',
                 0: '',
                 -3: 'm',
                 -6: '$\mu$',
                 -9: 'n',
                 -12: 'p',
                 -15: 'f'}
    
    magnitude = 3*np.floor(np.log10(np.abs(max_value))/3)   # find closest multiple of 3 power  of 10
    
    multiplier = []  # value to multiply data by
    unit_prefix = [] # str to add to units
    
    if magnitude in SI_table.keys():
        multiplier = 10**-magnitude
        unit_prefix = SI_table[magnitude]
    else:
        multiplier = 1
        unit_prefix = ''
    
    new_units = unit_prefix + units
    
    return multiplier, new_units

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    
    x = np.linspace(0, 10000)
    y = (x - 333)**2
    
    fig = plt.figure('SI_axes_testing')
    axes = plt.gca()
    
    plt.plot(x, y)
    
    plt.xlabel('x / m')
    plt.ylabel('y / s')
    
    rescale_all_axes(axes)
    
    
    
    
    
    
    
    