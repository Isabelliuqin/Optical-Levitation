# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:54:15 2019

@author: wrk10

functions to control figure layout on screen and saving

"""
import matplotlib.pyplot as plt
import numpy as np
from subprocess import call

def move_figure(width = None, width_fraction=2/5, height = None, height_fraction = 1, align = 'right'):
    # set figure dimensions
    manager = plt.get_current_fig_manager()

    manager.full_screen_toggle()  # primitive but works to get screen size
    py = manager.canvas.height()
    px = manager.canvas.width()
    manager.full_screen_toggle()

    # assumes QT backend
    geom = manager.window.geometry()

    x, y, wx, wy = geom.getRect()
    
    if not width:
        width = int(px*width_fraction)
    
    if not height:
        height = int(py*height_fraction)-10
    
    if align == 'left':
        newx = 0
    elif align == 'right':
        newx = px-width-10
    else:
        raise ValueError(f"align must be 'left' or 'right' not {align}")
        
    newy = int(py/2-height/2)+25

    manager.window.setGeometry(newx, newy, width, height)
    
def move_figure_to_second_screen():
    manager = plt.get_current_fig_manager()
    
    # assumes QT backend
    geom = manager.window.geometry()
    
    manager.window.showNormal()             # reset window so no clash with maximising if already maximised
    
    x, y, wx, wy = geom.getRect()           # get window geometry
    
    manager.window.setGeometry(-900, y, wx, wy)     # set window position to off left of main screen, (0,0) is top left of main screen
    
    manager.window.showMaximized()          # maximise window in current screen
    
    
def scale_fonts_and_layout():
    
    ax = plt.gca()
    
    ax.xaxis.label.set_size('x-large')
    ax.yaxis.label.set_size('x-large')

    plt.tight_layout()

    plt.draw()
    return
    

    
def save_figure(path, file_format):
    path = str(path)
    
    if file_format != 'emf':
        plt.savefig(path + '.' + file_format,format = file_format)
    else:
    
        temp_svg_file = path + '.svg'
        
        plt.savefig(temp_svg_file, format = 'svg')
        
        call(["C:\Program Files\Inkscape\inkscape.exe", "--file", str(temp_svg_file), "--export-emf", str(temp_svg_file)[:-4] + ".emf"])
    
    
## data formatting

def moduluo_remove_jumps(x, y, modulo, threshold = 0.9):
    # inserts NaN values to remove line jumps between discontinuities
    # x = x data points
    # y = y data points
    # modulo = value for (y % modulo)
    # threshold = fractional jump in modulo value to mask, can be tuned for less smooth data
    #
    # returns (x, y)
    
    x = np.array(x)
    y = np.array(y)
    
    y = np.mod(y, modulo)
    
    nan_insert_indices = np.where(np.abs(np.diff(y)) > modulo*threshold)[0] + 1    # np.where returns tuple of numpy array, we need array. +1 offsets correctly for np.insert
    
    y = np.insert(y, nan_insert_indices, np.NaN)
    x = np.insert(x, nan_insert_indices, np.NaN)
    
    return x, y
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    