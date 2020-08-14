# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 18:30:33 2020

@author: liuqi
"""


import matplotlib.pyplot as plt
from platform import python_version as pythonversion
from matplotlib import __version__ as matplotlibversion

print('python: '+pythonversion())
print('matplotlib: '+matplotlibversion)


w_0 = 2.67
w = w_0
rho_0x = 0
rho_0y = 0
rho_0z = 0
rho = 30

def table_parameter(w, rho_0z, rho_0y, rho_0x, rho):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    y = [1, 2, 3, 4, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1]    
    col_labels = ['VALUE']
    row_labels = ['wavelength(um)', 'w_0(um)', 'n_0', 'n_s real', 'n_s imag', 'beam radius w(um)', 'rho_0z(um)', 'rho_0y(um)', 'rho_0x(um)', 'target radius rho(um)']
    table_vals = [[1.064], [10], [1], [0.04], [-7.6097], [w], [rho_0z], [rho_0y], [rho_0x], [rho]]

# Draw table
    the_table = plt.table(cellText=table_vals,
                          colWidths=[0.08] * 8,
                          rowLabels=row_labels,
                          colLabels=col_labels,
                          loc='center')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(24)
    the_table.scale(4, 4)

    # Removing ticks and spines enables you to get the figure only with table
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    for pos in ['right','top','bottom','left']:
        plt.gca().spines[pos].set_visible(False)
        plt.savefig('matplotlib-table.png', bbox_inches='tight', pad_inches=0.05)
        
        