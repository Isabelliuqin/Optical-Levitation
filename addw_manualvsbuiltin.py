# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:26:22 2020

@author: liuqi
"""


import scipy as sp
import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pylab as plt
import scipy.integrate as spi
from scipy.integrate import quad
import seaborn
from scipy.integrate import odeint
from scipy.integrate import dblquad
import Will_Module_addwdep as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM

import time

#######################################################################
# Figure 3.2, Comparison between manual integration method and dblquad
#######################################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 400

plt.close('all')

###########################


#Our sphere

g = 9.8 
c = 3 * 10**8

w_0 = 10 * 10 ** (-6)

Lambda = 1.064 * 10**(-6)

z_R = np.pi* w_0 ** 2 / Lambda

rho = 30 * 10 ** (-6)
 
n_0 = 1

n_s_n = 0.04

k = 7.6097

n_s = n_s_n - k*1j

sig_s = 10.49 * 10 ** 3 * ( (3 ** 3 - 2.25 ** 3) / 3 ** 3 ) #density of sphere in kg/m^3

sig_0 = 0 #density of medium in kg/m^3

m = 4/3 * np.pi * rho ** 3 * ( sig_s - sig_0 )

Permittivity = 8.85 * 10**(-12)

P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam



##################Qy vs rho_0y#############################


######################################################
#plot Qy vs rho_0y at rho = 15w
######################################################

w = 10 * 10 ** (-6)

rho_0 = [0 , 0]

rho_0[1] = np.linspace(-30 * w/np.sqrt(2), 30*w/np.sqrt(2), 100)

rho = 15*w

radial_flisty4 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'manual', grid_size = grid_size))

Q_y4 = radial_flisty4 * c / ( n_0 * P )


plt.figure(15)

plt.plot(np.sqrt(2)*rho_0[1] / w, Q_y4, lw=2, c="c", label="rho = 150um, w = 10um")


new_ticks1 = np.linspace(-30, 30, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 1, 5),fontsize=20)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-30))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=20)

plt.xlabel('rho_0y/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qy',fontsize=20)


plt.grid()
plt.show()

