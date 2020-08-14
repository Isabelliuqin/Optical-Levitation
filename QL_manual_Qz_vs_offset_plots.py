# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 10:58:43 2020

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
import Will_Module_Trialplot_Q as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM

import time

##########################################################
# integration method control: Qz vs w plots
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 200

plt.close('all')

###########################

t_start = time.time()


#Our sphere

g = 9.8 #gravitational acceleration
c = 3 * 10**8
#m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)

#TEM01* reflective target Table 6 matching

w_0 = 10* 10 ** (-6)


Lambda = 1.064 * 10**(-6)
#P = 1
z_R = np.pi* w_0 ** 2 / Lambda
rho = 30 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

n_0 = 1

n_s_n = 0.04

k = 7.6097

n_s = n_s_n - k*1j

sig_s = 10.49 * 10 ** 3 * ( 3 ** 3 - 2.9 ** 3 / 3 ** 3 ) #density of sphere in kg/m^3

sig_0 = 0 #density of medium in kg/m^3



m = 4/3 * np.pi * rho ** 3 * ( sig_s - sig_0 )

Permittivity = 8.85 * 10**(-12)

P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam


'''
######################################
#5 plot of Q_z vs rho_0x for various w
######################################

rho = 30 * 10 ** (-6)

rho_0 = [0,0]   #no offset

rho_0[0] = np.linspace(-rho, rho, 100)

w = [np.sqrt(2)*rho, 2*rho, 2.5*rho]


Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )

Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )

plt.figure(1)


plt.plot(rho_0[0] * 10 ** 6, Q_z0 , lw=2, c="c", label="w = sqrt(2)*rho")
plt.plot(rho_0[0] * 10 ** 6, Q_z1 , lw=2, c="r", label="w = 2rho")
plt.plot(rho_0[0] * 10 ** 6, Q_z2 , lw=2, c="g", label="w = 2.5rho")


new_ticks1 = np.linspace(-30, 30, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.4, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=13)

plt.xlabel('rho_0x(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30, 60, 75', 'related to w', '0','x-aixs -30 to 30', '30')
'''



########################################
#6 plot of Q_z vs rho_0x for various rho
########################################

a = 30 * 10 ** (-6)

w = np.sqrt(2)*rho

rho_0 = [0,0]   #no offset

rho_0[0] = np.linspace(-rho, rho, 100)

rho = [0.5*a, 0.75*a, a]


Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )

Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )

plt.figure(2)


plt.plot(rho_0[0]/a, Q_z0 , lw=2, c="c", label="rho = 15um")
plt.plot(rho_0[0]/a, Q_z1 , lw=2, c="r", label="rho = 22.5um")
plt.plot(rho_0[0]/a, Q_z2 , lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(-1, 1, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.4, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=13)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30', 'related to w', '0','x-aixs -30 to 30', '15, 22.5, 30')
