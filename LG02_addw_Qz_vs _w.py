# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 12:32:30 2020

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
import Will_Module_addwdep_LG02 as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM

import time

############################
# integration method control
############################

integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 400

plt.close('all')

###########################

#Our sphere

g = 9.8 

c = 3 * 10**8

w_0 = 0.73 * 10 ** (-6)

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




#############################################################
#plot of Q_z vs displacement w for various offset rho_0x
#############################################################

rho_0 = [0,0]   #no offset

rho_0[0] = [0, 0.25*rho, 0.5*rho, 0.75*rho, rho, 1.25*rho, 2*rho]

w = np.linspace(w_0, 1.5*rho, 500)

Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )
Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][3],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )

Axial_flist_vs_d4 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][4],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z4 = Axial_flist_vs_d4 * c / ( n_0 * P )

Axial_flist_vs_d5 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][5],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z5 = Axial_flist_vs_d5 * c / ( n_0 * P )

Axial_flist_vs_d6 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][6],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z6 = Axial_flist_vs_d6 * c / ( n_0 * P )

plt.figure(1)


plt.plot( w / (rho), Q_z0 , lw=2, c="c", label="rho_0x/rho = 0")
plt.plot( w / (rho), Q_z1 , lw=2, c="r", label="rho_0x/rho = 0.25")
plt.plot( w / (rho), Q_z2 , lw=2, c="g", label="rho_0x/rho = 0.5")
plt.plot( w / (rho), Q_z3 , lw=2, c="y", label="rho_0x/rho = 0.75")
plt.plot( w / (rho), Q_z4 , lw=2, c="b", label="rho_0x/rho = 1")
plt.plot( w / (rho), Q_z5 , lw=2, c="k", label="rho_0x/rho = 1.25")
plt.plot( w / (rho), Q_z6 , lw=2, c="m", label="rho_0x/rho = 2")


new_ticks1 = np.linspace(0, 1.5, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=13)
plt.title('rho = 30um, w0 = 1um',fontsize=20)

plt.xlabel('w/rho',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()






'''
#1 plot of Q_z vs displacement d


rho_0 = [0,0]   #no offset


w = np.linspace(w_0, 60 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


Axial_flist_vs_d =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(3)


plt.plot( d * 10 ** (6), Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 1800, 10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0z(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','0',  '30')
'''