# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:23:55 2020

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
import Will_Module_addwdep_LG03 as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM

import Module_LG03_Fx_stiffness_at_parameter as FTAP

import time

##########################################################
# integration method control: Qx vs w plots
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 400

plt.close('all')

###########################

#Our sphere

g = 9.8

c = 3 * 10**8

w_0 = 0.7* 10 ** (-6)


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

#P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam

P = 23.21



#########################################
#11 plot Qx vs rho_0x for various w
#########################################




rho_0 = [0,0]   #no offset



w = [w_0, 0.5*np.sqrt(2/3) *rho, 0.55 *np.sqrt(2/3) * rho, 0.6 *np.sqrt(2/3) * rho, 0.75 *np.sqrt(2/3) * rho, np.sqrt(2/3) *rho, 1.25*np.sqrt(2/3) *rho, 1.5*np.sqrt(2/3) *rho]

rho_00 = np.linspace(-4*np.sqrt(3/2)*w[0], 4*np.sqrt(3/2)*w[0], 100)

rho_01 = np.linspace(-4*np.sqrt(3/2)*w[1], 4*np.sqrt(3/2)*w[1], 100)

rho_02 = np.linspace(-4*np.sqrt(3/2)*w[2], 4*np.sqrt(3/2)*w[2], 100)

rho_03 = np.linspace(-4*np.sqrt(3/2)*w[3], 4*np.sqrt(3/2)*w[3], 100)

rho_04 = np.linspace(-4*np.sqrt(3/2)*w[4], 4*np.sqrt(3/2)*w[4], 100)

rho_05 = np.linspace(-4*np.sqrt(3/2)*w[5], 4*np.sqrt(3/2)*w[5], 100)

rho_06 = np.linspace(-4*np.sqrt(3/2)*w[6], 4*np.sqrt(3/2)*w[6], 100)

rho_07 = np.linspace(-4*np.sqrt(3/2)*w[7], 4*np.sqrt(3/2)*w[7], 100)


radial_flistx0 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_00,rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0 * c / ( n_0 * P )


radial_flistx1 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_01,rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x1 = radial_flistx1* c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_02,rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x2 = radial_flistx2 * c / ( n_0 * P )

radial_flistx3 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_03,rho_0[1], rho, n_0, n_s, w_0, w[3], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x3 = radial_flistx3 * c / ( n_0 * P )

radial_flistx4 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_04,rho_0[1], rho, n_0, n_s, w_0, w[4], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x4 = radial_flistx4* c / ( n_0 * P )

radial_flistx5 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_05,rho_0[1], rho, n_0, n_s, w_0, w[5], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x5 = radial_flistx5 * c / ( n_0 * P )

radial_flistx6 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_06,rho_0[1], rho, n_0, n_s, w_0, w[6], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x6 = radial_flistx6 * c / ( n_0 * P )

radial_flistx7 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_07,rho_0[1], rho, n_0, n_s, w_0, w[7], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x7 = radial_flistx7 * c / ( n_0 * P )



plt.figure(11)
#plt.plot(rho_00 / (np.sqrt(3/2)*w[0]), Q_x0, lw=2, c="c", label="w = w0")
plt.plot(rho_01 / (np.sqrt(3/2)*w[1]), Q_x1, lw=2, c="r", label="sqrt(3/2)*w/rho = 0.5")
plt.plot(rho_02 / (np.sqrt(3/2)*w[2]), Q_x2, lw=2, c="g", label="sqrt(3/2)*w/rho = 0.55")
plt.plot(rho_03 / (np.sqrt(3/2)*w[3]), Q_x3, lw=2, c="y", label="sqrt(3/2)*w/rho = 0.6")
plt.plot(rho_04 / (np.sqrt(3/2)*w[4]), Q_x4, lw=2, c="b", label="sqrt(3/2)*w/rho = 0.75")
plt.plot(rho_05 / (np.sqrt(3/2)*w[5]), Q_x5, lw=2, c="k", label="sqrt(3/2)*w/rho = 1")
plt.plot(rho_06 / (np.sqrt(3/2)*w[6]), Q_x6, lw=2, c="m", label="sqrt(3/2)*w/rho = 1.25")
plt.plot(rho_07 / (np.sqrt(3/2)*w[7]), Q_x7, lw=2, c="c", label="sqrt(3/2)*w/rho = 1.5")

new_ticks1 = np.linspace(-4, 4, 9) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.3, 0.3, 7),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-4))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/(sqrt(3/2)*w)',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('rho = 30um, w0 = 2um',fontsize=15)
plt.grid()
plt.show()




#####################################
#gradient of Qx at (0,0) for various w
#####################################



rho_0 = [0 , 0]

rho_0x = 0

resolution = 0.01 * 10 ** (-6)

w = np.linspace(0.1 *np.sqrt(2/3) *rho, 1.5 *np.sqrt(2/3) *rho, 1000)

grad_Fx = FTAP.Fx_stiffness_vs_w_plots(rho_0x,rho_0[0], rho, n_0, n_s, w_0, w, z_R, P, resolution, target = "reflective", integration_method = integration_method, grid_size = grid_size)

plt.figure(13)
plt.plot(np.sqrt(3/2)*w/rho, np.array(grad_Fx) * 10 ** 4, lw=2, c="c", label="rho_0x = 0")

print ((np.sqrt(3/2)*w/rho)[np.argmin(abs(np.array(grad_Fx)))]) #print the inflection point


new_ticks1 = np.linspace(0.4, 1.6, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-10, 5, 4),fontsize=20)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data', 0.4))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('sqrt(3/2)w/rho',fontsize=20)
plt.ylabel('grad_x(stiffness)10^(-8)',fontsize=20)

plt.title('rho = 30um, w0 = 2um',fontsize=15)
plt.grid()
plt.show()
