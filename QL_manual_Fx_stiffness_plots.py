# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:25:13 2020

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
# integration method control: stiffness plot of Qx, Qz
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
###########################################
#13 plot gradient Fx vs rho_0x for various rho
###########################################

a = 30 * 10 ** (-6)
w = a*np.sqrt(2)

rho_0 = [0 , 0]

rho_0[0] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_x0 = np.asarray(TQ.Fx_total_gradient(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


grad_x1 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method  = integration_method, grid_size = grid_size)['Fx_grad'])

grad_x2 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


rho_0x = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))


plt.figure(13)
plt.plot( rho_0x, grad_x0*10**8, lw=2, c="c", label="rho = 15um")

plt.plot( rho_0x, grad_x1*10**8, lw=2, c="r", label="rho = 22.5um")

plt.plot( rho_0x, grad_x2*10**8, lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-3, 3, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('grad_Fx(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30', 'related to w', rho_0[1]* 10 ** 6, 'x-axis -30 to 30', '15, 22.5, 30')
'''



#############################################
#13 plot gradient Fx vs rho_0x for various w
#############################################
a = 30 * 10 ** (-6)
w = [a*np.sqrt(2), 2*a, 2.5*a]

rho_0 = [0 , 0]

rho_0[0] = np.linspace(-a, a, 100)

rho = a


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_x0 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


grad_x1 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])

grad_x2 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


rho_0x = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))


plt.figure(13)
plt.plot( rho_0x, grad_x0*10**8, lw=2, c="c", label="w = sqrt(2)*30um")

plt.plot( rho_0x, grad_x1*10**8, lw=2, c="r", label="w = 60um")

plt.plot( rho_0x, grad_x2*10**8, lw=2, c="g", label="w = 75um")


new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-3, 3, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('grad_Fx(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30, 60, 75', 'related to w', rho_0[1]* 10 ** 6, 'x-axis -30 to 30', '30')
