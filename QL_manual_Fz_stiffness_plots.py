# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 19:07:25 2020

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

import Module_Fz_stiffness_at_parameter as FTAPz

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

w_0 = 10 * 10 ** (-6)


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
#13 plot gradient Fz vs rho_0z for various rho
###########################################

a = 30 * 10 ** (-6)


rho_0 = [0 , 0]

w = np.linspace(w_0, 60 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho = [0.5*a, 0.75*a, a]


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_z0 = np.asarray(TQ.Fz_total_gradient(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


grad_z1 = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method  = integration_method, grid_size = grid_size)['Fz_grad'])

grad_z2 = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


rho_0z = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0z'])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))


plt.figure(13)
plt.plot( rho_0z * 10 ** 6, grad_z0*10**8, lw=2, c="c", label="rho = 15um")

plt.plot( rho_0z * 10 ** 6, grad_z1*10**8, lw=2, c="r", label="rho = 22.5um")

plt.plot( rho_0z * 10 ** 6, grad_z2*10**8, lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(0,1800,10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 2, 7),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho_0z(um)',fontsize=20)
plt.ylabel('grad_Fz(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('x-axis w0 to 60', 'related to w', '0', '0', '15, 22.5, 30')


'''

#############################################
#13 plot gradient Fz vs w for various offsets
#############################################
a = 30 * 10 ** (-6)


w = np.linspace(w_0, 3*np.sqrt(2) * rho, 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0 , 0]

#rho_0[0] = [-w_0/np.sqrt(2), -0.5*w_0, 0, 0.5*w_0, w_0/np.sqrt(2)]

rho = a


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_z0 = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


#grad_z1 = np.asarray(TQ.Fz_total_gradient(rho_0[0][1], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])

#grad_z2 = np.asarray(TQ.Fz_total_gradient(rho_0[0][2], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])

#grad_z3 = np.asarray(TQ.Fz_total_gradient(rho_0[0][3], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


rho_0z = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0z'])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))

w_in = w_0 * np.sqrt(1 + (rho_0z/z_R)**2)


plt.figure(13)
plt.plot( rho_0z * 10 ** 6, grad_z0*10**8, lw=2, c="c", label="rho_0x = 0")#

#plt.plot( w * 10 ** 6, grad_z1*10**8, lw=2, c="r", label="rho_0x = -5um")

#plt.plot( w * 10 ** 6, grad_z2*10**8, lw=2, c="g", label="rho_0x = 5um")

#plt.plot( w * 10 ** 6, grad_z3*10**8, lw=2, c="y", label="rho_0x = 7.07um")


new_ticks1 = np.linspace(0, 2000, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.5, 6, 14),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('w/(sqrt(2)*rho)',fontsize=20)
plt.ylabel('grad_Fz(stiffness)10^(-8)',fontsize=20)

plt.title('rho = 30um, w0 = 2um',fontsize=20)
plt.grid()
plt.show()


MTP.table_parameter('x-axs w0 to 1.5*sqrt(2)rho', 'related to w', rho_0[1]* 10 ** 6, '0', '30')



'''
####################################################
#Fx_stiffness at certain rho_0x and plotted over Rs
####################################################

a = 30 * 10 ** (-6)


w = 20 * 10** (-6)

#d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0 , 0]

resolution = 0.1 * 10 ** (-6)

#rho_0[0] = [-w_0/np.sqrt(2), -0.5*w_0, 0, 0.5*w_0, w_0/np.sqrt(2)]

rho = np.linspace(10 * 10 ** (-6), a, 100)

Fz_grad = np.asarray(FTAPz.Fz_stiffness_vs_rho_plots(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, resolution, target = 'reflective', integration_method = integration_method, grid_size = grid_size))
                     
                     
plt.figure(13)
plt.plot( rho * 10 ** (6), Fz_grad*10**8, lw=2, c="c", label="w = 20um")




new_ticks1 = np.linspace(10, 30, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, 3, 4),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)


ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',10))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho(um)',fontsize=20)
plt.ylabel('grad_Fz(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('20', 'related to w', rho_0[1]* 10 ** 6, '0', 'x-axis 10 to 30')
'''