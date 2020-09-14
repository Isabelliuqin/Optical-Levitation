# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 13:10:25 2020

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

import Module_LG02_Fz_stiffness_at_parameter as FTAPz
import Module_LG02_Fx_stiffness_at_parameter as FTAPx

import time

##########################################################
# integration method control: stiffness plot of Qx, Qz
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 400

plt.close('all')

###########################

#Our sphere

g = 9.8 

c = 3 * 10**8

w_0 = 0.894 * 10 ** (-6)

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

#P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam

P = 15.6


#############################################
#13 plot gradient Fz vs w for various offsets
#############################################

w = np.linspace(w_0, 3*np.sqrt(2) * rho, 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0 , 0]

#rho_0[0] = [-w_0/np.sqrt(2), -0.5*w_0, 0, 0.5*w_0, w_0/np.sqrt(2)]




#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_z0 = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


#grad_z1 = np.asarray(TQ.Fz_total_gradient(rho_0[0][1], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])

#grad_z2 = np.asarray(TQ.Fz_total_gradient(rho_0[0][2], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])

#grad_z3 = np.asarray(TQ.Fz_total_gradient(rho_0[0][3], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fz_grad'])


rho_0z = np.asarray(TQ.Fz_total_gradient(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0z'])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))

w_in = w_0 * np.sqrt(1 + (rho_0z/z_R)**2)


plt.figure(13)
plt.plot( w_in/(np.sqrt(2) * rho), -grad_z0*10**8, lw=2, c="c", label="rho_0x = 0")#

#plt.plot( w * 10 ** 6, grad_z1*10**8, lw=2, c="r", label="rho_0x = -5um")

#plt.plot( w * 10 ** 6, grad_z2*10**8, lw=2, c="g", label="rho_0x = 5um")

#plt.plot( w * 10 ** 6, grad_z3*10**8, lw=2, c="y", label="rho_0x = 7.07um")


new_ticks1 = np.linspace(0, 3, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(5, -20, 6),fontsize=20)
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('w/(sqrt(2)*rho)',fontsize=20)
plt.ylabel('kz(10^(-8)N/m)',fontsize=20)

plt.title('rho = 30um, w0 = 0.894um(stiffness matched)',fontsize=20)
plt.grid()
plt.show()



###################################################
#Fz gradient vs w0 for various w/sqrt(2) = rho
###################################################

w = rho

rho_0 = [0 , 0]

resolution = 0.01 * 10 ** (-6)

w_0 = np.linspace(0.5*10**(-6), 0.73*10**(-6),100)
#z_R = np.pi* w_0 ** 2 / Lambda

grad_z0 = np.asarray(FTAPz.Fz_stiffness_vs_w0_plots(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, resolution, target = 'reflective', integration_method = integration_method, grid_size = grid_size))

grad_x0 = np.asarray(FTAPx.Fx_stiffness_vs_w0_plots(rho_0[0], rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, resolution, target = 'reflective', integration_method = integration_method, grid_size = grid_size))
    

plt.plot( w_0 * 10 ** 6, (-grad_z0*10**4 / (-grad_x0*10**4)) - 1, lw=2, c="c", label="w/(sqrt(2)rho) = 1")
print ((w_0 * 10 ** 6)[np.argmin(abs((-grad_z0*10**8 / (-grad_x0*10**8)) - 1))] ) 

new_ticks1 = np.linspace(0, 20 , 6) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1.5, 1.5, 7),fontsize=20)
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('w_0(um)',fontsize=20)
plt.ylabel('kz / kx - 1',fontsize=20)

plt.title('rho = 30um',fontsize=20)
plt.grid()
plt.show()

