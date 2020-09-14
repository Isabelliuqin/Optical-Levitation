# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 19:11:59 2020

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
import Module_Fx_stiffness_at_parameter as FTAP

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





#############################################
#13 plot gradient Fx vs rho_0x for various w
#############################################


rho_0 = [0,0]   #no offset



w = [0.5*np.sqrt(2)*rho, 0.65*np.sqrt(2)*rho, 0.8*np.sqrt(2)*rho, np.sqrt(2)*rho, 1.25*np.sqrt(2)*rho, 1.5*np.sqrt(2)*rho]

rho_00 = np.linspace(-4*w[0]/np.sqrt(2), 4*w[0]/np.sqrt(2), 100)

rho_01 = np.linspace(-4*w[1]/np.sqrt(2), 4*w[1]/np.sqrt(2), 100)

rho_02 = np.linspace(-4*w[2]/np.sqrt(2), 4*w[2]/np.sqrt(2), 100)

rho_03 = np.linspace(-4*w[3]/np.sqrt(2), 4*w[3]/np.sqrt(2), 100)

rho_04 = np.linspace(-4*w[4]/np.sqrt(2), 4*w[4]/np.sqrt(2), 100)

rho_05 = np.linspace(-4*w[5]/np.sqrt(2), 4*w[5]/np.sqrt(2), 100)


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_x0 = np.asarray(TQ.Fx_total_gradient(rho_00, rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


grad_x1 = np.asarray(TQ.Fx_total_gradient(rho_01, rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])

grad_x2 = np.asarray(TQ.Fx_total_gradient(rho_02, rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])

grad_x3 = np.asarray(TQ.Fx_total_gradient(rho_00, rho_0[1], rho, n_0, n_s, w_0, w[3], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])


grad_x4 = np.asarray(TQ.Fx_total_gradient(rho_01, rho_0[1], rho, n_0, n_s, w_0, w[4], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])

grad_x5 = np.asarray(TQ.Fx_total_gradient(rho_02, rho_0[1], rho, n_0, n_s, w_0, w[5], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['Fx_grad'])




rho_0x0 = np.asarray(TQ.Fx_total_gradient(rho_00, rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

rho_0x1 = np.asarray(TQ.Fx_total_gradient(rho_01, rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

rho_0x2 = np.asarray(TQ.Fx_total_gradient(rho_02, rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

rho_0x3 = np.asarray(TQ.Fx_total_gradient(rho_03, rho_0[1], rho, n_0, n_s, w_0, w[3], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

rho_0x4 = np.asarray(TQ.Fx_total_gradient(rho_04, rho_0[1], rho, n_0, n_s, w_0, w[4], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])

rho_0x5 = np.asarray(TQ.Fx_total_gradient(rho_05, rho_0[1], rho, n_0, n_s, w_0, w[5], z_R, P, target = 'reflective', integration_method = integration_method, grid_size = grid_size)['rho_0x'])



plt.figure(13)
plt.plot( rho_0x0, -grad_x0*10**8, lw=2, c="c", label="w/(sqrt(2)rho) = 0.5")

plt.plot( rho_0x1, -grad_x1*10**8, lw=2, c="r", label="w/(sqrt(2)rho) = 0.65")

plt.plot( rho_0x2, -grad_x2*10**8, lw=2, c="g", label="w/(sqrt(2)rho) = 0.8")
         
plt.plot( rho_0x3, -grad_x3*10**8, lw=2, c="y", label="w/(sqrt(2)rho) = 1")

plt.plot( rho_0x4, -grad_x4*10**8, lw=2, c="m", label="w/(sqrt(2)rho) = 1.25")
         
plt.plot( rho_0x5, -grad_x5*10**8, lw=2, c="k", label="w/(sqrt(2)rho) = 1.5")

new_ticks1 = np.linspace(-4, 4, 9) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-6, 10, 9),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-4))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('grad_Fx(stiffness)10^(-8)',fontsize=20)

plt.title('rho = 60um, w0 = 10um',fontsize=20)
plt.grid()
plt.show()


MTP.table_parameter('rho/(w/sqrt(2)) = 1,1/2, 1/3', 'related to w', rho_0[1]* 10 ** 6, 'x-axis rho_0x/(w/sqrt(2)) = -4 to 4', '30')


