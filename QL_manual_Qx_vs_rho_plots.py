# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:53:15 2020

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
# integration method control: Qx vs rho plots
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
##################################
#12 plot Qx vs Rs for various offsets
##################################

a = 30 * 10 ** (-6)

w =  a * np.sqrt(2)

rho_0 = [0.5*a, 0 ]

rho_0[0] = [0.25* a, 0.5*a, 0.75 * a, a]

rho = np.linspace(10 ** (-6), 3*a, 100)

radial_flistx0 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )

radial_flistx1 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x1 = radial_flistx1  * c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x2 = radial_flistx2  * c / ( n_0 * P )

radial_flistx3 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0][3],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x3 = radial_flistx3  * c / ( n_0 * P )

plt.figure(12)
plt.plot( rho * np.sqrt(2)/w, Q_x0, lw=2, c="c", label="rho_0x = 7.5um")

plt.plot( rho * np.sqrt(2)/w, Q_x1, lw=2, c="r", label="rho_0x = 15um")
plt.plot( rho * np.sqrt(2)/w, Q_x2, lw=2, c="g", label="rho_0x = 22.5um")
plt.plot( rho * np.sqrt(2)/w, Q_x3, lw=2, c="y", label="rho_0x = 30um")

new_ticks1 = np.linspace(0, 3, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.5, 0.1, 7),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('sqrt(2)rho/w',fontsize=20)
plt.ylabel('Qx',fontsize=20)

#plt.title('Qx vs radius/beam waist w0 at w = 10w0 and rho_0x =10 * w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30', 'related to w', rho_0[1]* 10 ** 6, '7.5,15,22.5,30', 'x-axis 1 to 90')
'''



##################################
#12 plot Qx vs Rs for various w
##################################

a = 30 * 10 ** (-6)

w = [a* np.sqrt(2), 2*a, 2.5*a, 3*a]

rho_0 = [7.07* 10 ** (-6), 0 ]



rho = np.linspace(10 ** (-6), 3*a, 100)

radial_flistx0 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )

radial_flistx1 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x1 = radial_flistx1  * c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x2 = radial_flistx2  * c / ( n_0 * P )

radial_flistx3 = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[3], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x3 = radial_flistx3  * c / ( n_0 * P )

print (rho[np.argmax(Q_x0)])
print (rho[np.argmax(Q_x1)])
print (rho[np.argmax(Q_x2)])

plt.figure(12)
plt.plot( rho * 10 ** 6, Q_x0, lw=2, c="c", label="w = sqrt(2)*30um")

plt.plot( rho * 10 ** 6, Q_x1, lw=2, c="r", label="w = 60um")
plt.plot( rho * 10 ** 6, Q_x2, lw=2, c="g", label="w = 75um")
plt.plot( rho * 10 ** 6, Q_x3, lw=2, c="y", label="w = 90um")

new_ticks1 = np.linspace(0, 90, 10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.3, 0.1, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=15)

plt.xlabel('rho(um)',fontsize=20)
plt.ylabel('Qx',fontsize=20)

#plt.title('Qx vs radius/beam waist w0 at w = 10w0 and rho_0x =10 * w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30, 60, 75, 90', 'related to w', rho_0[1]* 10 ** 6, '15', 'x-axis 1 to 90')

